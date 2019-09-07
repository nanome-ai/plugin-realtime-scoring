import nanome
from nanome.util import Logs
from nanome.util.enums import NotificationTypes

import os
import subprocess
import tempfile
import itertools
import stat
from timeit import default_timer as timer

# TMP
from nanome._internal._structure._io._pdb.save import Options as PDBOptions
from nanome._internal._structure._io._sdf.save import Options as SDFOptions

SDF_OPTIONS = nanome.api.structure.Complex.io.SDFSaveOptions()
SDF_OPTIONS.write_bonds = True
PDB_OPTIONS = nanome.api.structure.Complex.io.PDBSaveOptions()
PDB_OPTIONS.write_bonds = True

BFACTOR_MIN = 0
BFACTOR_MAX = 50

BFACTOR_GAP = BFACTOR_MAX - BFACTOR_MIN
BFACTOR_MID = BFACTOR_GAP / 2
BFACTOR_HALF = BFACTOR_MAX - BFACTOR_MID

DIR = os.path.dirname(__file__)
RESULTS_PATH = os.path.join(DIR, 'dsx', 'results.txt')

class RealtimeScoring(nanome.PluginInstance):
    def benchmark_start(self, fn_name):
        self._benchmarks[fn_name] = timer()
    
    def benchmark_stop(self, fn_name):
        print('%s : %f' % (fn_name, timer() - self._benchmarks[fn_name]))

    def start(self):
        self._benchmarks = {}
        self._protein_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._site_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._ligands_converted = tempfile.NamedTemporaryFile(delete=False, suffix=".mol2")

        def button_pressed(button):
            if self._is_running:
                self.stop_scoring()
            else:
                self.start_scoring()

        menu = nanome.ui.Menu.io.from_json(os.path.join(DIR, 'menu.json'))
        self.menu = menu

        self._p_selection = menu.root.find_node("Selection Panel", True)
        self._p_results = menu.root.find_node("Results Panel", True)
        self._ls_receptors = menu.root.find_node("Receptor List", True).get_content()
        self._ls_ligands = menu.root.find_node("Ligands List", True).get_content()
        self._ls_results = self._p_results.get_content()
        self._btn_score = menu.root.find_node("Button", True).get_content()
        self._btn_score.register_pressed_callback(button_pressed)

        self._pfb_complex = nanome.ui.LayoutNode()
        self._pfb_complex.add_new_button()

        self._pfb_result = nanome.ui.LayoutNode()
        self._pfb_result.add_new_label()

        self._is_running = False
        self._obabel_running = False
        self._dsx_running = False
        self._ligands = None

        self._selected_receptor = None
        self._selected_ligands = []

        menu.enabled = True
        self.update_menu(menu)

        self.request_complex_list(self.update_lists)

    def on_run(self):
        self.menu.enabled = True
        self.update_menu(self.menu)

    def on_stop(self):
        os.remove(self._protein_input.name)
        os.remove(self._ligands_input.name)
        os.remove(self._site_input.name)
        os.remove(self._ligands_converted.name)

    def update(self):
        if not self._is_running:
            return

        if self._obabel_running: 
            if self._obabel_ligands_process.poll() is not None:
                self.benchmark_stop("obabel")
                self._obabel_running = False
                self.dsx_start()

        elif self._dsx_running:
            if self._dsx_process.poll() is not None:
                self._dsx_running = False
                self.get_updated_complexes()
                
                output, _ = self._dsx_process.communicate()
                self.parse_scores(output)
                self.display_results()
    
    def on_complex_added(self):
        self.request_complex_list(self.update_lists)

    def on_complex_removed(self):
        self.request_complex_list(self.update_lists)

    def start_scoring(self):
        if self._selected_receptor is None:
            self.send_notification(NotificationTypes.error, "Please select a receptor")
            return
            
        if len(self._selected_ligands) == 0:
            self.send_notification(NotificationTypes.error, "Please select at least one ligand")
            return
            
        self._is_running = True

        self._btn_score.set_all_text("Stop scoring")
        self._p_selection.enabled = False
        self._p_results.enabled = True
        self.update_menu(self.menu)
        
        self._color_stream = None
        self._scale_stream = None
        
        self.benchmark_start("scoring_loop")
        self.get_updated_complexes()

    def stop_scoring(self):
        self._is_running = False
        self._obabel_running = False
        self._dsx_running = False
        
        self._btn_score.set_all_text("Start scoring")
        self._p_selection.enabled = True
        self._p_results.enabled = False
        self.update_menu(self.menu)

        self.request_complex_list(self.update_lists)

    def get_updated_complexes(self):
        self.benchmark_stop("scoring_loop")
        print('*' * 32)
        self.benchmark_start("scoring_loop")

        index_list = [self._selected_receptor] + self._selected_ligands
        self.benchmark_start("get_full_complexes")
        self.request_complexes(index_list, self.setup_streams)

    def setup_streams(self, complex_list):
        self.benchmark_stop("get_full_complexes")
        if self._color_stream == None or self._scale_stream == None:
            indices = []
            for complex in complex_list[1:]:
                for atom in complex.atoms:
                    indices.append(atom.index)
                    atom.atom_mode = nanome.api.structure.Atom.AtomRenderingMode.Point
            def on_color_stream_ready(stream, error):
                self._color_stream = stream
                self.on_stream_ready(complex_list)
            def on_scale_stream_ready(stream, error):
                self._scale_stream = stream
                self.on_stream_ready(complex_list)
            def on_update_structure_done():
                self._struct_updated = True
                self.on_stream_ready(complex_list)
            self._struct_updated = False
            self.create_atom_stream(indices, nanome.api.streams.Stream.Type.color, on_color_stream_ready)
            self.create_atom_stream(indices, nanome.api.streams.Stream.Type.scale, on_scale_stream_ready)
            self.update_structures_deep(complex_list[1:], on_update_structure_done)
        else:
            self.obabel_start(complex_list)

    def on_stream_ready(self, complex_list):
        if self._color_stream != None and self._scale_stream != None and self._struct_updated == True:
            self.obabel_start(complex_list)

    def obabel_start(self, complex_list):
        receptor = complex_list[0]
        mat = receptor.get_complex_to_workspace_matrix()
        for atom in receptor.atoms:
            atom.position = mat * atom.position

        site = nanome.structure.Complex()
        site_molecule = nanome.structure.Molecule()
        site_chain = nanome.structure.Chain()
        site_residue = nanome.structure.Residue()
        site.add_molecule(site_molecule)
        site_molecule.add_chain(site_chain)
        site_chain.add_residue(site_residue)
        ligands = nanome.structure.Complex()

        for complex in complex_list[1:]:
            mat = complex.get_complex_to_workspace_matrix()
            for molecule in complex.molecules:
                index = molecule.index
                ligands.add_molecule(molecule)
                molecule.index = index
                for atom in molecule.atoms:
                    atom.position = mat * atom.position
                    index = atom.index
                    site_residue.add_atom(atom)
                    atom.index = index

        receptor.io.to_pdb(self._protein_input.name, PDB_OPTIONS)
        Logs.debug("Receptor:", self._protein_input.name)
        ligands.io.to_sdf(self._ligands_input.name, SDF_OPTIONS)
        Logs.debug("Ligands:", self._ligands_input.name)
        site.io.to_sdf(self._site_input.name, SDF_OPTIONS)

        obabel_ligands_args = ['obabel', '-isdf', self._ligands_input.name, '-omol2', '-O' + self._ligands_converted.name]

        try:
            self._obabel_running = True
            self._obabel_ligands_process = subprocess.Popen(obabel_ligands_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.benchmark_start("obabel")
        except:
            nanome.util.Logs.error("Couldn't execute obabel, please check if packet 'openbabel' is installed")
            return
        self._ligands = ligands
    
    def dsx_start(self):
        Logs.debug("Starting DSX")
        Logs.debug("Protein:", self._protein_input.name)
        Logs.debug("Ligand:", self._ligands_converted.name)

        dsx_path = os.path.join(DIR, 'dsx/dsx_linux_64.lnx')
        dsx_args = [dsx_path, '-P', self._protein_input.name, '-L', self._ligands_converted.name, '-D', 'pdb_pot_0511', '-pp', '-F', 'results.txt']
        try:
            self._dsx_process = subprocess.Popen(dsx_args, cwd=os.path.join(DIR, 'dsx'), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            self._dsx_running = True
        except:
            nanome.util.Logs.error("Couldn't execute dsx, please check if executable is in the plugin folder and has permissions. Try executing chmod +x " + dsx_path)
            return

    def parse_scores(self, output):
        lines = output.splitlines()
        count = len(lines)
        i = 0

        def find_next_ligand():
            nonlocal i
            while i < count:
                if lines[i].startswith("# Receptor-Ligand:"):
                    i += 1
                    return True
                i += 1
            return False
        
        if not find_next_ligand():
            Logs.error("Couldn't parse DSX scores")
            return

        ligand_index = 0
        has_next_ligand = True
        while has_next_ligand:
            scores = dict()
            last_tuple = None
            last_arr = None
            score_min = None
            score_max = None
            
            while i < count:
                line = lines[i]
                if line.startswith("# End of pair potentials"):
                    has_next_ligand = find_next_ligand()
                    break

                line_items = line.split("__")
                atom_items = line_items[1].split("_")
                score = float(line_items[2])
                tup = (int(atom_items[1]), int(atom_items[2]))
                
                if last_tuple != tup:    
                    if tup in scores:
                        last_arr = scores[tup]
                    else:
                        last_arr = []
                        scores[tup] = last_arr
                last_tuple = tup
                last_arr.append(score)

                if score_min == None or score < score_min:
                    score_min = score
                if score_max == None or score > score_max:
                    score_max = score
                i += 1

            if score_min != None and score_max != None:
                score_gap = max(score_max - score_min, 0.01)
                for atom, score_arr in scores.items():
                    score = sum(score_arr) / len(score_arr)
                    bfactor = ((score - score_min) / score_gap) * BFACTOR_GAP + BFACTOR_MIN
                    molecule = self._ligands._molecules[atom[0] - 1 + ligand_index]
                    atom_data = next(itertools.islice(molecule.atoms, atom[1] - 1, atom[1]))
                    atom_data._bfactor = bfactor
            else:
                for atom in self._ligands.atoms:
                    atom._bfactor = BFACTOR_MID

            ligand_index += 1

        colors = []
        scales = []
        for atom in self._ligands.atoms:
            colors.append(int(255 * ((atom._bfactor - BFACTOR_MIN) / BFACTOR_GAP)))
            colors.append(255 - int(255 * ((atom._bfactor - BFACTOR_MIN) / BFACTOR_GAP)))
            colors.append(0)
            scale = (abs(BFACTOR_MID - atom._bfactor) / BFACTOR_HALF) * 1.0 + 0.1
            scales.append(scale)

        self._color_stream.update(colors)
        self._scale_stream.update(scales)

    def display_results(self):
        scores = []
        with open(RESULTS_PATH) as results_file:
            results = results_file.readlines()
            count = len(results)
            i = results.index('@RESULTS\n') + 4

            while i < count:
                if results[i] == '\n':
                    break
                result = results[i].split('|')
                name = result[1].strip()
                score = result[5].strip()
                scores.append('%s - %s' % (name, score))
                i += 1
        
        self._ls_results.items = []
        for score in scores:
            clone = self._pfb_result.clone()
            lbl = clone.get_content()
            lbl.text_value = score
            self._ls_results.items.append(clone)
        self.update_content(self._ls_results)

    def update_lists(self, complex_list):
        if self._is_running:
            return

        def update_selected_ligands():
            self._selected_ligands = []
            for item in self._ls_ligands.items:
                btn = item.get_content()
                if btn.selected:
                    self._selected_ligands.append(btn.index)

        def receptor_pressed(receptor):
            for item in self._ls_receptors.items:
                item.get_content().selected = False

            receptor.selected = True
            self._selected_receptor = receptor.index

            for item in self._ls_ligands.items:
                ligand = item.get_content()
                ligand.unusable = receptor.index == ligand.index
                if ligand.selected and ligand.unusable:
                    ligand.selected = False
            update_selected_ligands()
            
            self.update_menu(self.menu)

        def ligand_pressed(ligand):
            ligand.selected = not ligand.selected
            self.update_content(ligand)
            update_selected_ligands()
        
        def populate_list(ls, cb):
            ls.items = []
            for complex in complex_list:
                clone = self._pfb_complex.clone()
                btn = clone.get_content()
                btn.set_all_text(complex.name)
                btn.index = complex.index
                btn.register_pressed_callback(cb)
                ls.items.append(clone)
            self.update_content(ls)

        self._selected_receptor = None
        self._selected_ligands = []
        populate_list(self._ls_receptors, receptor_pressed)
        populate_list(self._ls_ligands, ligand_pressed)

def main():
    plugin = nanome.Plugin("RTS II", "(Development Version) Display realtime scoring information about a selected ligand", "Docking", False)
    plugin.set_plugin_class(RealtimeScoring)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()