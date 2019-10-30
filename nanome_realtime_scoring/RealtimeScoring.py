import nanome
from nanome.util import Logs
from nanome.util.enums import NotificationTypes

import os
import functools
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
        if not fn_name in self._benchmarks:
            self._benchmarks[fn_name] = [0, 0, 0]
        self._benchmarks[fn_name][0] = timer()
    
    def benchmark_stop(self, fn_name):
        entry = self._benchmarks[fn_name]
        time = timer() - entry[0]
        entry[1] += time
        entry[2] += 1
        avg = entry[1] / entry[2]

        nanome.util.Logs.debug('{:>10}    {:.2f}    (avg {:.2f})'.format(fn_name, time, avg))

    def start(self):
        self._benchmarks = {}

        self._protein_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
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

        self._receptor_index = None
        self._ligand_indices = []
        self._complexes = []

        menu.enabled = True
        self.update_menu(menu)

        self.request_complex_list(self.update_lists)

    def on_run(self):
        self.menu.enabled = True
        self.update_menu(self.menu)

    def on_stop(self):
        os.remove(self._protein_input.name)
        os.remove(self._ligands_input.name)
        os.remove(self._ligands_converted.name)

    def dsx_start(self):
        dsx_path = os.path.join(DIR, 'dsx', 'dsx_linux_64.lnx')
        dsx_args = [dsx_path, '-P', self._protein_input.name, '-L', self._ligands_converted.name, '-D', 'pdb_pot_0511', '-pp', '-F', 'results.txt', '-v']
        try:
            self._dsx_process = subprocess.Popen(dsx_args, cwd=os.path.join(DIR, 'dsx'), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            self._dsx_running = True
            nanome.util.Logs.debug("dsx running: " + str(self._dsx_running))
            self.benchmark_start("dsx")
        except:
            nanome.util.Logs.error("Couldn't execute dsx, please check if executable is in the plugin folder and has permissions. Try executing chmod +x " + dsx_path)
            return

    def update(self):
        if not self._is_running:
            return

        if self._obabel_running: 
            if self._obabel_process.poll() is not None:
                self.benchmark_stop("obabel")
                self._obabel_running = False
                self.dsx_start()
        elif self._dsx_running:
            if self._dsx_process.poll() is not None:
                self.benchmark_stop("dsx")
                self._dsx_running = False
                
                self.get_updated_complexes()
                dsx_output, _ = self._dsx_process.communicate()
                self.parse_scores(dsx_output)
                self.display_results()
    
    def on_complex_added(self):
        self.request_complex_list(self.update_lists)

    def on_complex_removed(self):
        self.request_complex_list(self.update_lists)

    def start_scoring(self):
        if self._receptor_index is None:
            self.send_notification(NotificationTypes.error, "Please select a receptor")
            return
            
        if len(self._ligand_indices) == 0:
            self.send_notification(NotificationTypes.error, "Please select at least one ligand")
            return
            
        self._is_running = True
        self._obabel_running = False
        self._dsx_running = False

        self._btn_score.set_all_text("Stop scoring")
        self._p_selection.enabled = False
        self._p_results.enabled = True
        self.update_menu(self.menu)
        
        self._color_stream = None
        self._scale_stream = None
        
        self.benchmark_start("total")
        self.get_full_complexes()

    def stop_scoring(self):
        self._is_running = False
        self._obabel_running = False
        self._dsx_running = False
        
        self._btn_score.set_all_text("Start scoring")
        self._p_selection.enabled = True
        self._p_results.enabled = False
        self.update_menu(self.menu)

        # for complex in self._complexes[1:]:
        #     for atom in complex.atoms:
        #         atom.atom_mode = atom._old_atom_mode
        # self.update_structures_deep(self._complexes[1:])

        self.request_complex_list(self.update_lists)

    def get_full_complexes(self):
        def set_complexes(complex_list):
            self._complexes = complex_list

            for complex_i in range(0, len(complex_list)):
                complex = complex_list[complex_i].convert_to_frames()
                complex.index = complex_list[complex_i].index
                complex_list[complex_i] = complex
                for atom in complex.atoms:
                    atom._old_position = atom.position
                    if complex_i > 0:
                        atom.labeled = True
                        if " ({})".format(atom.symbol) not in atom.label_text:
                            atom.label_text += " ({})".format(atom.symbol)
            self.update_structures_deep(complex_list[1:], functools.partial(self.request_complexes, index_list, self.setup_streams))
        
        index_list = [self._receptor_index] + self._ligand_indices
        self.request_complexes(index_list, set_complexes)

    def get_updated_complexes(self):
        self.benchmark_stop("total")
        nanome.util.Logs.debug('*' * 32)
        self.benchmark_start("total")
    
        def update_complexes(complex_list):
            # update self._complexes positions from shallow list
            for complex in self._complexes:
                for thing in complex_list:
                    if thing.index == complex.index:
                        complex.position = thing.position
                        complex.rotation = thing.rotation
                        break

            self.benchmark_stop("update")
            self.prepare_complexes(self._complexes)
        
        self.benchmark_start("update")
        self.request_complex_list(update_complexes)

    def setup_streams(self, complex_list):
        if self._color_stream == None or self._scale_stream == None:
            indices = []
            for complex in complex_list[1:]:
                for atom in complex.atoms:
                    indices.append(atom.index)
                    # atom._old_atom_mode = atom.atom_mode
                    atom.atom_mode = nanome.api.structure.Atom.AtomRenderingMode.Point
            def on_stream_ready(complex_list):
                if self._color_stream != None and self._scale_stream != None and self._struct_updated == True:
                    self.get_updated_complexes()
            def on_color_stream_ready(stream, error):
                self._color_stream = stream
                on_stream_ready(complex_list)
            def on_scale_stream_ready(stream, error):
                self._scale_stream = stream
                on_stream_ready(complex_list)
            def on_update_structure_done():
                self._struct_updated = True
                on_stream_ready(complex_list)
            self._struct_updated = False
            self.create_atom_stream(indices, nanome.api.streams.Stream.Type.color, on_color_stream_ready)
            self.create_atom_stream(indices, nanome.api.streams.Stream.Type.scale, on_scale_stream_ready)
            self.update_structures_deep(complex_list[1:], on_update_structure_done)
        else:
            self.get_updated_complexes()

    def prepare_complexes(self, complex_list):
        receptor = complex_list[0]
        for complex in complex_list:
            mat = complex.get_complex_to_workspace_matrix()
            for atom in complex.atoms:
                atom.position = mat * atom._old_position

        ligands = nanome.structure.Complex()
        self._ligands = ligands

        for complex in complex_list[1:]:
            complex = complex.convert_to_frames()
            for molecule in complex.molecules:
                index = molecule.index
                ligands.add_molecule(molecule)
                molecule.index = index
        
        receptor.io.to_pdb(self._protein_input.name, PDB_OPTIONS)
        ligands.io.to_sdf(self._ligands_input.name, SDF_OPTIONS)

        self.obabel_start()

    def obabel_start(self):
        obabel_args = ['obabel', '-isdf', self._ligands_input.name, '-omol2', '-O' + self._ligands_converted.name]
        try:
            self._obabel_process = subprocess.Popen(obabel_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self._obabel_running = True
            self.benchmark_start("obabel")
        except:
            nanome.util.Logs.error("Couldn't execute obabel, please check if packet 'openbabel' is installed")
            return

    def parse_scores(self, dsx_output):
        with open(os.path.join(DIR, 'dsx_output.txt'), 'w') as dsx_score_file:
            dsx_score_file.write(dsx_output)

        lines = dsx_output.splitlines()
        number_of_lines = len(lines)
        
        line_index = 0
        def find_next_ligand():
            nonlocal line_index
            while line_index < number_of_lines:
                if lines[line_index].startswith("# Receptor-Ligand:"):
                    line_index += 1
                    return True
                line_index += 1
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
            
            while line_index < number_of_lines:
                line = lines[line_index]
                if line.startswith("# End of pair potentials"):
                    has_next_ligand = find_next_ligand()
                    break
                line_items = line.split("__")
                atom_items = line_items[1].split("_")
                score = float(line_items[2])
                ireceptor_iligand = (int(atom_items[1]), int(atom_items[2]))
                
                if last_tuple != ireceptor_iligand:    
                    if ireceptor_iligand in scores:
                        last_arr = scores[ireceptor_iligand]
                    else:
                        last_arr = []
                        scores[ireceptor_iligand] = last_arr
                last_tuple = ireceptor_iligand
                last_arr.append(score)

                if score_min == None or score < score_min:
                    score_min = score
                if score_max == None or score > score_max:
                    score_max = score
                line_index += 1
            
            if score_max is None or score_min is None:
                continue
            score_gap = max(score_max - score_min, 0.01)
            try:
                for atom_tuple, score_arr in scores.items():
                    number_of_lines += 1
                    score = sum(score_arr) / len(score_arr)
                    bfactor = ((score - score_min) / score_gap) * BFACTOR_GAP + BFACTOR_MIN
                    bfactor_two = score / (-score_min if score < 0 else score_max)
                    molecule = self._ligands._molecules[atom_tuple[0] - 1 + ligand_index]
                    if not hasattr(molecule, "atom_score_limits"):
                        molecule.atom_score_limits = [float('inf'), float('-inf')]
                    if score < molecule.atom_score_limits[0]:
                        molecule.atom_score_limits[0] = score
                    elif score > molecule.atom_score_limits[1]:
                        molecule.atom_score_limits[1] = score

                    atom = next(itertools.islice(molecule.atoms, atom_tuple[1] - 1, atom_tuple[1]))
                    atom.score = score
            except:
                err_msg = "Error parsing ligand scores. Are your ligands missing bonds?"
                self._send_notification(NotificationTypes.error, err_msg)
                nanome.util.Logs.error(err_msg)
                self.running = False
                return

            ligand_index += 1

        colors = []
        scales = []
        for atom in self._ligands.atoms:
            red = 0
            green = 0
            blue = 0
            if hasattr(atom, "score"):
                ligand = atom.molecule
                denominator = -ligand.atom_score_limits[0] if atom.score < 0 else ligand.atom_score_limits[1]
                norm_score = atom.score / denominator
                red = 255 if norm_score > 0 else 0
                green = 255 if -norm_score >= 0 else 0
                red_scale = int(max(norm_score*255, 0))
                green_scale = int(max(-norm_score*255, 0))
                scale = max(green_scale,red_scale) / 255. * 1.5
            else:
                green = 0
                red = 255
                scale = 1
                
            colors.append(red)
            colors.append(green)
            colors.append(blue)
            scales.append(scale)
        
        try:
            self._color_stream.update(colors)
            self._scale_stream.update(scales)
        except:
            print(colors)

    def display_results(self):
        scores = []
        with open(RESULTS_PATH) as results_file:
            results = results_file.readlines()
            number_of_lines = len(results)
            line_index = results.index('@RESULTS\n') + 4

            while line_index < number_of_lines:
                if results[line_index] == '\n':
                    break
                result = results[line_index].split('|')
                name = result[1].strip()
                score = result[5].strip()
                scores.append('%s: %s' % (name, score))
                line_index += 1
        
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
            self._ligand_indices = []
            for item in self._ls_ligands.items:
                btn = item.get_content()
                if btn.selected:
                    self._ligand_indices.append(btn.index)

        def receptor_pressed(receptor):
            for item in self._ls_receptors.items:
                item.get_content().selected = False

            receptor.selected = True
            self._receptor_index = receptor.index

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

        self._receptor_index = None
        self._ligand_indices = []
        populate_list(self._ls_receptors, receptor_pressed)
        populate_list(self._ls_ligands, ligand_pressed)

def main():
    plugin = nanome.Plugin("Realtime Scoring", "Display realtime scoring information about a selected ligand", "Docking", False)
    plugin.set_plugin_class(RealtimeScoring)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()
