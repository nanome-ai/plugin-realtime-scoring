import nanome
from nanome.util import Logs

import os
import subprocess
import tempfile
import itertools

# TMP
from nanome._internal._structure._io._pdb.save import Options as PDBOptions
from nanome._internal._structure._io._sdf.save import Options as SDFOptions

SDF_OPTIONS = nanome.api.structure.Complex.io.SDFSaveOptions()
SDF_OPTIONS.write_bonds = True
PDB_OPTIONS = nanome.api.structure.Complex.io.PDBSaveOptions()
PDB_OPTIONS.write_bonds = True

MIN_BFACTOR = 0
MAX_BFACTOR = 50

bfactor_gap = MAX_BFACTOR - MIN_BFACTOR


class RealtimeScoring(nanome.PluginInstance):
    def start(self):
        def button_pressed(button):
            if self._is_running:
                self.stop_scoring()
            else:
                if self._selected_receptor == None:
                    self.send_notification(nanome.util.enums.NotificationTypes.error, "Please select a receptor")
                    return
                self.start_scoring()

        menu = nanome.ui.Menu.io.from_json(os.path.join(os.path.dirname(__file__), '_scoring_menu.json'))
        self.menu = menu

        self._list = menu.root.find_node("List", True).get_content()
        self._button = menu.root.find_node("Button", True).get_content()
        self._button.register_pressed_callback(button_pressed)

        self._complex_item_prefab = nanome.ui.LayoutNode()
        child = self._complex_item_prefab.create_child_node()
        child.add_new_button()

        self._score_item_prefab = nanome.ui.LayoutNode()
        child = self._score_item_prefab.create_child_node()
        child.add_new_label()

        self._is_running = False
        self._request_workspace = False
        self._smina_running = False
        self._dsx_running = False
        self.request_complex_list(self.on_complex_list_received)

        menu.enabled = True
        self._menu = menu
        self.update_menu(menu)

    def on_run(self):
        self._menu.enabled = True
        self.update_menu(self._menu)

    def update(self):
        if self._is_running == False:
            return

        if self._request_workspace:
            self._request_workspace = False
            self.request_complex_list(self.on_complex_list_received_scoring)
        elif self._smina_running:
            if self._smina_process.poll() is not None and self._obabel_protein_process.poll() is not None and self._obabel_ligands_process.poll() is not None:
                self.scoring_done()
        elif self._dsx_running:
            if self.check_dsx():
                self.dsx_done()

    def check_dsx(self):
        poll_result = self._dsx_process.poll()
        if poll_result is None:
            output, _ = self._dsx_process.communicate()
            self.parse_scores(output)
        return poll_result is not None

    def parse_scores(self, output):
        lines = output.splitlines()
        count = len(lines)
        i = 0
    
        while i < count:
            if lines[i].startswith("# Receptor-Ligand:"):
                break
            i += 1
        if i >= count:
            Logs.error("Couldn't parse DSX scores")
            return
        i += 1

        scores = dict()
        last_tuple = None
        last_arr = None
        score_min = None
        score_max = None
        while i < count:
            line = lines[i]
            if line.startswith("# End of pair potentials"):
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

        score_gap = score_max - score_min
        for atom, score_arr in scores.items():
            score = sum(score_arr) / len(score_arr)
            bfactor = ((score - score_min) / score_gap) * bfactor_gap + MIN_BFACTOR
            molecule = self._ligands._molecules[atom[0] - 1]
            atom = next(itertools.islice(molecule.atoms, atom[1] - 1, atom[1]))
            atom._bfactor = bfactor
        
        for complex in self._to_update:
            mat = complex.transform.get_workspace_to_complex_matrix()
            for atom in complex.atoms:
                atom.molecular.position = mat * atom.molecular.position

        self.update_structures_deep(self._to_update, self._update_done)

    def on_complex_added(self):
        self.request_complex_list(self.on_complex_list_received)

    def on_complex_removed(self):
        self.request_complex_list(self.on_complex_list_received)

    def start_scoring(self):
        self._is_running = True
        Logs.debug("Start scoring")
        self._button.set_all_text("Stop scoring")
        self.update_content(self._button)
        self._request_workspace = True

    def stop_scoring(self):
        self._is_running = False
        Logs.debug("Stop scoring")
        self._button.set_all_text("Start scoring")
        self.update_content(self._button)
        self.request_complex_list(self.on_complex_list_received)

    def on_complex_list_received_scoring(self, complex_list):
        index_list = [self._selected_receptor.complex.index]
        for complex in complex_list:
            if complex.index == self._selected_receptor.complex.index:
                continue
            # TODO: only add selected complexes
            index_list.append(complex.index)
        self.request_complexes(index_list, self.on_full_complexes_received)

    def scoring_done(self):
        docked_ligands = nanome.structure.Complex.io.from_sdf(path=self._ligand_output.name)

        Logs.debug("Scoring done")
        self._list.items = []
        for molecule in docked_ligands.molecules:
            clone = self._score_item_prefab.clone()
            lbl = clone.get_children()[0].get_content()
            lbl.text_value = molecule.molecular.name + " - " + molecule._associated["> <minimizedAffinity>"]
            self._list.items.append(clone)

        self.update_menu(self._menu)
        Logs.debug("Starting DSX")
        Logs.debug("Protein:", self._protein_converted.name)
        Logs.debug("Ligand:", self._ligands_converted.name)

        dsx_path = os.path.join(os.path.dirname(__file__), 'dsx/dsx_linux_64.lnx')
        dsx_args = [dsx_path, '-P', self._protein_converted.name, '-L', self._ligands_converted.name, '-D', 'pdb_pot_0511', '-pp']
        try:
            self._dsx_process = subprocess.Popen(dsx_args, cwd=os.path.join(os.path.dirname(__file__), 'dsx'), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        except:
            nanome.util.Logs.error("Couldn't execute dsx, please check if executable is in the plugin folder and has permissions. Try executing chmod +x " + dsx_path)
            return

        self._smina_running = False
        self._dsx_running = True
    
    def dsx_done(self):
        os.remove(self._protein_input.name)
        os.remove(self._ligands_input.name)
        os.remove(self._site_input.name)
        os.remove(self._ligand_output.name)
        os.remove(self._protein_converted.name)
        os.remove(self._ligands_converted.name)

        self._dsx_running = False

    def _update_done(self):
        self._request_workspace = True

    def on_full_complexes_received(self, complex_list):
        receptor = complex_list[0]
        mat = receptor.transform.get_complex_to_workspace_matrix()
        for atom in receptor.atoms:
            atom.molecular.position = mat * atom.molecular.position

        site = nanome.structure.Complex()
        site_molecule = nanome.structure.Molecule()
        site_chain = nanome.structure.Chain()
        site_residue = nanome.structure.Residue()
        site.add_molecule(site_molecule)
        site_molecule.add_chain(site_chain)
        site_chain.add_residue(site_residue)
        ligands = nanome.structure.Complex()

        self._to_update = complex_list[1:]

        for complex in self._to_update:
            mat = complex.transform.get_complex_to_workspace_matrix()
            for molecule in complex.molecules:
                index = molecule.index
                ligands.add_molecule(molecule)
                molecule.index = index
                for atom in molecule.atoms:
                    atom.molecular.position = mat * atom.molecular.position
                    index = atom.index
                    site_residue.add_atom(atom)
                    atom.index = index

        self._protein_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._site_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._ligand_output = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._protein_converted = tempfile.NamedTemporaryFile(delete=False, suffix=".mol2")
        self._ligands_converted = tempfile.NamedTemporaryFile(delete=False, suffix=".mol2")
        receptor.io.to_pdb(self._protein_input.name, PDB_OPTIONS)
        Logs.debug("Receptor:", self._protein_input.name)
        ligands.io.to_sdf(self._ligands_input.name, SDF_OPTIONS)
        Logs.debug("Ligands:", self._ligands_input.name)
        site.io.to_sdf(self._site_input.name, SDF_OPTIONS)

        smina_path = os.path.join(os.path.dirname(__file__), 'smina')
        smina_args = [smina_path, '--autobox_ligand', self._site_input.name, '--score_only', '-r', self._protein_input.name, '--ligand', self._ligands_input.name, '--out', self._ligand_output.name]
        obabel_protein_args = ['obabel', '-ipdb', self._protein_input.name, '-omol2', '-O' + self._protein_converted.name]
        obabel_ligands_args = ['obabel', '-isdf', self._ligands_input.name, '-omol2', '-O' + self._ligands_converted.name]

        try:
            self._smina_process = subprocess.Popen(smina_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except:
            nanome.util.Logs.error("Couldn't execute smina, please check if executable is in the plugin folder and has permissions. Try executing chmod +x " + smina_path)
            return
        try:
            self._obabel_protein_process = subprocess.Popen(obabel_protein_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self._obabel_ligands_process = subprocess.Popen(obabel_ligands_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except:
            nanome.util.Logs.error("Couldn't execute obabel, please check if packet 'openbabel' is installed")
            return
        self._smina_running = True
        self._ligands = ligands

    def on_complex_list_received(self, complex_list):
        if self._is_running:
            return
        
        def complex_pressed(button):
            lastSelected = self._selected_receptor
            if lastSelected != None:
                lastSelected.selected = False
                self.update_content(lastSelected)
            button.selected = True
            self._selected_receptor = button
            self.update_content(button)

        self._selected_receptor = None
        self._list.items = []

        for complex in complex_list:
            clone = self._complex_item_prefab.clone()
            ln_btn = clone.get_children()[0]
            btn = ln_btn.get_content()
            btn.set_all_text(complex.molecular.name)
            btn.complex = complex
            btn.register_pressed_callback(complex_pressed)
            self._list.items.append(clone)

        self.update_menu(self._menu)

def main():
    plugin = nanome.Plugin("Realtime Scoring", "Display realtime scoring information about a selected ligand", "Docking", False)
    plugin.set_plugin_class(RealtimeScoring)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()
