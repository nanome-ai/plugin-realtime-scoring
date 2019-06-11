import nanome
from nanome.util import Logs

import os
import subprocess
import tempfile

# TMP
from nanome._internal._structure._io._pdb.save import Options as PDBOptions
from nanome._internal._structure._io._sdf.save import Options as SDFOptions

SDF_OPTIONS = nanome.api.structure.Complex.io.SDFSaveOptions()
SDF_OPTIONS.write_bonds = True
PDB_OPTIONS = nanome.api.structure.Complex.io.PDBSaveOptions()
PDB_OPTIONS.write_bonds = True


class RealtimeScoring(nanome.PluginInstance):
    def start(self):
        def button_pressed(button):
            if self._is_running:
                self.stop_scoring()
            else:
                if self._selected_receptor == None:
                    self.send_notification(nanome.util.NotificationTypes.error, "Please select a receptor")
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
            Logs.message("Output:", output.decode('utf-8').splitlines())
        return poll_result is not None

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

        self._list.items = []
        for molecule in docked_ligands.molecules:
            clone = self._score_item_prefab.clone()
            lbl = clone.get_children()[0].get_content()
            lbl.text_value = molecule.molecular.name + " - " + molecule.molecular._associated["> <minimizedAffinity>"]
            self._list.items.append(clone)

        self.update_menu(self._menu)

        dsx_path = os.path.join(os.path.dirname(__file__), 'dsx/dsx_linux_64.lnx')
        dsx_args = [dsx_path, '-P', self._protein_converted.name, '-L', self._ligands_converted.name, '-D', 'pdb_pot_0511', '-pp']
        try:
            self._dsx_process = subprocess.Popen(dsx_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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

        for complex in complex_list[1:]:
            mat = complex.transform.get_complex_to_workspace_matrix()
            for molecule in complex.molecules:
                ligands.add_molecule(molecule)
                for atom in molecule.atoms:
                    atom.molecular.position = mat * atom.molecular.position
                    site_residue.add_atom(atom)

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
