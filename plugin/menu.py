import nanome
from os import path
from nanome.api import ui, shapes
from nanome.util import Logs, async_callback
from nanome.util.enums import NotificationTypes

BASE_PATH = path.dirname(f'{path.realpath(__file__)}')
MENU_PATH = path.join(BASE_PATH, 'menu_json', 'menu.json')


class MainMenu:

    def __init__(self, plugin_instance):
        self._menu = ui.Menu.io.from_json(MENU_PATH)
        self.plugin = plugin_instance
        self.ln_selection: ui.LayoutNode = self._menu.root.find_node("Selection Panel", True)
        self.ln_selection.enabled = True
        self.ln_results: ui.LayoutNode = self._menu.root.find_node("Results Panel", True)
        self._ls_receptors: ui.UIList = self._menu.root.find_node("Receptor List", True).get_content()
        self._ls_ligands: ui.UIList = self._menu.root.find_node("Ligands List", True).get_content()
        self.btn_score: ui.Button = self._menu.root.find_node("btn_score", True).get_content()
        self.btn_score.register_pressed_callback(self.on_scoring_button_pressed)
        self.btn_score.toggle_on_press = True

        self._pfb_complex = nanome.ui.LayoutNode()
        pfb_btn = self._pfb_complex.add_new_button()
        pfb_btn.toggle_on_press = True
        self._pfb_result = nanome.ui.LayoutNode()
        self._pfb_result.add_new_label()

    @async_callback
    async def on_scoring_button_pressed(self, button):
        if button.selected or not self.plugin.realtime_enabled:
            await self.start_scoring()
        else:
            self.stop_scoring()
        self.plugin.update_content(button)

    async def start_scoring(self):
        receptor_index = self.receptor_index
        residue_indices = self.ligand_residue_indices
        Logs.message("Start Scoring")
        Logs.debug(f"Residue Count: {len(residue_indices)}")
        if self.plugin.realtime_enabled:
            # Don't switch panels if realtime is enabled
            self.ln_selection.enabled = False
            self.ln_results.enabled = True

        if receptor_index is None:
            self.plugin.send_notification(
                NotificationTypes.error, "Please select a receptor")
            return

        if len(residue_indices) == 0:
            self.plugin.send_notification(
                NotificationTypes.error, "Please select at least one ligand")
            return

        # Add loading message to results panel
        results_list = self.ln_results.get_content()
        results_list.items = []
        clone = self._pfb_result.clone()
        lbl = clone._get_content()
        lbl.text_value = 'Loading...'
        results_list.items.append(clone)
        if self.plugin.realtime_enabled:
            self.plugin.update_menu(self._menu)

        await self.plugin.setup_receptor_and_ligands(
            receptor_index, residue_indices)
        await self.plugin.score_ligands()
        if self.plugin.realtime_enabled:
            self._menu.title = "Scores"

    def stop_scoring(self):
        Logs.message("Stopping Scoring Streams")
        self.plugin.stop_scoring()
        self.ln_selection.enabled = True
        self.ln_results.enabled = False
        self.plugin.update_menu(self._menu)

    @async_callback
    async def render(self, force_enable=False):
        if not self.plugin.realtime_enabled:
            # Change behavior of button if realtime is disabled
            self.btn_score.toggle_on_press = False
            self.btn_score.text.value.set_all("Start Scoring")
            self.btn_score.text.value.unusable = "Scoring..."
            self.btn_score.toggle_on_press = False
            self.btn_score.disable_on_press = True
            self.plugin.update_content(self.btn_score)
        complex_list = self.plugin.complex_cache
        self.populate_list(self._ls_receptors, complex_list, self.on_receptor_pressed)
        self.populate_list(self._ls_ligands, complex_list)
        if force_enable:
            self._menu.enabled = True
        self.plugin.update_menu(self._menu)

    @property
    def receptor_index(self):
        """Get the complex index from the currently selected button."""
        for item in self._ls_receptors.items:
            if item.get_content().selected:
                return item.get_content().index

    @property
    def ligand_residue_indices(self):
        """Get the list of residues from the currently selected buttons."""
        residues = []
        for item in self._ls_ligands.items:
            btn = item.get_content()
            if btn.selected and hasattr(btn, 'residue_indices'):
                residues.extend(btn.residue_indices)
        return residues

    @async_callback
    async def on_receptor_pressed(self, receptor_btn):
        # Deselect every other receptor button
        for item in self._ls_receptors.items:
            item_btn = item.get_content()
            item_btn.selected = \
                item_btn._content_id == receptor_btn._content_id

        # Remove previously extracted ligand items from list
        # Iterate in reverse order for simpler deletions
        for i in range(len(self._ls_ligands.items) - 1, -1, -1):
            item = self._ls_ligands.items[i]
            item_btn = item.get_content()
            if hasattr(item_btn, 'extracted_ligand'):
                self._ls_ligands.items.remove(item)

        # Extract ligands from receptor, and add as entry to ligand list
        receptor_index = receptor_btn.index
        complex_list = self.plugin.complex_cache
        receptor = next(
            cmp for cmp in complex_list
            if cmp.index == receptor_index)
        mol = next(
            ml for i, ml in enumerate(receptor.molecules)
            if i == receptor.current_frame)
        ligands = await mol.get_ligands()
        # Create new button for every ligand.
        for lig in ligands:
            clone = self._pfb_complex.clone()
            btn = clone.get_content()
            btn.text.value.set_all(lig.name)
            btn.index = receptor.index
            btn.residue_indices = [res.index for res in lig.residues]
            btn.extracted_ligand = True
            # make sure structure tree is stored on residue, we will need it later
            for residue in lig.residues:
                # Find the chain that this residue belongs to, and set parent
                res_chain = next(
                    chain for chain in mol.chains
                    if chain.name == residue.chain.name
                )
                residue._parent = res_chain
            self._ls_ligands.items.append(clone)

        # Disable corresponding ligand button for selected receptor
        for item in self._ls_ligands.items:
            ligand_btn = item.get_content()
            receptor_btn_txt = receptor_btn.text.value.selected
            ligand_btn_txt = ligand_btn.text.value.selected
            ligand_btn.unusable = receptor_btn_txt == ligand_btn_txt
        self.plugin.update_menu(self._menu)

    def populate_list(self, ui_list, complex_list, callback=None):
        ui_list.items = []
        for comp in complex_list:
            clone = self._pfb_complex.clone()
            btn = clone.get_content()
            btn.text.value.set_all(comp.full_name)
            btn.index = comp.index
            btn.residue_indices = [res.index for res in comp.residues]
            if callback:
                btn.register_pressed_callback(callback)
            ui_list.items.append(clone)
        self.plugin.update_content(ui_list)

    def update_ligand_scores(self, aggregate_score_list):
        results_list = self.ln_results.get_content()
        results_list.items = []
        scores_set = aggregate_score_list[0]
        if not scores_set:
            Logs.warning("No aggregate scores returned by scoring algorithm.")
            return
        scores = scores_set[0]
        for name, score in scores.items():
            clone = self._pfb_result.clone()
            lbl = clone._get_content()
            lbl.text_value = '{}: {}'.format(name, score)
            results_list.items.append(clone)
        self.plugin.update_content(results_list)
