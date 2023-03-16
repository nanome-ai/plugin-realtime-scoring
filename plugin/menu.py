import nanome
from os import path
from nanome.api import ui
from nanome.util import Logs, async_callback
from nanome.util.enums import NotificationTypes

BASE_PATH = path.dirname(f'{path.realpath(__file__)}')
MENU_PATH = path.join(BASE_PATH, 'menu_json', 'menu.json')


class MainMenu:

    def __init__(self, plugin_instance):
        super().__init__()
        self._menu = ui.Menu.io.from_json(MENU_PATH)
        self.plugin = plugin_instance
        self.ln_selection: ui.LayoutNode = self._menu.root.find_node("Selection Panel", True)
        self.ln_selection.enabled = True
        self.ln_results: ui.LayoutNode = self._menu.root.find_node("Results Panel", True)
        self._ls_receptors: ui.UIList = self._menu.root.find_node("Receptor List", True).get_content()
        self._ls_ligands: ui.UIList = self._menu.root.find_node("Ligands List", True).get_content()
        self._ls_results: ui.UIList = self.ln_results.get_content()
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
        receptor_index = self.receptor_index
        ligand_indices = self.ligand_indices
        if button.selected or not self.plugin.realtime_enabled:
            await self.start_scoring(receptor_index, ligand_indices)
        else:
            self.stop_scoring()
        self.plugin.update_content(button)

    async def start_scoring(self, receptor_index, ligand_indices):
        Logs.message("Start Scoring")
        if self.plugin.realtime_enabled:
            # Don't switch panels if realtime is enabled
            self.ln_selection.enabled = False
            self.ln_results.enabled = True

        if receptor_index is None:
            self.send_notification(NotificationTypes.error, "Please select a receptor")
            return

        if len(ligand_indices) == 0:
            self.send_notification(NotificationTypes.error, "Please select at least one ligand")
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
        await self.plugin.setup_receptor_and_ligands(receptor_index, ligand_indices)
        await self.plugin.score_ligands()
        if self.plugin.realtime_enabled:
            self._menu.title = "Scores"

    def stop_scoring(self):
        self.plugin.stop_streams()
        self.ln_selection.enabled = True
        self.ln_results.enabled = False
        self.plugin.update_menu(self._menu)

    @async_callback
    async def render(self, complex_list, force_enable=False):
        if not self.plugin.realtime_enabled:
            # Change behavior of button if realtime is disabled
            self.btn_score.toggle_on_press = False
            self.btn_score.text.value.set_all("Start Scoring")
            self.btn_score.text.value.unusable = "Scoring..."
            self.btn_score.toggle_on_press = False
            self.btn_score.disable_on_press = True
            self.plugin.update_content(self.btn_score)
        self.populate_list(self._ls_receptors, complex_list, self.on_receptor_pressed)
        self.populate_list(self._ls_ligands, complex_list)
        self.complex_list = complex_list
        if force_enable:
            self._menu.enabled = True
        self.plugin.update_menu(self._menu)

    @property
    def receptor_index(self):
        for item in self._ls_receptors.items:
            if item.get_content().selected:
                return item.get_content().index

    @property
    def ligand_indices(self):
        indices = []
        for item in self._ls_ligands.items:
            btn = item.get_content()
            if btn.selected:
                indices.append(btn.index)
        return indices

    def on_receptor_pressed(self, btn):
        for item in self._ls_receptors.items:
            item.get_content().selected == False
        btn.selected = True

        for item in self._ls_ligands.items:
            ligand_btn = item.get_content()
            ligand_btn.unusable = btn.index == ligand_btn.index
            if ligand_btn.selected and ligand_btn.unusable:
                ligand_btn.selected = False
        self.plugin.update_menu(self._menu)

    def populate_list(self, ui_list, complex_list, callback=None):
        ui_list.items = []
        for complex in complex_list:
            clone = self._pfb_complex.clone()
            btn = clone.get_content()
            btn.text.value.set_all(complex.full_name)
            btn.index = complex.index
            if callback:
                btn.register_pressed_callback(callback)
            ui_list.items.append(clone)
        self.plugin.update_content(ui_list)

    def hide_scores(self, show_ligand_names=False):
        self._to_display = False
        self._ls_results.items = []
        if show_ligand_names:
            for i in range(len(self._ligand_names)):
                clone = self._pfb_result.clone()
                lbl = clone._get_content()
                lbl.text_value = '{}: {}'.format(self._ligand_names[i], "Loading scores")
                self._ls_results.items.append(clone)
            self.plugin.update_content(self._ls_results)
        else:
            clone = self._pfb_result.clone()
            lbl = clone._get_content()
            lbl.text_value = ''
            self._ls_results.items.append(clone)
            self.plugin.update_content(self._ls_results)

    def update_ligand_scores(self, aggregate_score_list):
        results_list = self.ln_results.get_content()
        results_list.items = []
        scores = aggregate_score_list[0][0]
        for name, score in scores.items():
            clone = self._pfb_result.clone()
            lbl = clone._get_content()
            lbl.text_value = '{}: {}'.format(name, score)
            results_list.items.append(clone)
        self.plugin.update_content(results_list)