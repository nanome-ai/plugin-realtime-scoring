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
        self._p_selection: ui.LayoutNode = self._menu.root.find_node("Selection Panel", True)
        self._p_selection.enabled = True
        self._p_results: ui.LayoutNode = self._menu.root.find_node("Results Panel", True)
        self._ls_receptors: ui.UIList = self._menu.root.find_node("Receptor List", True).get_content()
        self._ls_ligands: ui.UIList = self._menu.root.find_node("Ligands List", True).get_content()
        self._ls_results: ui.UIList = self._p_results.get_content()
        self.btn_score: ui.Button = self._menu.root.find_node("btn_score", True).get_content()
        self.btn_score.register_pressed_callback(self.button_pressed)
        self.btn_score.toggle_on_press = True

        self._pfb_complex = nanome.ui.LayoutNode()
        self._pfb_complex.add_new_button()

        self._pfb_result = nanome.ui.LayoutNode()
        self._pfb_result.add_new_label()

    @async_callback
    async def button_pressed(self, button):
        receptor_index = self._receptor_index
        ligand_indices = self._ligand_indices
        if button.selected:
            await self.start_scoring(receptor_index, ligand_indices)
        else:
            self.stop_scoring()

    async def start_scoring(self, receptor_index, ligand_indices):
        Logs.message("Start Scoring")
        self._p_selection.enabled = False
        self._p_results.enabled = True

        if receptor_index is None:
            self.send_notification(NotificationTypes.error, "Please select a receptor")
            return

        if len(self._ligand_indices) == 0:
            self.send_notification(NotificationTypes.error, "Please select at least one ligand")
            return

        # Add loading message to results panel
        results_list = self._p_results.get_content()
        results_list.items = []
        clone = self._pfb_result.clone()
        lbl = clone._get_content()
        lbl.text_value = 'Loading...'
        results_list.items.append(clone)
        self.plugin.update_menu(self._menu)
        
        await self.plugin.setup_receptor_and_ligands(receptor_index, ligand_indices)
        await self.plugin.score_ligands()
        self._menu.title = "Scores"
        # self.hide_scores()

    def stop_scoring(self):
        self.plugin.stop_streams()
        self._p_selection.enabled = True
        self._p_results.enabled = False
        self.plugin.update_menu(self._menu)

    @async_callback
    async def render(self, complex_list, force_enable=False):
        self.update_lists(complex_list)
        if force_enable:
            self._menu.enabled = True
        self.plugin.update_menu(self._menu)

    def update_lists(self, complex_list):

        self._shallow_complexes = complex_list

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

            self.plugin.update_menu(self._menu)

        def ligand_pressed(ligand):
            ligand.selected = not ligand.selected
            self.plugin.update_content(ligand)
            update_selected_ligands()

        def populate_list(ls, cb):
            ls.items = []
            for complex in complex_list:
                clone = self._pfb_complex.clone()
                btn = clone.get_content()
                btn.text.value.set_all(complex.full_name)
                btn.index = complex.index
                btn.register_pressed_callback(cb)
                ls.items.append(clone)
            self.plugin.update_content(ls)

        self._receptor_index = None
        self._ligand_indices = []
        populate_list(self._ls_receptors, receptor_pressed)
        populate_list(self._ls_ligands, ligand_pressed)

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
        results_list = self._p_results.get_content()
        results_list.items = []
        scores = aggregate_score_list[0][0]
        for name, score in scores.items():
            clone = self._pfb_result.clone()
            lbl = clone._get_content()
            lbl.text_value = '{}: {}'.format(name, score)
            results_list.items.append(clone)
        self.plugin.update_content(results_list)