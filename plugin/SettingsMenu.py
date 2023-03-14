import nanome
import os
from nanome.api import ui


BASE_PATH = os.path.dirname(f'{os.path.realpath(__file__)}')
MENU_PATH = os.path.join(BASE_PATH, 'menu_json', 'settings.json')


class SettingsMenu():

    def __init__(self, plugin):
        self._plugin = plugin
        self._menu = nanome.ui.Menu.io.from_json(MENU_PATH)
        self._menu.index = 120  # arbitrary

        self._btn_labels: ui.Button = self._menu.root.find_node('AtomLabelsButton').get_content()
        self._btn_labels.toggle_on_press = True
        self._btn_labels.selected = False

        self._btn_score_all_frames: ui.Button = self._menu.root.find_node('AllFramesButton').get_content()
        self._btn_score_all_frames.toggle_on_press = True
        self._btn_score_all_frames.selected = False

        self._btn_total: ui.Button = self._menu.root.find_node("Total Button", True).get_content()
        self._btn_total.toggle_on_press = True
        self._btn_total.selected = False

        self._btn_pcs: ui.Button = self._menu.root.find_node("PCS Button", True).get_content()
        self._btn_pcs.toggle_on_press = True
        self._btn_pcs.selected = False

    def open_menu(self, menu=None):
        self._plugin.menu = self._menu
        self._plugin.menu.enabled = True
        self._plugin.update_menu(self._plugin.menu)

    @property
    def update_labels(self):
        return self._btn_labels.selected

    @property
    def score_all_frames(self):
        return self._btn_score_all_frames.selected
    
    @property
    def show_total_scores(self):
        return self._btn_total.selected
    
    @property
    def per_contact_scores(self):
        return self._btn_pcs.selected
    
