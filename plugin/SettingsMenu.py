import nanome
import os
from nanome.api.ui.image import Image as Image


BASE_PATH = os.path.dirname(f'{os.path.realpath(__file__)}')
MENU_PATH = os.path.join(BASE_PATH, 'menu_json', 'settings.json')


class SettingsMenu():

    def __init__(self, plugin):
        self._plugin = plugin
        self._labels = True
        self._score_all_frames = False
        self.show_total = True
        self.show_pcs = True
        self._menu = nanome.ui.Menu.io.from_json(MENU_PATH)
        self._menu.index = 120  # arbitrary
        # self._menu.register_closed_callback(closed_callback)

        self._btn_labels = self._menu.root.find_node('AtomLabelsButton').get_content()
        self._btn_labels.register_pressed_callback(self.toggle_labels)
        self._btn_score_all_frames = self._menu.root.find_node('AllFramesButton').get_content()
        self._btn_score_all_frames.register_pressed_callback(self.toggle_score_all_frames)

        def total_pressed(button):
            if self.show_total:
                button.set_all_text("off")
                self.show_total = False
            else:
                button.set_all_text("on")
                self.show_total = True
            self._plugin.update_content(button)

        def pcs_pressed(button):
            if self.show_total:
                button.set_all_text("off")
                self.show_pcs = False
            else:
                button.set_all_text("on")
                self.show_pcs = True
            self._plugin.update_content(button)

        self._btn_total = self._menu.root.find_node("Total Button", True).get_content()
        self._btn_total.register_pressed_callback(total_pressed)
        self._btn_pcs = self._menu.root.find_node("PCS Button", True).get_content()
        self._btn_pcs.register_pressed_callback(pcs_pressed)

    def open_menu(self, menu=None):
        self._plugin.menu = self._menu
        self._plugin.menu.enabled = True
        self._plugin.update_menu(self._plugin.menu)

    def score_all_frames(self):
        return self._score_all_frames

    def toggle_labels(self, button):
        self._labels = not self._labels
        text = "on" if self._labels else "off"
        self._btn_labels.set_all_text(text)
        self._plugin.update_content(self._btn_labels)

    def toggle_score_all_frames(self, button):
        self._score_all_frames = not self._score_all_frames
        text = "on" if self._score_all_frames else "off"
        self._btn_score_all_frames.set_all_text(text)
        self._plugin.update_content(self._btn_score_all_frames)
