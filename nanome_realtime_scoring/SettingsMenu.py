import nanome
from nanome.api.ui.image import Image as Image

import os

class SettingsMenu():
    def __init__(self, docking_plugin, closed_callback):
        self._plugin = docking_plugin
        self._score_all_frames = False

        self._menu = nanome.ui.Menu.io.from_json(os.path.join(os.path.dirname(__file__), 'settings.json'))
        self._menu.register_closed_callback(closed_callback)

        self._btn_score_all_frames = self._menu.root.find_node('Button').get_content()
        self._btn_score_all_frames.register_pressed_callback(self.toggle_score_all_frames)
    
    def open_menu(self, menu=None):
        self._plugin.menu = self._menu
        self._plugin.menu.enabled = True
        self._plugin.update_menu(self._plugin.menu)
    
    def score_all_frames(self):
        return self._score_all_frames

    def toggle_score_all_frames(self, button):
        self._score_all_frames = not self._score_all_frames
        text = "on" if self._score_all_frames else "off"
        self._btn_score_all_frames.set_all_text(text)
        self._plugin.update_content(self._btn_score_all_frames)