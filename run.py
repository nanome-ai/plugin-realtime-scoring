import nanome
from plugin.RealtimeScoring import RealtimeScoring
from nanome.util import Color


if __name__ == "__main__":
    custom_data = {
        'color_positive_score': Color(255, 0, 0, 200),
        'color_negative_score': Color(0, 0, 255, 200),
        'realtime_enabled': True
    }
    plugin_name = 'Realtime Scoring'
    description = "Display realtime scoring info about a selected ligand."
    tag = 'Scoring'
    has_advanced_settings = True
    plugin = nanome.Plugin(plugin_name, description, tag, has_advanced_settings)
    plugin.set_plugin_class(RealtimeScoring)
    plugin.set_custom_data(custom_data)
    plugin.run()
