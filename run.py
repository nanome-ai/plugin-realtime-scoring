import nanome
from plugin.RealtimeScoring import RealtimeScoring
# from nanome.util import Color
from plugin import __version__


if __name__ == "__main__":
    # Uncomment to customize the plugin behavior.
    # custom_data = {
    #     'color_positive_score': Color.Red(),
    #     'color_negative_score': Color.Blue(),
    #     'realtime_enabled': True
    # }
    plugin_name = 'Realtime Scoring'
    description = "Display realtime scoring info about a selected ligand."
    tag = 'Scoring'
    has_advanced_settings = True
    plugin = nanome.Plugin(plugin_name, description, tag, has_advanced_settings, version=__version__)
    plugin.set_plugin_class(RealtimeScoring)
    # plugin.set_custom_data(custom_data)
    plugin.run()
