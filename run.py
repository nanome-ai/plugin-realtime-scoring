import nanome
from plugin.RealtimeScoring import RealtimeScoring
from nanome.util import Color

def main(custom_data):
    description = "Display realtime scoring info about a selected ligand."
    plugin = nanome.Plugin("Realtime Scoring", description, "Scoring", True)
    plugin.set_plugin_class(RealtimeScoring)
    plugin.set_custom_data(custom_data)
    plugin.run()


if __name__ == "__main__":
    custom_data = {
        'color_positive_score': Color.White(),
        'color_negative_score': Color.Yellow(),
    }
    main(custom_data)