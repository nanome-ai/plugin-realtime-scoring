# Nanome - Realtime Scoring

Displays the docking score of selected molecules in the workspace, in realtime using DSX. Support for scoring multiple docking results concurrently is experimental, but available through the options menu.

Linux Support Only

### Installation

```sh
$ pip install nanome-realtime-scoring
```

### Usage

To start the plugin:

```sh
$ nanome-realtime-scoring -a plugin_server_address
```

In Nanome:

- Activate Plugin
- In the plugin window, select a receptor and ligands, then start scoring
- Plugin will display a list of all other complexes, with their docking score
- Moving a complex around will update its score

### License

MIT
