# Nanome - Realtime Scoring

Displays the docking score of all molecules in the workspace, in realtime.

Runs only on Linux for now

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
- In the plugin window, select a receptor, and start scoring
- Plugin will display a list of all other complexes, with their docking score
- Moving a complex around will update its score

### License

MIT