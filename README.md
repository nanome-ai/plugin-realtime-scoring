# Nanome - Realtime Scoring

Displays the docking score of selected molecules in the workspace, in realtime using DSX. Support for scoring multiple docking results concurrently is experimental, but available through the options menu.

Linux Support Only

## Dependencies

[Docker](https://docs.docker.com/get-docker/)

## Usage

To run Realtime Scoring in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh -a <plugin_server_address> [optional args]
```

---

In Nanome:

- Activate Plugin
- In the plugin window, select a receptor and ligands, then start scoring
- Plugin will display a list of all other complexes, with their docking score
- Moving a complex around will update its score

### Plugin Details

The following scores from DSX can be displayed:

- Total score: Total docking score for the ligand including possible torsion, sas and intramolecular contributions.
- Per contact score: Total score divided by the number of atom-atom-interactions having any contribution to the total score. (Do not confuse with a per atom score)

## Development

To run Realtime Scoring with autoreload:

```sh
$ python3 -m pip install -r requirements.txt
$ python3 run.py -r -a <plugin_server_address> [optional args]
```

## License

MIT
