import nanome
from nanome.api.shapes import Shape, Sphere
from nanome.util import Logs, Color, enums
from nanome.util.enums import NotificationTypes

import os
import shlex
import functools
import subprocess
import tempfile

from .SettingsMenu import SettingsMenu
from .menu import MainMenu
from .dsx_parser import dsx_parse
from nanome.util import async_callback

SDF_OPTIONS = nanome.api.structure.Complex.io.SDFSaveOptions()
SDF_OPTIONS.write_bonds = True
PDB_OPTIONS = nanome.api.structure.Complex.io.PDBSaveOptions()
PDB_OPTIONS.write_bonds = True

DIR = os.path.dirname(__file__)
RESULTS_PATH = os.path.join(DIR, 'dsx', 'results.txt')


class RealtimeScoring(nanome.AsyncPluginInstance):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._is_running = False
        self.settings = SettingsMenu(self)
        self.temp_dir = tempfile.TemporaryDirectory()

    async def score_ligand(self, receptor_index, ligand_indices):
        # Get latest version of receptor and ligand
        deep_comps = await self.request_complexes([receptor_index, *ligand_indices])
        for i in range(0, len(deep_comps)):
            comp = deep_comps[i]
            comp_index = comp.index
            comp.locked = True
            # deep_comps[i] = comp.convert_to_frames()
            deep_comps[i].index = comp_index
        # self.update_structures_deep(deep_comps)
        receptor_comp = deep_comps[0]
        ligand_comps = deep_comps[1:]
        # Generate sphere streams
        spheres = self.generate_spheres(ligand_comps)
        receptor_pdb = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
        receptor_comp.io.to_pdb(receptor_pdb.name, PDB_OPTIONS)
        # For each ligand, generate a PDB file and run DSX
        for ligand_comp in ligand_comps:
            ligand_sdf = tempfile.NamedTemporaryFile(suffix='.sdf')
            # Convert ligand pdb to a mol2.
            ligand_mol2 = tempfile.NamedTemporaryFile(suffix='.mol2')
            ligand_comp.io.to_sdf(ligand_sdf.name, SDF_OPTIONS)
            dsx_output_txt = tempfile.NamedTemporaryFile(suffix='.txt')
            self.nanobabel_convert(ligand_sdf.name, ligand_mol2.name)
            dsx_popen = self.dsx_start(receptor_pdb.name, ligand_mol2.name, dsx_output_txt.name)
            dsx_popen.wait()
            # Make sure scores are added to ligand atoms
            dsx_parse(dsx_output_txt.name, ligand_comp)
            assert any(hasattr(atom, 'score') for atom in ligand_comp.atoms), 'No scores found in DSX output'

    @staticmethod
    def generate_spheres(ligand_comps):
        spheres = []
        for comp in ligand_comps:
            molecule_list = list(comp.molecules)
            curr_atoms = molecule_list[comp.current_frame].atoms
            for atom in curr_atoms:
                sphere = Sphere()
                sphere.color = Color(100, 100, 100, 120)
                sphere.radius = 1.1
                anchor = sphere.anchors[0]
                anchor.anchor_type = enums.ShapeAnchorType.Atom
                anchor.target = atom.index
                spheres.append(sphere)
        return spheres

    @async_callback
    async def start(self):
        self.menu = MainMenu(self)

    @async_callback
    async def on_run(self):
        complex_list = await self.request_complex_list()
        self.menu.render(complex_list)

    def start_scoring(self, receptor_index, ligand_indices):
        self._is_running = True
        self._nanobabel_running = False
        self._dsx_running = False

        self._scores_ready = False
        self._creating_streams = False
        self._stop_after_deleting_spheres = False

        self._respond_to_update = False
        self.get_full_complexes(receptor_index, ligand_indices)

    def on_advanced_settings(self):
        self.settings.open_menu()

    def on_stop(self):
        self.temp_dir.cleanup()

    def open_menu(self, menu=None):
        self.menu = self._menu
        self.menu.enabled = True
        self.update_menu(self.menu)

    @staticmethod
    def dsx_start(receptor_pdb, ligands_pdb, output_file):
        dsx_path = os.path.join(DIR, 'dsx', 'dsx_linux_64.lnx')
        dsx_args = [
            dsx_path,
            '-P', receptor_pdb,
            '-L', ligands_pdb,
            '-D', 'pdb_pot_0511',
            '-pp',
            '-F',
            output_file
        ]
        try:
            dsx_process = subprocess.Popen(dsx_args, cwd=os.path.join(DIR, 'dsx'), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        except Exception:
            nanome.util.Logs.error("Couldn't execute dsx, please check if executable is in the plugin folder and has permissions. Try executing chmod +x " + dsx_path)
            return
        return dsx_process

    def on_complex_added(self):
        self.request_complex_list(self.update_lists)

    def on_complex_removed(self):
        self.request_complex_list(self.update_lists)

    def complex_updated(self, complex_list, complex):
        if not self._is_running:
            return
        if not self._respond_to_update:
            return
        molecule_list = list(complex.molecules)
        atom_count = len(list(molecule_list[complex.current_frame].atoms))
        if self._ligand_frames[complex.index] != complex.current_frame or \
           self._ligand_atom_counts[complex.index] != atom_count:
            if self._uploading_spheres:
                self._respond_to_update = False
                self._delete_spheres = True
                # self.menu.hide_scores(True)
            elif self._creating_streams:
                self._respond_to_update = False
                self._creating_streams = False
                self.clear_sphere_streams()
                # self.menu.hide_scores(True)
                self.get_full_complexes()
            else:
                self._respond_to_update = True

    def nanobabel_convert(self, input_file, output_file):
        cmd = f'nanobabel convert -i {input_file} -o {output_file}'
        args = shlex.split(cmd)
        try:
            popen = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            popen.wait()
        except Exception:
            nanome.util.Logs.error("Couldn't execute nanobabel, please check if packet 'openbabel' is installed")
            return


def main():
    description = "Display realtime scoring info about a selected ligand."
    plugin = nanome.Plugin("Realtime Scoring", description, "Scoring", True)
    plugin.set_plugin_class(RealtimeScoring)
    plugin.run()


if __name__ == "__main__":
    main()
