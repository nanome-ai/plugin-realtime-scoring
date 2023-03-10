import nanome
from nanome.api.shapes import Shape, Sphere
from nanome.util import Logs, Color, enums
from nanome.util.enums import NotificationTypes

import os
import shlex
import subprocess
import tempfile

from .SettingsMenu import SettingsMenu
from .menu import MainMenu
from .dsx_parser import dsx_parse
from nanome.util import async_callback
from nanome.api import streams
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
            deep_comps[i].index = comp_index
            # Update coordinates to be relative to workspace
            mat = comp.get_complex_to_workspace_matrix()
            for atom in comp.atoms:
                atom.position = mat * atom.position
        receptor_comp = deep_comps[0]
        ligand_comps = deep_comps[1:]
        await self.calculate_scores(receptor_comp, ligand_comps)

    async def calculate_scores(self, receptor_comp, ligand_comps):
        # Generate sphere and attach streams
        spheres = self.generate_spheres(ligand_comps)
        await Shape.upload_multiple(spheres)
        self.sphere_indices = [sphere.index for sphere in spheres]
        # if self.settings._labels:
        #    self.create_writing_stream(self._atom_indices, enums.StreamType.label)
        self.color_stream, _ = await self.create_writing_stream(self.sphere_indices, enums.StreamType.shape_color)
        # create streams
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
            self.run_dsx(receptor_pdb.name, ligand_mol2.name, dsx_output_txt.name)
            # Make sure scores are added to ligand atoms
            dsx_parse(dsx_output_txt.name, ligand_comp)
            color_stream_data = self.get_color_stream_data(ligand_comp)
            self.color_stream.update(color_stream_data)
            Logs.message("Updated color stream")

    @staticmethod
    def get_color_stream_data(comp):
        """Generate color stream data based on atom scores.
        
        This assumes score values have been added to each atom,
        and atom_score_limits on the molecule
        """
        data = []
        for atom in comp.atoms:
            atom_score = getattr(atom, 'score', False)
            red = 0
            green = 0
            blue = 0
            alpha = 0
            ligand = atom.molecule
            if atom_score:
                denominator = -ligand.atom_score_limits[0] if atom.score < 0 else ligand.atom_score_limits[1]
                norm_score = atom.score / denominator
                red = 255 if norm_score > 0 else 0
                green = 255 if norm_score < 0 else 0
                alpha = int(140 * abs(norm_score))
                Logs.debug(f'{atom_score}: {alpha}')
            data.append(red)
            data.append(green)
            data.append(blue)
            data.append(alpha)
        return data

    @staticmethod
    def generate_spheres(ligand_comps):
        """Create a sphere for each atom on each ligand"""
        spheres = []
        for comp in ligand_comps:
            molecule_list = list(comp.molecules)
            curr_atoms = molecule_list[comp.current_frame].atoms
            for atom in curr_atoms:
                sphere = Sphere()
                sphere.color = Color(100, 100, 100, 0)
                sphere.radius = 1.0
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

    def on_advanced_settings(self):
        self.settings.open_menu()

    def on_stop(self):
        self.temp_dir.cleanup()

    def open_menu(self, menu=None):
        self.menu = self._menu
        self.menu.enabled = True
        self.update_menu(self.menu)

    @staticmethod
    def run_dsx(receptor_pdb, ligands_mol2, output_file):
        "Run DSX and write output to provided output_file."
        dsx_path = os.path.join(DIR, 'dsx', 'dsx_linux_64.lnx')
        dsx_args = [
            dsx_path,
            '-P', receptor_pdb,
            '-L', ligands_mol2,
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
        dsx_process.wait()
        dsx_output, _ = dsx_process.communicate()
        with open(output_file, 'w') as f:
            f.write(dsx_output)

    def on_complex_added(self):
        self.request_complex_list(self.update_lists)

    def on_complex_removed(self):
        self.request_complex_list(self.update_lists)

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
