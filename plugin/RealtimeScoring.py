import nanome
from nanome.api.shapes import Shape, Sphere
from nanome.util import Logs, Color, enums

import os
import shlex
import subprocess
import tempfile

from .SettingsMenu import SettingsMenu
from .menu import MainMenu
from .dsx_parser import dsx_parse
from nanome.util import async_callback
from nanome.api import streams
from nanome.api import structure


SDF_OPTIONS = structure.Complex.io.SDFSaveOptions()
SDF_OPTIONS.write_bonds = True
PDB_OPTIONS = structure.Complex.io.PDBSaveOptions()
PDB_OPTIONS.write_bonds = True

DIR = os.path.dirname(__file__)
RESULTS_PATH = os.path.join(DIR, 'dsx', 'results.txt')


class RealtimeScoring(nanome.AsyncPluginInstance):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.settings = SettingsMenu(self)
        self.temp_dir = tempfile.TemporaryDirectory()

    async def start_ligand_streams(self, ligand_atoms):
        """Set up streams and Shapes used for rendering scoring results."""
        spheres = self.generate_spheres(ligand_atoms)
        await Shape.upload_multiple(spheres)
        atom_indices = [atom.index for atom in ligand_atoms]
        self.label_stream, _ = await self.create_writing_stream(atom_indices, enums.StreamType.label)
        sphere_indices = [sphere.index for sphere in spheres]
        self.color_stream, _ = await self.create_writing_stream(sphere_indices, enums.StreamType.shape_color)
    
    @staticmethod
    def get_atoms(comp_list):
        """Return a chain object whose .__next__() method returns elements from the first iterable until it is exhausted,
        then elements from the next iterable, until all of the iterables are exhausted."""
        for comp in comp_list:
            for atom in comp.atoms:
                yield atom

    async def score_ligands(self, receptor_index, ligand_indices):
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
                atom._old_position = atom.position
                atom.position = mat * atom.position
        receptor_comp = deep_comps[0]
        ligand_comps = deep_comps[1:]

        comp_atoms = self.get_atoms(deep_comps)
        await self.start_ligand_streams(comp_atoms)
        await self.calculate_scores(receptor_comp, ligand_comps)
        # Update the sphere color around each atom
        comp_atoms = self.get_atoms(deep_comps)
        color_stream_data = self.get_color_stream_data(comp_atoms)
        self.color_stream.update(color_stream_data)
        # If update labels is turned on, update the label stream
        Logs.message("Updated color stream")
        if self.settings._labels:
            comp_atoms = self.get_atoms(deep_comps)
            label_stream_data = self.get_label_stream_data(comp_atoms)
            self.label_stream.update(label_stream_data)
            Logs.message("Updated label stream")

        # Reset ligand coordinates and ensure labels are enabled if settings say so
        for comp in deep_comps:
            for atom in comp.atoms:
                atom.position = atom._old_position
                if comp in ligand_comps and self.settings._labels:
                    atom.labeled = atom.label_text != 'N/A'
        self.update_structures_deep(deep_comps)

    async def calculate_scores(self, receptor_comp, ligand_comps):
        # Generate sphere and attach streams
        receptor_pdb = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
        receptor_comp.io.to_pdb(receptor_pdb.name, PDB_OPTIONS)
        # For each ligand, generate a PDB file and run DSX
        for ligand_comp in ligand_comps:
            ligand_sdf = tempfile.NamedTemporaryFile(suffix='.sdf')
            # Convert ligand pdb to a mol2.
            ligand_mol2 = tempfile.NamedTemporaryFile(suffix='.mol2')
            ligand_comp.io.to_sdf(ligand_sdf.name, SDF_OPTIONS)
            # Convert ligand from sdf to mol2. Not sure if we really need this
            self.nanobabel_convert(ligand_sdf.name, ligand_mol2.name)
            # Run DSX and retreive data from the subprocess.
            dsx_output_txt = tempfile.NamedTemporaryFile(suffix='.txt')
            self.run_dsx(receptor_pdb.name, ligand_mol2.name, dsx_output_txt.name)
            # Add new data to the ligand_comp.
            dsx_parse(dsx_output_txt.name, ligand_comp)

    @staticmethod
    def get_color_stream_data(atom_generator):
        """Generate color stream data based on atom scores.
        
        This assumes score values have been added to each atom,
        and atom_score_limits on the molecule
        """
        data = []
        for atom in atom_generator:
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
    def get_label_stream_data(atom_generator):
        """Generate color stream data based on atom scores.
        
        This assumes score values have been added to each atom,
        and atom_score_limits on the molecule
        """
        data = []
        for atom in atom_generator:
            atom_score = getattr(atom, 'score', None)
            label_text = str(round(atom_score, 2)) if atom_score else 'N/A'
            data.append(label_text)            
        return data

    @staticmethod
    def generate_spheres(ligand_atoms):
        """Create a sphere for each atom on each ligand"""
        spheres = []
        for atom in ligand_atoms:
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
            dsx_path, '-P', receptor_pdb, '-L', ligands_mol2, '-D', 'pdb_pot_0511',
            '-pp', '-F', output_file
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