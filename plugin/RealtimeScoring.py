import nanome
from nanome.api.shapes import Shape, Sphere
from nanome.util import Logs, Color, enums

import os
import tempfile

from .SettingsMenu import SettingsMenu
from .menu import MainMenu
from . import dsx_scoring
from nanome.util import async_callback
from nanome.api import structure


SDF_OPTIONS = structure.Complex.io.SDFSaveOptions()
SDF_OPTIONS.write_bonds = True
PDB_OPTIONS = structure.Complex.io.PDBSaveOptions()
PDB_OPTIONS.write_bonds = True

DIR = os.path.dirname(__file__)
RESULTS_PATH = os.path.join(DIR, 'dsx', 'results.txt')


class RealtimeScoring(nanome.AsyncPluginInstance):

    scoring_algorithm = dsx_scoring.score_ligands

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.settings = SettingsMenu(self)
        self.temp_dir = tempfile.TemporaryDirectory()

    async def start_ligand_streams(self, ligand_atoms):
        """Set up streams and Shapes used for rendering scoring results."""
        spheres = self.generate_spheres(ligand_atoms)
        await Shape.upload_multiple(spheres)
        atom_indices = [atom.index for atom in ligand_atoms]
        sphere_indices = [sphere.index for sphere in spheres]
        self.label_stream, _ = await self.create_writing_stream(atom_indices, enums.StreamType.label)
        self.color_stream, _ = await self.create_writing_stream(sphere_indices, enums.StreamType.shape_color)

    @staticmethod
    def get_atoms(comp_list):
        """Return a chain object whose .__next__() method returns elements from the first iterable until it is exhausted,
        then elements from the next iterable, until all of the iterables are exhausted."""
        for comp in comp_list:
            for atom in comp.atoms:
                yield atom
    
    @property
    def ligand_atoms(self):
        """Return a chain object whose .__next__() method returns elements from the first iterable until it is exhausted,
        then elements from the next iterable, until all of the iterables are exhausted."""
        return self.get_atoms(self.ligand_comps)

    async def setup_receptor_and_ligands(self, receptor_index, ligand_indices):
        deep_comps = await self.request_complexes([receptor_index, *ligand_indices])
        for comp in deep_comps:
            comp.locked = True
            # Update coordinates to be relative to workspace
            mat = comp.get_complex_to_workspace_matrix()
            for atom in comp.atoms:
                atom._old_position = atom.position
                atom.position = mat * atom.position
        self.receptor_comp = deep_comps[0]
        self.ligand_comps = deep_comps[1:]
        await self.start_ligand_streams(self.ligand_atoms)
        await self.score_ligands()

    async def score_ligands(self):
        if not getattr(self, 'receptor_comp', None):
            Logs.error("Receptor not set")
            return
        if not getattr(self, 'ligand_comps', None):
            Logs.error("Ligands not set")
            return
        score_data = await self.calculate_scores(self.receptor_comp, self.ligand_comps)
        await self.render_atom_scores(score_data)
        # Reset ligand coordinates and ensure labels are enabled if settings say so
        # for comp in deep_comps:
        #     for atom in comp.atoms:
        #         atom.position = atom._old_position
        #         if comp in ligand_comps and self.settings._labels:
        #             atom.labeled = atom.label_text != 'N/A'
        # self.update_structures_deep(deep_comps)

    @classmethod
    async def calculate_scores(cls, receptor_comp, ligand_comps):
        atom_scores = dsx_scoring.score_ligands(receptor_comp, ligand_comps)
        return atom_scores

    async def render_atom_scores(self, score_data):
        # Update the sphere color around each atom
        color_stream_data = self.get_color_stream_data(score_data, self.ligand_atoms)
        self.color_stream.update(color_stream_data)
        # If update labels is turned on, update the label stream
        Logs.message("Updated color stream")
        if self.settings._labels:
            label_stream_data = self.get_label_stream_data(score_data, self.ligand_atoms)
            self.label_stream.update(label_stream_data)
            Logs.message("Updated label stream")
    
    @staticmethod
    def get_color_stream_data(score_data, ligand_atoms):
        """Generate color stream data based on atom scores.
        
        This assumes score values have been added to each atom,
        and atom_score_limits on the molecule
        score_data: list of tuples of (atom, atom_score)
        """
        data = []
        scored_indices = {atom.index: atom_score for atom, atom_score in score_data}
        min_score = min(scored_indices.values())
        max_score = max(scored_indices.values())
        for atom in ligand_atoms:
            red = 0
            green = 0
            blue = 0
            alpha = 0
            atom_score = scored_indices.get(atom.index, False)
            if atom_score:
                denominator = -min_score if atom_score < 0 else max_score
                norm_score = atom_score / denominator
                red = 255 if norm_score > 0 else 0
                green = 255 if norm_score < 0 else 0
                alpha = int(140 * abs(norm_score))
                Logs.debug(f'{atom_score}: {alpha}')
            data.extend((red, green, blue, alpha))
        return data

    @staticmethod
    def get_label_stream_data(score_data, ligand_atoms):
        """Generate color stream data based on atom scores.
        
        This assumes score values have been added to each atom,
        and atom_score_limits on the molecule
        """
        scored_indices = {atom.index: atom_score for atom, atom_score in score_data}
        data = []
        for atom in ligand_atoms:
            atom_score = scored_indices.get(atom.index, False)
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

    def on_complex_added(self):
        self.request_complex_list(self.update_lists)

    def on_complex_removed(self):
        self.request_complex_list(self.update_lists)


def main():
    description = "Display realtime scoring info about a selected ligand."
    plugin = nanome.Plugin("Realtime Scoring", description, "Scoring", True)
    plugin.set_plugin_class(RealtimeScoring)
    plugin.run()


if __name__ == "__main__":
    main()