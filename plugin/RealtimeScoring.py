import nanome
import os
import tempfile
from datetime import datetime, timedelta
from nanome.api import structure
from nanome.api.shapes import Shape, Sphere
from nanome.util import async_callback, Logs, Color, enums

from dsx import scoring_algo
from scoring_schema import ScoringOutputSchema
from .SettingsMenu import SettingsMenu
from .menu import MainMenu


SDF_OPTIONS = structure.Complex.io.SDFSaveOptions()
SDF_OPTIONS.write_bonds = True
PDB_OPTIONS = structure.Complex.io.PDBSaveOptions()
PDB_OPTIONS.write_bonds = True

DIR = os.path.dirname(__file__)
RESULTS_PATH = os.path.join(DIR, 'dsx', 'results.txt')


class RealtimeScoring(nanome.AsyncPluginInstance):

    scoring_algorithm = scoring_algo.score_ligands

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.main_menu = MainMenu(self)
        self.settings = SettingsMenu(self)
        self.temp_dir = tempfile.TemporaryDirectory()
        self.last_update = datetime.now()
        self.is_updating = False

        try:
            custom_data = self._custom_data[0]
        except IndexError:
            custom_data = {}
        default_negative_color = Color.Blue()
        default_positive_color = Color.Red()
        self.color_positive_score = custom_data.get('color_positive_score', default_positive_color)
        self.color_negative_score = custom_data.get('color_negative_score', default_negative_color)

    @async_callback
    async def on_run(self):
        complex_list = await self.request_complex_list()
        self.main_menu.render(complex_list)

    @async_callback
    async def update(self):
        update_time_secs = 3
        has_receptor = getattr(self, 'receptor_comp', None)
        has_ligands = getattr(self, 'ligand_comps', None)
        has_color_stream = getattr(self, 'color_stream', None)
        has_label_stream = getattr(self, 'label_stream', None)
        due_for_update = datetime.now() - self.last_update > timedelta(seconds=update_time_secs)
        if all([
            has_receptor, has_ligands, has_color_stream,
                has_label_stream, due_for_update, not self.is_updating]):
            Logs.debug("Updating cached ligands.")
            self.is_updating = True
            comp_indices = [self.receptor_comp.index] + [cmp.index for cmp in self.ligand_comps]
            updated_comps = await self.request_complexes(comp_indices)
            self.set_atoms_to_workspace_positions(updated_comps)
            self.last_update = datetime.now()
            # Check if positions have changed in workspace
            needs_rescore = False
            needs_stream_update = False
            cached_comps = [self.receptor_comp, *self.ligand_comps]
            for cached_comp, updated_comp in zip(cached_comps, updated_comps):
                position_changed = cached_comp.position.unpack() != updated_comp.position.unpack()
                rotation_changed = str(cached_comp.rotation) != str(updated_comp.rotation)
                atoms_changed = sum(1 for _ in cached_comp.atoms) != sum(1 for _ in updated_comp.atoms)
                if position_changed or rotation_changed:
                    needs_rescore = True
                if atoms_changed:
                    needs_rescore = True
                    needs_stream_update = True
                    break
            self.receptor_comp = updated_comps[0]
            self.ligand_comps = updated_comps[1:]
            if needs_stream_update:
                Logs.debug("Ligand has changed. Updating streams.")
                await self.start_ligand_streams(self.ligand_atoms)
            if needs_rescore:
                Logs.debug("Rescoring Ligands.")
                await self.score_ligands()    
            self.is_updating = False

    async def start_ligand_streams(self, ligand_atoms):
        """Set up streams and Shapes used for rendering scoring results."""
        if hasattr(self, 'spheres') and getattr(self, 'spheres', None):
            await Shape.destroy_multiple(self.spheres)
        self.spheres = self.generate_spheres(self.ligand_atoms)
        await Shape.upload_multiple(self.spheres)
        atom_indices = [atom.index for atom in ligand_atoms]
        sphere_indices = [sphere.index for sphere in self.spheres]
        self.label_stream, _ = await self.create_writing_stream(atom_indices, enums.StreamType.label)
        self.color_stream, _ = await self.create_writing_stream(sphere_indices, enums.StreamType.shape_color)
        self.size_stream, _ = await self.create_writing_stream(sphere_indices, enums.StreamType.sphere_shape_radius)

    @staticmethod
    def get_atoms(comp_list):
        """Yield all atoms from a list of complexes."""
        for comp in comp_list:
            for atom in comp.atoms:
                yield atom

    @property
    def ligand_atoms(self):
        """Yield all atoms from all ligands."""
        return self.get_atoms(self.ligand_comps)

    @staticmethod
    def set_atoms_to_workspace_positions(comp_list):
        """Set all atoms in a list of complexes to their positions in the workspace."""
        for comp in comp_list:
            # Update coordinates to be relative to workspace
            mat = comp.get_complex_to_workspace_matrix()
            for atom in comp.atoms:
                atom._old_position = atom.position
                atom.position = mat * atom.position

    async def setup_receptor_and_ligands(self, receptor_index, ligand_indices):
        deep_comps = await self.request_complexes([receptor_index, *ligand_indices])
        self.set_atoms_to_workspace_positions(deep_comps)
        self.receptor_comp = deep_comps[0]
        self.ligand_comps = deep_comps[1:]
        await self.start_ligand_streams(self.ligand_atoms)

    async def score_ligands(self):
        if not getattr(self, 'receptor_comp', None):
            Logs.error("Receptor not set")
            return
        if not getattr(self, 'ligand_comps', None):
            Logs.error("Ligands not set")
            return
        score_data = await self.calculate_scores(self.receptor_comp, self.ligand_comps)

        # Concatenate all ligand results into a single list
        all_atom_scores = []
        for ligand_scores in score_data:
            all_atom_scores += ligand_scores['atom_scores']
        
        aggregate_scores = [score['aggregate_scores'] for score in score_data]
        await self.render_atom_scores(all_atom_scores)
        self.main_menu.update_ligand_scores(aggregate_scores)

    @classmethod
    async def calculate_scores(cls, receptor_comp, ligand_comps):
        ligand_scores = cls.scoring_algorithm(receptor_comp, ligand_comps)
        validation_errors = ScoringOutputSchema(many=True).validate(ligand_scores)
        if validation_errors:
            # Logs.error("Validation errors: ", validation_errors)
            raise ValueError("Validation errors: ", validation_errors)
        return ligand_scores

    async def render_atom_scores(self, score_data):
        # Update the sphere color around each atom
        radius_data = self.get_radius_stream_data(score_data, self.ligand_atoms)
        self.size_stream.update(radius_data)
        Logs.message("Updated radius stream")
        color_stream_data = self.get_color_stream_data(score_data, self.ligand_atoms)
        self.color_stream.update(color_stream_data)
        # If update labels is turned on, update the label stream
        Logs.message("Updated color stream")
        if self.settings.update_labels:
            label_stream_data = self.get_label_stream_data(score_data, self.ligand_atoms)
            self.label_stream.update(label_stream_data)
            Logs.message("Updated label stream")

    def get_color_stream_data(self, score_data, ligand_atoms):
        """Generate color stream data based on atom scores.

        Note: Technically the stream is tied to a Sphere, but
        there is a sphere corresponding to every ligand atom, 
        so using the atom index works just as well.
        """
        data = []
        score_dict = {atom_index: atom_score for atom_index, atom_score in score_data}
        try:
            min_score = min(score_dict.values())
            max_score = max(score_dict.values())
        except ValueError:
            # No scored atoms, still update colors
            pass

        alpha = 200
        for atom in ligand_atoms:
            atom_score = score_dict.get(atom.index, False)
            atom_color = Color(0, 0, 0, alpha)
            if atom_score:
                denominator = -min_score if atom_score < 0 else max_score
                norm_score = atom_score / denominator
                if norm_score > 0:
                    atom_color = self.color_positive_score
                elif norm_score < 0:
                    atom_color = self.color_negative_score
            data.extend(atom_color.rgba)
        return data

    @staticmethod
    def get_label_stream_data(score_data, ligand_atoms):
        """Generate color stream data based on atom scores."""
        score_dict = {
            atom_index: atom_score
            for atom_index, atom_score in score_data
        }
        data = []
        for atom in ligand_atoms:
            atom_score = score_dict.get(atom.index, False)
            label_text = str(round(atom_score, 2)) if atom_score else ''
            data.append(label_text)
        return data

    @staticmethod
    def get_radius_stream_data(score_data, ligand_atoms):
        """Generate data to send to sphere radius stream."""
        data = []
        score_dict = {
            atom_index: atom_score
            for atom_index, atom_score in score_data
        }
        try:
            min_score = min(score_dict.values())
            max_score = max(score_dict.values())
        except ValueError:
            # No scored atoms, still update colors
            pass

        max_radius = 0.9
        min_radius = 0.4
        for atom in ligand_atoms:
            atom_score = score_dict.get(atom.index, False)
            if not atom_score:
                data.append(0)
                continue
            denominator = min_score if atom_score < 0 else max_score
            norm_score = abs(atom_score / denominator)
            radius = max(norm_score * max_radius, min_radius)
            data.append(radius)
        return data

    @staticmethod
    def generate_spheres(ligand_atoms):
        """Create a sphere for each atom on each ligand"""
        spheres = []
        for atom in ligand_atoms:
            sphere = Sphere()
            sphere.color = Color(100, 100, 100, 0)
            sphere.radius = .75
            anchor = sphere.anchors[0]
            anchor.anchor_type = enums.ShapeAnchorType.Atom
            anchor.target = atom.index
            spheres.append(sphere)
        return spheres

    def stop_streams(self):
        self.color_stream.destroy()
        self.label_stream.destroy()
        self.size_stream.destroy()
        self.color_stream = None
        self.label_stream = None
        self.size_stream = None

    def on_advanced_settings(self):
        self.settings.open_menu()

    def on_stop(self):
        self.temp_dir.cleanup()

    def open_menu(self, menu=None):
        self.main_menu = self._menu
        self.main_menu.enabled = True
        self.update_menu(self.main_menu)

    def on_complex_added(self):
        self.request_complex_list(self.update_lists)

    def on_complex_removed(self):
        self.request_complex_list(self.update_lists)
