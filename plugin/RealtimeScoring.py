import nanome
from datetime import datetime, timedelta
from nanome.api import structure
from nanome.api.shapes import Shape, Sphere
from nanome.util import async_callback, Logs, Color, enums

from dsx import scoring_algo
from plugin.utils import ScoringOutputSchema
from plugin.SettingsMenu import SettingsMenu
from plugin.menu import MainMenu
from plugin import utils


class RealtimeScoring(nanome.AsyncPluginInstance):

    scoring_algorithm = scoring_algo.score_ligands

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.menu = MainMenu(self)
        self.settings = SettingsMenu(self)
        self.last_update = datetime.now()
        self.is_updating = False
        # Configurable settings
        self.color_negative_score = Color(0, 0, 255, 200)  # Blue
        self.color_positive_score = Color(255, 0, 0, 200)  # Red
        self.realtime_enabled = True
        # api structures
        self.receptor_index = None
        self.ligand_residue_indices = []
        self.complex_cache = {}

    def start(self):
        # Update settings based on custom data added at runtime.
        custom_data = {}
        if self.custom_data:
            try:
                custom_data = self.custom_data[0]
            except (IndexError, AttributeError):
                pass

        if custom_data.get('color_negative_score'):
            self.color_negative_score = custom_data.get('color_negative_score')
        if custom_data.get('color_positive_score'):
            self.color_positive_score = custom_data.get('color_positive_score')
        if custom_data.get('realtime_enabled'):
            self.realtime_enabled = custom_data.get('realtime_enabled')

    @async_callback
    async def on_run(self):
        complex_list = await self.request_complex_list()
        comp_indices = set([comp.index for comp in complex_list])
        self.complex_cache = await self.request_complexes(comp_indices)
        self.menu.render()

    @async_callback
    async def update(self):
        if not self.realtime_enabled:
            return
        update_time_secs = 3
        has_receptor = getattr(self, 'receptor_comp', None)
        has_ligands = getattr(self, 'ligand_residues', None)
        has_color_stream = getattr(self, 'color_stream', None)
        has_label_stream = getattr(self, 'label_stream', None)
        due_for_update = datetime.now() - self.last_update > timedelta(seconds=update_time_secs)
        if all([
            has_receptor, has_ligands, has_color_stream,
                has_label_stream, due_for_update, not self.is_updating]):
            Logs.debug("Updating cached ligands.")
            self.is_updating = True

            # Get updated complexes:
            lig_comp_indices = set()
            for res in self.ligand_residues:
                lig_comp_indices.add(res.complex.index)

            comp_indices = set([self.receptor_comp.index] + list(lig_comp_indices))
            updated_comps = await self.request_complexes(comp_indices)
            # If any of the complexes were deleted, destroy the streams
            if any([comp is None for comp in updated_comps]):
                Logs.message("Receptor or ligand deleted. Stopping streams.")
                # Use the score button so that the UI is updated.
                btn_score = self.menu.btn_score
                btn_score.selected = False
                self.menu.on_scoring_button_pressed(btn_score)
                self.receptor_index = None
                self.ligand_residue_indices = []
                return

            self.set_atoms_to_workspace_positions(updated_comps)
            self.last_update = datetime.now()
            # Check if positions have changed in workspace
            needs_rescore = False
            needs_stream_update = False
            for cached_comp, updated_comp in zip(self.complex_cache, updated_comps):
                position_changed = cached_comp.position.unpack() != updated_comp.position.unpack()
                rotation_changed = str(cached_comp.rotation) != str(updated_comp.rotation)
                atoms_changed = sum(1 for _ in cached_comp.atoms) != sum(1 for _ in updated_comp.atoms)
                if position_changed or rotation_changed:
                    needs_rescore = True
                if atoms_changed:
                    needs_rescore = True
                    needs_stream_update = True
                    break

            # Update cached complexes
            for i in range(len(self.complex_cache) - 1, -1, -1):
                cached_comp = self.complex_cache[i]
                if cached_comp.index in comp_indices:
                    updated_version = next(cmp for cmp in updated_comps if cmp.index == cached_comp.index)
                    if updated_version:
                        Logs.debug("Updating cached complex")
                        self.complex_cache[i] = updated_version
            # Update ligand residues with updated complexes
            if needs_stream_update:
                Logs.message("Receptor or ligand modified. Recreating streams.")
                await self.start_ligand_streams(self.ligand_atoms)
            if needs_rescore:
                Logs.message("Complex Positions changed. Rescoring Ligands.")
                await self.score_ligands()
            self.is_updating = False

    async def start_ligand_streams(self, ligand_atoms):
        """Set up streams and Shapes used for rendering scoring results."""
        if getattr(self, 'spheres', None):
            await Shape.destroy_multiple(self.spheres)
        self.spheres = self.generate_spheres(self.ligand_atoms)
        await Shape.upload_multiple(self.spheres)
        atom_indices = [atom.index for atom in ligand_atoms]
        sphere_indices = [sphere.index for sphere in self.spheres]
        self.label_stream, _ = await self.create_writing_stream(atom_indices, enums.StreamType.label)
        self.color_stream, _ = await self.create_writing_stream(sphere_indices, enums.StreamType.shape_color)
        self.size_stream, _ = await self.create_writing_stream(sphere_indices, enums.StreamType.sphere_shape_radius)

    @staticmethod
    def get_atoms(struct_list):
        """Yield all atoms from a list of structures."""
        for struc in struct_list:
            for atom in struc.atoms:
                yield atom

    @property
    def ligand_atoms(self):
        """Yield all atoms from all ligands."""
        return self.get_atoms(self.ligand_residues)

    @property
    def receptor_comp(self) -> structure.Complex:
        """Get the receptor complex."""
        return next((
            comp for comp in self.complex_cache
            if comp.index == self.receptor_index
        ), None)

    @property
    def ligand_residues(self):
        """Get the indices of all ligands."""
        residues = []
        for comp in self.complex_cache:
            for res in comp.residues:
                if res.index in self.ligand_residue_indices:
                    residues.append(res)
        return residues

    @staticmethod
    def set_atoms_to_workspace_positions(comp_list):
        """Set all atoms in a list of complexes to their positions in the workspace."""
        for comp in comp_list:
            # Update coordinates to be relative to workspace
            mat = comp.get_complex_to_workspace_matrix()
            for atom in comp.atoms:
                atom._old_position = atom.position
                atom.position = mat * atom.position

    async def setup_receptor_and_ligands(self, receptor_index, residue_indices):
        # Let's make sure we have deep receptor and ligand complexes
        self.receptor_index = receptor_index
        self.ligand_residue_indices = residue_indices
        await self.start_ligand_streams(self.ligand_atoms)

    async def score_ligands(self):
        if not getattr(self, 'receptor_comp', None):
            Logs.error("Receptor not set")
            return
        if not getattr(self, 'ligand_residues', None):
            Logs.error("Ligands not set")
            return
        if not getattr(self, 'ligand_residues', None):
            Logs.warning("Ligand Residues not specified")
            return
        score_data = await self.calculate_scores(self.receptor_comp, self.ligand_residues)

        # Concatenate all ligand results into a single list
        all_atom_scores = []
        for ligand_scores in score_data:
            all_atom_scores += ligand_scores['atom_scores']

        aggregate_scores = [score['aggregate_scores'] for score in score_data]
        await self.render_atom_scores(all_atom_scores)
        self.menu.update_ligand_scores(aggregate_scores)

    @classmethod
    async def calculate_scores(cls, receptor_comp, ligand_residues):
        # write ligand residues to separate complex
        ligand_comps = list(set([lig.complex for lig in ligand_residues]))
        for i, lig in enumerate(ligand_comps):
            ligand_comps[i] = utils.extract_residues_from_complex(lig, ligand_residues)
        ligand_scores = cls.scoring_algorithm(receptor_comp, ligand_comps)
        validation_errors = ScoringOutputSchema(many=True).validate(ligand_scores)
        if validation_errors:
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

        for atom in ligand_atoms:
            atom_score = score_dict.get(atom.index, False)
            atom_color = Color(0, 0, 0, 0)
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

    @async_callback
    async def on_complex_list_changed(self):
        comp_list = await self.request_complex_list()
        await self.menu.render(comp_list)
