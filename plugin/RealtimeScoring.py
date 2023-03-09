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

BFACTOR_MIN = 0
BFACTOR_MAX = 50

BFACTOR_GAP = BFACTOR_MAX - BFACTOR_MIN
BFACTOR_MID = BFACTOR_GAP / 2
BFACTOR_HALF = BFACTOR_MAX - BFACTOR_MID

DIR = os.path.dirname(__file__)
RESULTS_PATH = os.path.join(DIR, 'dsx', 'results.txt')


class RealtimeScoring(nanome.AsyncPluginInstance):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._is_running = False
        self.settings = SettingsMenu(self)

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
        receptor_pdb = tempfile.NamedTemporaryFile(suffix='.pdb')
        receptor_comp.io.to_pdb(receptor_pdb.name, PDB_OPTIONS)
        # For each ligand, generate a PDB file and run DSX
        for ligand_comp in ligand_comps:
            ligand_sdf = tempfile.NamedTemporaryFile(suffix='.sdf')
            ligand_mol2 = tempfile.NamedTemporaryFile(suffix='.mol2')
            dsx_output_txt = tempfile.NamedTemporaryFile(suffix='.txt')
            ligand_comp.io.to_sdf(ligand_sdf.name, SDF_OPTIONS)
            # Convert ligand pdb to a mol2.
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
        self._protein_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._ligands_converted = tempfile.NamedTemporaryFile(delete=False, suffix=".mol2")
        self._is_running = False
        self._menu = MainMenu(self)
        self.menu = self._menu
        self._nanobabel_running = False
        self._dsx_running = False

        self._ligands = None
        self._receptor_index = None
        self._ligand_indices = []
        self._complexes = []
        self._spheres = []
        self._sphere_count = 0
        self._atom_count = 0
        self._ligand_atom_counts = {}
        self._ligand_frames = {}
        self._label_stream = None
        self._color_stream = None
        self._scale_stream = None

        self._is_button_loading = False
        self._streams_ready = False
        self._creating_streams = False
        self._uploading_spheres = False
        self._delete_spheres = False
        self._stop_after_deleting_spheres = False
        self._respond_to_update = True
        self.complex_list = await self.request_complex_list()

    def on_run(self):
        self.menu.render(self.complex_list)

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
        os.remove(self._protein_input.name)
        os.remove(self._ligands_input.name)
        os.remove(self._ligands_converted.name)

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

    def stop_scoring(self):
        self.freeze_button()
        self._receptor_index = None
        self._ligand_indices = []

        self._is_running = False
        self._nanobabel_running = False
        self._dsx_running = False

        if self._uploading_spheres:
            self._delete_spheres = True
            self._stop_after_deleting_spheres = True
        else:
            self.clear_sphere_streams()

        self.menu.title = "Realtime Scoring"
        self._p_results.enabled = False
        self._is_button_loading = False
        self._btn_score.text.value.set_all("Start scoring")
        self.request_complex_list(self.update_lists)
        self.update_menu(self.menu)
        self.menu.unfreeze_button("Start scoring")

    def get_full_complexes(self, receptor_index, ligand_indices):

        def set_complexes(complex_list):
            self._respond_to_update = True
            self._complexes = complex_list
            self._ligand_atom_counts = {}
            self._ligand_frames = {}
            atom_counts = [0 for _ in range(len(complex_list))]
            for complex_i in range(0, len(complex_list)):
                complex = complex_list[complex_i]
                complex.locked = True
                complex = complex_list[complex_i].convert_to_frames()
                complex.index = complex_list[complex_i].index
                complex_list[complex_i] = complex
                molecule_list = list(complex.molecules)
                atom_counts[complex_i] = len(list(molecule_list[complex.current_frame].atoms))

                for atom in complex_list[complex_i].atoms:
                    atom._old_position = atom.position

                if complex_i > 0:
                    if (atom_counts[complex_i] > atom_counts[0]):
                        err_msg = "Error with receptor/ligand combination. Ligand cannot be larger than receptor."
                        Logs.error(err_msg)
                        self.send_notification(NotificationTypes.error, err_msg)
                        self.stop_scoring()
                        return
                    self._ligand_frames[complex.index] = complex.current_frame
                    self._ligand_atom_counts[complex.index] = atom_counts[complex_i]
                    for residue in complex_list[complex_i].residues:
                        residue.labeled = False
                    for atom in complex_list[complex_i].atoms:
                        atom.labeled = self.settings._labels

                    complex_updated_callback = functools.partial(self.complex_updated, complex_list[1:])
                    complex_list[complex_i].register_complex_updated_callback(complex_updated_callback)

            # self.unfreeze_button("Stop scoring")
            # self.menu.hide_scores(True)
            self._respond_to_update = False
            self.update_structures_deep(complex_list[0:], functools.partial(self.request_complexes, index_list[1:], self.setup_spheres))

        self._respond_to_update = True
        index_list = [receptor_index] + ligand_indices
        self.request_complexes(index_list, set_complexes)

    async def reassign_complexes_and_setup_streams(self, complexes):
        for complex in complexes:
            for atom in complex.atoms:
                atom._old_position = atom.position
        self._complexes = [self._complexes[0]] + complexes
        self.setup_spheres(complexes)

    async def get_last_n_complexes(self, n, complexes_shallow):
        complex_indices = [complex.index for complex in complexes_shallow[-n:]]
        complex_list = await self.request_complexes(complex_indices)
        self.reassign_complexes_and_setup_streams(complex_list)

    def get_updated_complexes(self):
        def update_complexes(complex_list):
            # update self._complexes positions from shallow list
            for complex in self._complexes:
                for thing in complex_list:
                    if thing.index == complex.index:
                        complex.position = thing.position
                        complex.rotation = thing.rotation
                        break
            self.prepare_complexes(self._complexes)
        self.request_complex_list(update_complexes)

    def clear_sphere_streams(self):
        streams_ready = self._streams_ready
        self._streams_ready = False

        scales = []
        if self._spheres != []:
            if streams_ready:
                # make all spheres transparent so user can't see deletion process
                for i in range(len(self._spheres)):
                    scales.append(0)
                try:
                    self._scale_stream.update(scales)
                except Exception:
                    Logs.error("Trying to update stream w/ incorrect size")
            Shape.destroy_multiple(self._spheres)
            self._spheres.clear()
        self._sphere_count = 0
        self._atom_count = 0
        self._color_stream = None
        self._scale_stream = None

    def setup_spheres(self, complex_list):
        if not self._is_running:
            return
        self._respond_to_update = True
        self._uploading_spheres = True
        self._delete_spheres = False
        self._curr_complex_list = complex_list

        for complex in complex_list:
            has_multiple_frames = len(list(complex.molecules)) > 1
            if has_multiple_frames:
                molecule_list = list(complex.molecules)
                curr_atoms = molecule_list[complex.current_frame].atoms
            else:
                curr_atoms = complex.atoms
            for atom in curr_atoms:
                self._atom_count += 1
                sphere = Sphere()
                sphere.color = nanome.util.Color(100, 100, 100, 120)
                sphere.radius = 1.3
                anchor = sphere.anchors[0]
                anchor.anchor_type = nanome.util.enums.ShapeAnchorType.Atom
                anchor.target = atom.index
                self._spheres.append(sphere)
        Shape.upload_multiple(self._spheres, self.on_shape_created)

    def on_shape_created(self, success):
        self._sphere_count = len(self._spheres)
        if self._sphere_count >= self._atom_count:
            self._uploading_spheres = False
            if self._delete_spheres:
                self._creating_streams = False
                self.clear_sphere_streams()
                if not self._stop_after_deleting_spheres:
                    self.get_full_complexes()
            else:
                self._creating_streams = True
                self.setup_streams(self._curr_complex_list)

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

    def setup_streams(self, complex_list):
        self._sphere_indices = []
        self._atom_indices = []
        for i in range(len(self._spheres)):
            self._sphere_indices.append(self._spheres[i].index)
        if self._color_stream is None or self._scale_stream is None:
            for complex in complex_list:
                for atom in complex.atoms:
                    self._atom_indices.append(atom.index)

            def on_stream_ready(complex_list):
                if self._color_stream is not None and self._scale_stream is not None and not (self.settings._labels ^ (self._label_stream is not None)):
                    self._streams_ready = True
                    self._to_display = True
                    self.get_updated_complexes()

            def on_label_stream_ready(stream, error):
                self._label_stream = stream
                on_stream_ready(complex_list)

            def on_color_stream_ready(stream, error):
                self._color_stream = stream
                on_stream_ready(complex_list)

            def on_scale_stream_ready(stream, error):
                self._scale_stream = stream
                on_stream_ready(complex_list)
            if self.settings._labels:
                self.create_writing_stream(self._atom_indices, nanome.api.streams.Stream.Type.label, on_label_stream_ready)
            self.create_writing_stream(self._sphere_indices, nanome.api.streams.Stream.Type.shape_color, on_color_stream_ready)
            self.create_writing_stream(self._sphere_indices, nanome.api.streams.Stream.Type.sphere_shape_radius, on_scale_stream_ready)
        else:
            self.get_updated_complexes()

    def prepare_complexes(self, complex_list):
        receptor = complex_list[0]
        for complex in complex_list:
            mat = complex.get_complex_to_workspace_matrix()
            for atom in complex.atoms:
                atom.position = mat * atom._old_position

        ligands = nanome.structure.Complex()
        self._ligands = ligands

        for complex in complex_list[1:]:
            complex = complex.convert_to_frames()
            for frame, molecule in enumerate(list(complex.molecules)):
                index = molecule.index
                if frame == complex.current_frame:
                    ligands.add_molecule(molecule)
                molecule.index = index

        receptor.io.to_pdb(self._protein_input.name, PDB_OPTIONS)
        ligands.io.to_sdf(self._ligands_input.name, SDF_OPTIONS)

    def nanobabel_convert(self, input_file, output_file):
        cmd = f'nanobabel convert -i {input_file} -o {output_file}'
        args = shlex.split(cmd)
        try:
            popen = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            popen.wait()
        except Exception:
            nanome.util.Logs.error("Couldn't execute nanobabel, please check if packet 'openbabel' is installed")
            return

    def display_results(self):
        if (not self._to_display):
            return
        self._ls_results.items = []

        # scores = []
        scores = []
        pcss = []
        with open(RESULTS_PATH) as results_file:
            results = results_file.readlines()
            number_of_lines = len(results)
            res_line_i = results.index('@RESULTS\n') + 4

            while res_line_i < number_of_lines:
                if results[res_line_i] == '\n':
                    break
                result = results[res_line_i].split('|')
                # name = result[1].strip()
                score = result[3].strip()
                pcs = result[5].strip()
                scores.append(score)
                pcss.append(pcs)
                res_line_i += 1

        for i, score in enumerate(scores):
            clone = self._pfb_result.clone()
            lbl = clone.get_content()
            if self.settings.show_total and self.settings.show_pcs:
                format_str = '{}: {}={}; {}={}'
                lbl.text_value = format_str.format(self._ligand_names[i], "TS", scores[i], "PCS", pcss[i])
            elif self.settings.show_total:
                format_str = '{}: {}={}'
                lbl.text_value = format_str.format(self._ligand_names[i], "TS", scores[i])
            elif self.settings.show_pcs:
                format_str = '{}: {}={}'
                lbl.text_value = format_str.format(self._ligand_names[i], "PCS", pcss[i])
            else:
                format_str = '{}'
                lbl.text_value = format_str.format(self._ligand_names[i])
            self._ls_results.items.append(clone)
        self.update_content(self._ls_results)


def main():
    description = "Display realtime scoring info about a selected ligand."
    plugin = nanome.Plugin("Realtime Scoring", description, "Scoring", True)
    plugin.set_plugin_class(RealtimeScoring)
    plugin.run()


if __name__ == "__main__":
    main()
