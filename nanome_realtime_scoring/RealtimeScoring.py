import nanome
import nanome.api.shapes as shapes
from nanome.util import Logs
from nanome.util.enums import NotificationTypes

import os
import shlex
import functools
import subprocess
import tempfile
import itertools
import stat
import time
from timeit import default_timer as timer

from .SettingsMenu import SettingsMenu

# TMP
from nanome._internal._structure._io._pdb.save import Options as PDBOptions
from nanome._internal._structure._io._sdf.save import Options as SDFOptions


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


class RealtimeScoring(nanome.PluginInstance):
    def benchmark_start(self, fn_name):
        if not fn_name in self._benchmarks:
            self._benchmarks[fn_name] = [0, 0, 0]
        self._benchmarks[fn_name][0] = timer()

    def benchmark_stop(self, fn_name):
        entry = self._benchmarks[fn_name]
        time = timer() - entry[0]
        entry[1] += time
        entry[2] += 1
        avg = entry[1] / entry[2]

        #nanome.util.Logs.debug('{:>10}    {:.2f}    (avg {:.2f})'.format(fn_name, time, avg))

    def start(self):
        self._benchmarks = {}

        self._protein_input = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self._ligands_input = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        self._ligands_converted = tempfile.NamedTemporaryFile(delete=False, suffix=".mol2")

        def button_pressed(button):
            if not self._is_button_loading:
                if self._is_running:
                    self.stop_scoring()
                else:
                    self.start_scoring()

        self._menu = nanome.ui.Menu.io.from_json(os.path.join(DIR, 'menu.json'))
        self.menu = self._menu
        self.settings = SettingsMenu(self, self.open_menu)

        self._p_selection = self._menu.root.find_node("Selection Panel", True)
        self._p_results = self._menu.root.find_node("Results Panel", True)
        self._ls_receptors = self._menu.root.find_node("Receptor List", True).get_content()
        self._ls_ligands = self._menu.root.find_node("Ligands List", True).get_content()
        self._ls_results = self._p_results.get_content()
        self._btn_score = self._menu.root.find_node("Button", True).get_content()
        self._btn_score.register_pressed_callback(button_pressed)

        self._pfb_complex = nanome.ui.LayoutNode()
        self._pfb_complex.add_new_button()

        self._pfb_result = nanome.ui.LayoutNode()
        self._pfb_result.add_new_label()

        self._is_running = False
        self._nanobabel_running = False
        self._dsx_running = False

        self._ligands = None
        self._receptor_index = None
        self._ligand_indices = []
        self._ligand_names = []
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

        self.menu.enabled = True
        self.update_menu(self.menu)

        self.request_complex_list(self.update_lists)

    def on_run(self):
        self.menu.enabled = True
        self.update_menu(self.menu)

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

    def dsx_start(self):
        dsx_path = os.path.join(DIR, 'dsx', 'dsx_linux_64.lnx')
        dsx_args = [dsx_path, '-P', self._protein_input.name, '-L', self._ligands_converted.name, '-D', 'pdb_pot_0511', '-pp', '-F', 'results.txt']
        try:
            self._dsx_process = subprocess.Popen(dsx_args, cwd=os.path.join(DIR, 'dsx'), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            self._dsx_running = True
            self.benchmark_start("dsx")
        except:
            nanome.util.Logs.error("Couldn't execute dsx, please check if executable is in the plugin folder and has permissions. Try executing chmod +x " + dsx_path)
            return

    def update(self):
        if not self._is_running:
            return

        if self._nanobabel_running:
            if self._nanobabel_process.poll() is not None:
                self.benchmark_stop("nanobabel")
                self._nanobabel_running = False
                self.dsx_start()
        elif self._dsx_running:
            if self._dsx_process.poll() is not None:
                self.benchmark_stop("dsx")
                self._dsx_running = False

                self.get_updated_complexes()
                dsx_output, _ = self._dsx_process.communicate()
                self.parse_scores(dsx_output)
                self.display_results()

    def on_complex_added(self):
        self.request_complex_list(self.update_lists)

    def on_complex_removed(self):
        self.request_complex_list(self.update_lists)

    def hide_scores(self, show_ligand_names=False):
        self._to_display = False
        self._ls_results.items = []
        if show_ligand_names:
            for i in range(len(self._ligand_names)):
                clone = self._pfb_result.clone()
                lbl = clone._get_content()
                lbl.text_value = '{}: {}'.format(self._ligand_names[i], "Loading scores")
                self._ls_results.items.append(clone)
            self.update_content(self._ls_results)
        else:
            clone = self._pfb_result.clone()
            lbl = clone._get_content()
            lbl.text_value = ''
            self._ls_results.items.append(clone)
            self.update_content(self._ls_results)

    def freeze_button(self):
        self._is_button_loading = True
        self._btn_score.text.value.set_all("Loading...")
        self.update_menu(self.menu)

    def unfreeze_button(self, text):
        self._is_button_loading = False
        self._btn_score.text.value.set_all(text)
        self.update_menu(self.menu)

    def start_scoring(self):
        self.freeze_button()

        if self._receptor_index is None:
            self.send_notification(NotificationTypes.error, "Please select a receptor")
            self.unfreeze_button("Start scoring")
            return

        if len(self._ligand_indices) == 0:
            self.send_notification(NotificationTypes.error, "Please select at least one ligand")
            self.unfreeze_button("Start scoring")
            return

        self.menu.title = "Scores"
        self.hide_scores()

        self._is_running = True
        self._nanobabel_running = False
        self._dsx_running = False

        self._p_selection.enabled = False
        self._p_results.enabled = True

        self.update_menu(self.menu)

        self._scores_ready = False
        self._creating_streams = False
        self._stop_after_deleting_spheres = False

        self.benchmark_start("total")

        self._respond_to_update = False
        self.get_full_complexes()

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
        self._p_selection.enabled = True
        self._p_results.enabled = False
        self._is_button_loading = False
        self._btn_score.text.value.set_all("Start scoring")
        self.request_complex_list(self.update_lists)
        self.update_menu(self.menu)

        self.unfreeze_button("Start scoring")

    def get_full_complexes(self):

        def set_complexes(complex_list):
            self._respond_to_update = True
            self._complexes = complex_list
            self._ligand_names = []
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
                    self._ligand_names.append(complex_list[complex_i].full_name)
                    for residue in complex_list[complex_i].residues:
                        residue.labeled = False
                    for atom in complex_list[complex_i].atoms:
                        atom.labeled = self.settings._labels

                    complex_updated_callback = functools.partial(self.complex_updated, complex_list[1:])
                    complex_list[complex_i].register_complex_updated_callback(complex_updated_callback)

            self.unfreeze_button("Stop scoring")
            self.hide_scores(True)
            self._respond_to_update = False
            self.update_structures_deep(complex_list[0:], functools.partial(self.request_complexes, index_list[1:], self.setup_spheres))

        self._respond_to_update = True
        index_list = [self._receptor_index] + self._ligand_indices
        self.request_complexes(index_list, set_complexes)

    def reassign_complexes_and_setup_streams(self, complexes):
        for complex in complexes:
            for atom in complex.atoms:
                atom._old_position = atom.position
        self._complexes = [self._complexes[0]] + complexes
        self.setup_spheres(complexes)

    def get_last_n_complexes(self, n, complexes_shallow):
        complex_indices = [complex.index for complex in complexes_shallow[-n:]]
        self.request_complexes(complex_indices, self.reassign_complexes_and_setup_streams)

    def get_updated_complexes(self):
        self.benchmark_stop("total")
        self.benchmark_start("total")

        def update_complexes(complex_list):
            # update self._complexes positions from shallow list
            for complex in self._complexes:
                for thing in complex_list:
                    if thing.index == complex.index:
                        complex.position = thing.position
                        complex.rotation = thing.rotation
                        break

            self.benchmark_stop("update")
            self.prepare_complexes(self._complexes)

        self.benchmark_start("update")
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
                except:
                    Logs.error("Trying to update stream w/ incorrect size")
            # destroy spheres
            for i in range(len(self._spheres)):
                self._spheres[i].destroy()
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
                sphere = shapes.Sphere()
                sphere.color = nanome.util.Color(100, 100, 100, 120)
                sphere.radius = 1.3
                anchor = sphere.anchors[0]
                anchor.anchor_type = nanome.util.enums.ShapeAnchorType.Atom
                anchor.target = atom.index
                self._spheres.append(sphere)

        shapes.Sphere.upload_multiple(self._spheres)
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
                self.hide_scores(True)
            elif self._creating_streams:
                self._respond_to_update = False
                self._creating_streams = False
                self.clear_sphere_streams()
                self.hide_scores(True)
                self.get_full_complexes()
            else:
                self._respond_to_update = True

    def setup_streams(self, complex_list):
        self._sphere_indices = []
        self._atom_indices = []
        for i in range(len(self._spheres)):
            self._sphere_indices.append(self._spheres[i].index)
        if self._color_stream == None or self._scale_stream == None:
            indices = []
            for complex in complex_list:
                for atom in complex.atoms:
                    self._atom_indices.append(atom.index)

            def on_stream_ready(complex_list):
                if self._color_stream != None and self._scale_stream != None and not (self.settings._labels ^ (self._label_stream != None)):
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

        self.nanobabel_start()

    def nanobabel_start(self):
        cmd = f'nanobabel convert -i {self._ligands_input.name} -o {self._ligands_converted.name}'
        args = shlex.split(cmd)
        try:
            self._nanobabel_process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self._nanobabel_running = True
            self.benchmark_start("nanobabel")
        except:
            nanome.util.Logs.error("Couldn't execute nanobabel, please check if packet 'openbabel' is installed")
            return

    def parse_scores(self, dsx_output):
        lines = dsx_output.splitlines()
        number_of_lines = len(lines)

        line_index = 0

        def find_next_ligand():
            nonlocal line_index
            while line_index < number_of_lines - 1:
                if lines[line_index].startswith("# Receptor-Ligand:"):
                    line_index += 1
                    return True
                line_index += 1
            return False

        if not find_next_ligand():
            Logs.error("Couldn't parse DSX scores")
            Logs.error("Output:\n" + str(dsx_output))
            err_msg = "Error parsing scores. Are the ligand and receptor selected correctly?"
            self.send_notification(NotificationTypes.error, err_msg)
            self.stop_scoring()
            self.clear_sphere_streams()
            return

        ligand_index = 0
        has_next_ligand = True
        while has_next_ligand:
            scores = dict()
            last_tuple = None
            last_arr = None
            score_min = None
            score_max = None

            while line_index < number_of_lines - 1:
                line = lines[line_index]
                if line.startswith("# End of pair potentials"):
                    has_next_ligand = find_next_ligand()
                    break
                line_items = line.split("__")
                atom_items = line_items[1].split("_")
                score = float(line_items[2])
                ireceptor_iligand = (int(atom_items[1]), int(atom_items[2]))

                if last_tuple != ireceptor_iligand:
                    if ireceptor_iligand in scores:
                        last_arr = scores[ireceptor_iligand]
                    else:
                        last_arr = []
                        scores[ireceptor_iligand] = last_arr
                last_tuple = ireceptor_iligand
                last_arr.append(score)

                if score_min == None or score < score_min:
                    score_min = score
                if score_max == None or score > score_max:
                    score_max = score
                line_index += 1

            if score_max is None or score_min is None:
                continue
            score_gap = max(score_max - score_min, 0.01)
            try:
                for atom_tuple, score_arr in scores.items():
                    score = sum(score_arr) / len(score_arr)
                    bfactor = ((score - score_min) / score_gap) * BFACTOR_GAP + BFACTOR_MIN
                    bfactor_two = score / (-score_min if score < 0 else score_max)
                    molecule = self._ligands._molecules[atom_tuple[0] - 1 + ligand_index]
                    if not hasattr(molecule, "atom_score_limits"):
                        molecule.atom_score_limits = [float('inf'), float('-inf')]
                    if score < molecule.atom_score_limits[0]:
                        molecule.atom_score_limits[0] = score
                    elif score > molecule.atom_score_limits[1]:
                        molecule.atom_score_limits[1] = score

                    atom = next(itertools.islice(molecule.atoms, atom_tuple[1] - 1, atom_tuple[1]))
                    atom.score = score
            except:
                err_msg = "Error parsing ligand scores. Are your ligands missing bonds?"
                self.send_notification(NotificationTypes.error, err_msg)
                nanome.util.Logs.error(err_msg)
                self.stop_scoring()
                return

            ligand_index += 1

        labels = []
        colors = []
        scales = []
        for atom in self._ligands.atoms:
            red = 0
            green = 0
            blue = 0
            alpha = 141
            if hasattr(atom, "score"):
                ligand = atom.molecule
                denominator = -ligand.atom_score_limits[0] if atom.score < 0 else ligand.atom_score_limits[1]
                norm_score = atom.score / denominator
                red = 255 if norm_score > 0 else 0
                green = 255 if -norm_score >= 0 else 0
                red_scale = int(max(norm_score * 255, 0))
                green_scale = int(max(-norm_score * 255, 0))
                scale = max(green_scale, red_scale) / 255. + 1
            else:
                norm_score = 0.0
                green = 0
                red = 255
                scale = 1

            labels.append(f'{norm_score:.3f}')
            colors.append(red)
            colors.append(green)
            colors.append(blue)
            colors.append(alpha)
            scales.append(scale)

        try:
            if self._streams_ready:
                if self._label_stream:
                    self._label_stream.update(labels)
                self._color_stream.update(colors)
                self._scale_stream.update(scales)
        except:
            Logs.error("Error while updating sphere stream")

    def display_results(self):
        if (not self._to_display):
            return
        self._ls_results.items = []

        #scores = []
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
                name = result[1].strip()
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

    def update_lists(self, complex_list):
        if self._is_running:
            return

        self._shallow_complexes = complex_list

        def update_selected_ligands():
            self._ligand_indices = []
            for item in self._ls_ligands.items:
                btn = item.get_content()
                if btn.selected:
                    self._ligand_indices.append(btn.index)

        def receptor_pressed(receptor):
            for item in self._ls_receptors.items:
                item.get_content().selected = False

            receptor.selected = True
            self._receptor_index = receptor.index

            for item in self._ls_ligands.items:
                ligand = item.get_content()
                ligand.unusable = receptor.index == ligand.index
                if ligand.selected and ligand.unusable:
                    ligand.selected = False
            update_selected_ligands()

            self.update_menu(self.menu)

        def ligand_pressed(ligand):
            ligand.selected = not ligand.selected
            self.update_content(ligand)
            update_selected_ligands()

        def populate_list(ls, cb):
            ls.items = []
            for complex in complex_list:
                clone = self._pfb_complex.clone()
                btn = clone.get_content()
                btn.text.value.set_all(complex.full_name)
                btn.index = complex.index
                btn.register_pressed_callback(cb)
                ls.items.append(clone)
            self.update_content(ls)

        self._receptor_index = None
        self._ligand_indices = []
        populate_list(self._ls_receptors, receptor_pressed)
        populate_list(self._ls_ligands, ligand_pressed)


def main():
    description = "Display realtime scoring info about a selected ligand."
    plugin = nanome.Plugin("Realtime Scoring", description, "Scoring", True)
    plugin.set_plugin_class(RealtimeScoring)
    plugin.run()


if __name__ == "__main__":
    main()
