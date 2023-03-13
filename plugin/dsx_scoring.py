import itertools
import os
import subprocess
import shlex
import tempfile
from nanome.api import structure
from nanome.util import Logs


__all__ = ['score_ligands']


DIR = os.path.dirname(__file__)

SDF_OPTIONS = structure.Complex.io.SDFSaveOptions()
SDF_OPTIONS.write_bonds = True
PDB_OPTIONS = structure.Complex.io.PDBSaveOptions()
PDB_OPTIONS.write_bonds = True


def score_ligands(receptor: structure.Complex, ligand_comps: list[structure.Complex]):
    receptor_pdb = tempfile.NamedTemporaryFile(suffix='.pdb')
    receptor.io.to_pdb(receptor_pdb.name, PDB_OPTIONS)
    # For each ligand, generate a PDB file and run DSX
    for ligand_comp in ligand_comps:
        ligand_sdf = tempfile.NamedTemporaryFile(suffix='.sdf')
        ligand_comp.io.to_sdf(ligand_sdf.name, SDF_OPTIONS)
        # Convert ligand sdf to a mol2.
        ligand_mol2 = tempfile.NamedTemporaryFile(suffix='.mol2')
        # Convert ligand from sdf to mol2.
        nanobabel_convert(ligand_sdf.name, ligand_mol2.name)
        # Run DSX and retreive data from the subprocess.
        dsx_output_file = tempfile.NamedTemporaryFile(suffix='.txt')
        dsx_output = run_dsx(receptor_pdb.name, ligand_mol2.name, dsx_output_file.name)
        atom_scores = dsx_parse(dsx_output, ligand_comp)
        return atom_scores


def run_dsx(receptor_pdb, ligands_mol2, output_file_path):
    "Run DSX and write output to provided output_file."
    dsx_path = os.path.join(DIR, 'dsx', 'dsx_linux_64.lnx')
    dsx_args = [
        dsx_path, '-P', receptor_pdb, '-L', ligands_mol2, '-D', 'pdb_pot_0511',
        '-pp', '-F', output_file_path
    ]
    try:
        dsx_process = subprocess.Popen(dsx_args, cwd=os.path.join(DIR, 'dsx'), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    except Exception:
        Logs.error("Couldn't execute dsx, please check if executable is in the plugin folder and has permissions. Try executing chmod +x " + dsx_path)
        return
    dsx_process.wait()
    dsx_output, _ = dsx_process.communicate()
    return dsx_output


def nanobabel_convert(input_file, output_file):
    cmd = f'nanobabel convert -i {input_file} -o {output_file}'
    args = shlex.split(cmd)
    try:
        popen = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        popen.wait()
    except Exception:
        Logs.error("Couldn't execute nanobabel, please check if packet 'openbabel' is installed")
        return


def dsx_parse(dsx_output, ligand_comp):
    ligand_sets = dsx_output.split("# Receptor-Ligand:")[1:]
    scores = dict()
    atom_scores = list()
    for ligand_set in ligand_sets:
        for line in ligand_set.splitlines():
            if not line or line.startswith("#"):
                continue
            line_items = line.split("__")
            try:
                atom_items = line_items[1].split("_")
            except IndexError:
                print("Error parsing line: ", line)
                continue
            score = float(line_items[2])
            ireceptor_iligand = (int(atom_items[1]), int(atom_items[2]))
            if ireceptor_iligand not in atom_scores:
                scores[ireceptor_iligand] = [score]
            else:
                scores[ireceptor_iligand].append(score)
        if not scores:
            continue
        # Attach calculated score to each atom
        ligand_molecule = next(
            mol for i, mol in enumerate(ligand_comp.molecules)
            if i == ligand_comp.current_frame
        )
        if not hasattr(ligand_molecule, "atom_score_limits"):
            ligand_molecule.atom_score_limits = [float('inf'), float('-inf')]

        for atom_tuple, score_arr in scores.items():
            score = sum(score_arr) / len(score_arr)
            ligand_molecule.atom_score_limits[0] = min(ligand_molecule.atom_score_limits[0], score)
            ligand_molecule.atom_score_limits[1] = max(ligand_molecule.atom_score_limits[1], score)
            atom = next(itertools.islice(ligand_molecule.atoms, atom_tuple[1] - 1, atom_tuple[1]))
            atom_scores.append((atom, score))
    return atom_scores


def parse_ligand_lines(dsx_output_lines):
    ligand_lines = []
    receptor_ligand_line = '# Receptor-Ligand'
    while True:
        try:
            line = dsx_output_lines.pop(0)
        except IndexError:
            break
        if line.startswith(receptor_ligand_line):
            continue
        elif line.startswith(endline):
            # Pop empty line with '----------\n'
            try:
                ligand_lines.pop(-1)
            except IndexError:
                pass
            break
        else:
            ligand_lines.append(line)
    return ligand_lines


def find_next_ligand(line_list):
    ligand_line_start = "# Receptor-Ligand:"
    starting_line_index = 0
    ligand_found = False
    for i in range(0, len(line_list)):
        line = line_list[i]
        if line.startswith(ligand_line_start):
            starting_line_index = i
            ligand_found = True
            break
    del line_list[:starting_line_index]
    return ligand_found


def parse_scores(self, dsx_output):
    lines = dsx_output.splitlines()
    number_of_lines = len(lines)
    line_index = 0
    if not find_next_ligand():
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

            if score_min is None or score < score_min:
                score_min = score
            if score_max is None or score > score_max:
                score_max = score
            line_index += 1

        if score_max is None or score_min is None:
            continue
        # score_gap = max(score_max - score_min, 0.01)
        try:
            for atom_tuple, score_arr in scores.items():
                score = sum(score_arr) / len(score_arr)
                # bfactor = ((score - score_min) / score_gap) * BFACTOR_GAP + BFACTOR_MIN
                # bfactor_two = score / (-score_min if score < 0 else score_max)
                molecule = self._ligands._molecules[atom_tuple[0] - 1 + ligand_index]
                if not hasattr(molecule, "atom_score_limits"):
                    molecule.atom_score_limits = [float('inf'), float('-inf')]
                if score < molecule.atom_score_limits[0]:
                    molecule.atom_score_limits[0] = score
                elif score > molecule.atom_score_limits[1]:
                    molecule.atom_score_limits[1] = score

                atom = next(itertools.islice(molecule.atoms, atom_tuple[1] - 1, atom_tuple[1]))
                atom.score = score
        except Exception:
            # err_msg = "Error parsing ligand scores. Are your ligands missing bonds?"
            # self.send_notification(NotificationTypes.error, err_msg)
            # nanome.util.Logs.error(err_msg)
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
    except Exception:
        print("Error while updating sphere stream")
