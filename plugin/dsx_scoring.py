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
