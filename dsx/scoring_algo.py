import itertools
import io
import os
import tempfile
from nanome.api import structure
from nanome.util import Logs, Process

__all__ = ['score_ligands']


DIR = os.path.dirname(__file__)

SDF_OPTIONS = structure.Complex.io.SDFSaveOptions()
PDB_OPTIONS = structure.Complex.io.PDBSaveOptions()
PDB_OPTIONS.write_bonds = True


async def score_ligands(receptor: structure.Complex, ligand_comps: 'list[structure.Complex]'):
    output = []
    with tempfile.TemporaryDirectory() as dir:
        receptor_pdb = tempfile.NamedTemporaryFile(dir=dir, delete=False, suffix='.pdb')
        receptor.io.to_pdb(receptor_pdb.name, PDB_OPTIONS)
        # For each ligand, generate a PDB file and run DSX
        for ligand_comp in ligand_comps:
            ligand_sdf = tempfile.NamedTemporaryFile(dir=dir, delete=False, suffix='.sdf')
            ligand_comp.io.to_sdf(ligand_sdf.name, SDF_OPTIONS)
            # Convert ligand sdf to a mol2.
            ligand_mol2 = tempfile.NamedTemporaryFile(dir=dir, delete=False, suffix='.mol2')
            # Convert ligand from sdf to mol2.
            await nanobabel_convert(ligand_sdf.name, ligand_mol2.name)
            # Run DSX and retreive data from the subprocess.
            dsx_results_file = tempfile.NamedTemporaryFile(dir=dir, delete=False, suffix='.txt')
            dsx_output = await run_dsx(receptor_pdb.name, ligand_mol2.name, dsx_results_file.name)
            atom_scores = parse_output(dsx_output, ligand_comp)
            aggregate_scores = parse_results(dsx_results_file.name)
            ligand_data = {
                'complex_index': ligand_comp.index,
                'aggregate_scores': aggregate_scores,
                'atom_scores': atom_scores
            }
            output.append(ligand_data)
    return output


async def run_dsx(receptor_pdb, ligands_mol2, output_file_path) -> str:
    """Run DSX and write output to provided output_file."""
    dsx_path = os.path.join(DIR, 'bin', 'dsx_linux_64.lnx')
    pdb_pot_0511 = os.path.join(DIR, 'bin', 'pdb_pot_0511')
    dsx_args = [
        dsx_path, '-P', receptor_pdb, '-L', ligands_mol2, '-D', pdb_pot_0511,
        '-pp', '-F', output_file_path
    ]
    dsx_stdout = io.StringIO()
    try:
        dsx_process = Process(dsx_path, dsx_args, label="DSX", output_text=True)
        dsx_process.on_output = lambda output: dsx_stdout.write(output + '\n')
        await dsx_process.start()
    except Exception:
        Logs.error("Couldn't execute dsx, please check if executable is in the plugin folder and has permissions. Try executing chmod +x " + dsx_path)
        return
    return dsx_stdout.getvalue()


async def nanobabel_convert(input_file, output_file):
    nanobabel_path = 'nanobabel'
    cmd_args = ['convert', '-i', input_file, '-o', output_file]
    nanobabel_process = Process(nanobabel_path, cmd_args, label="nanobabel", output_text=True)
    await nanobabel_process.start()


def parse_output(dsx_output, ligand_comp):
    """Get per atom scores from output of DSX process."""
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

        for atom_tuple, score_arr in scores.items():
            score = sum(score_arr) / len(score_arr)
            atom = next(itertools.islice(ligand_molecule.atoms, atom_tuple[1] - 1, atom_tuple[1]))
            atom_scores.append((atom.index, score))
    return atom_scores


def parse_results(dsx_output_file):
    """Parse the output of DSX and return a list of total scores and per contact scores."""
    data = []
    with open(dsx_output_file) as results_file:
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
            data.append(
                {
                    'total_score': float(score),
                    'per_contact_score': float(pcs)
                }
            )
            res_line_i += 1
    return data
