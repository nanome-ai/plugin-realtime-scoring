import asyncio
import itertools
import os
import unittest

from nanome.api import structure
from dsx import scoring_algo
from random import randint


assets_dir = os.path.join(os.path.dirname(__file__), 'assets')


def run_awaitable(awaitable, *args, **kwargs):
    loop = asyncio.get_event_loop()
    if loop.is_running:
        loop = asyncio.new_event_loop()
    result = loop.run_until_complete(awaitable(*args, **kwargs))
    loop.close()
    return result


class DsxTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.receptor_pdb = os.path.join(assets_dir, '5ceo.pdb')
        cls.receptor_comp = structure.Complex.io.from_pdb(path=cls.receptor_pdb)
        cls.ligand_pdb = os.path.join(assets_dir, '50D.pdb')
        cls.ligand_comp = structure.Complex.io.from_pdb(path=cls.ligand_pdb)
        # Generate indices for receptor and ligand
        for residue in itertools.chain(cls.receptor_comp.residues, cls.ligand_comp.residues):
            residue.index = randint(1000000000, 9999999999)
            for atom in residue.atoms:
                atom.index = randint(1000000000, 9999999999)

    def test_parse_output(self):
        results_file = os.path.join(assets_dir, 'dsx_output.txt')
        with open(results_file, 'r') as f:
            dsx_output = f.read()
        expected_atoms_with_scores = 29
        atom_scores = scoring_algo.parse_output(dsx_output, self.ligand_comp)
        self.assertEqual(len(atom_scores), expected_atoms_with_scores)

    def test_parse_results(self):
        results_file = os.path.join(assets_dir, 'dsx_results.txt')
        expected_total_score = -127.995
        expected_per_contact_score = -0.2
        aggregate_scores = scoring_algo.parse_results(results_file)
        self.assertEqual(aggregate_scores[0]['total_score'], expected_total_score)
        self.assertEqual(aggregate_scores[0]['per_contact_score'], expected_per_contact_score)
