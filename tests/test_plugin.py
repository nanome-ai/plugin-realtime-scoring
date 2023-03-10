import asyncio
import nanome
import os
import unittest
import itertools
from nanome.api import structure, PluginInstance, shapes
from plugin.RealtimeScoring import RealtimeScoring
from random import randint
from unittest.mock import MagicMock
from plugin.dsx_parser import dsx_parse
from random import randint


assets_dir = os.path.join(os.path.dirname(__file__), 'assets')


def run_awaitable(awaitable, *args, **kwargs):
    loop = asyncio.get_event_loop()
    if loop.is_running:
        loop = asyncio.new_event_loop()
    result = loop.run_until_complete(awaitable(*args, **kwargs))
    loop.close()
    return result


class RealtimeScoringTestCase(unittest.TestCase):

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

    def setUp(self):
        self.plugin = RealtimeScoring()
        PluginInstance._instance = self.plugin
        self.plugin._network = MagicMock()
        nanome._internal.network.plugin_network.PluginNetwork._instance = MagicMock()

    def tearDown(self):
        self.plugin.temp_dir.cleanup()

    def test_score_ligand(self):
        async def validate_score_ligand(self):
            receptor_index = self.receptor_comp.index
            ligand_indices = [self.ligand_comp.index]

            # Mock sShapes upload_multiple call
            upload_multiple_fut = asyncio.Future()
            upload_multiple_fut.set_result(None)
            shapes.Shape.upload_multiple = MagicMock(return_value=upload_multiple_fut)
            # Mock request complexes call
            request_complexes_fut = asyncio.Future()
            request_complexes_fut.set_result([self.receptor_comp, self.ligand_comp])
            self.plugin.request_complexes = MagicMock(return_value=request_complexes_fut)

            # Mock create stream call
            create_writing_stream_fut = asyncio.Future()
            stream_mock = MagicMock()
            create_writing_stream_fut.set_result((stream_mock, MagicMock()))
            self.plugin.create_writing_stream = MagicMock(return_value=create_writing_stream_fut)
            # Run score ligand
            await self.plugin.score_ligand(receptor_index, ligand_indices)
            # Assert update called on stream
            stream_mock.update.assert_called_once()
        run_awaitable(validate_score_ligand, self)

    def test_dsx_parser(self):
        results_file = os.path.join(assets_dir, 'dsx_output.txt')
        for atom in self.ligand_comp.atoms:
            self.assertTrue(not hasattr(atom, 'score'))
        dsx_parse(results_file, self.ligand_comp)
        expected_atoms_with_scores = 28
        score_found = 0
        for atom in self.ligand_comp.atoms:
            if hasattr(atom, 'score'):
                score_found += 1
        self.assertEqual(score_found, expected_atoms_with_scores)
