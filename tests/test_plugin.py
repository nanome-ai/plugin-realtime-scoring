import asyncio
import nanome
from nanome.util import enums
import os
import unittest
import itertools
from nanome.api import structure, PluginInstance, shapes
from plugin.RealtimeScoring import RealtimeScoring
from random import randint
from unittest.mock import MagicMock
from plugin.dsx_scoring import dsx_parse_output
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

    def test_setup_receptor_and_ligands(self):
        async def validate_setup_receptor_and_ligands(self):
            receptor_index = self.receptor_comp.index
            ligand_indices = [self.ligand_comp.index]

            # Mock Shapes upload_multiple call
            upload_multiple_fut = asyncio.Future()
            upload_multiple_fut.set_result(None)
            shapes.Shape.upload_multiple = MagicMock(return_value=upload_multiple_fut)
            
            # Mock request complexes call
            request_complexes_fut = asyncio.Future()
            request_complexes_fut.set_result([self.receptor_comp, self.ligand_comp])
            self.plugin.request_complexes = MagicMock(return_value=request_complexes_fut)

            # Mock create stream calls
            create_stream_fut = asyncio.Future()
            stream_mock = MagicMock()
            create_stream_fut.set_result((stream_mock, unittest.mock.ANY))
            self.plugin.create_writing_stream = MagicMock(return_value=create_stream_fut)
            self.assertEqual(self.plugin.create_writing_stream.call_count, 0)
            # Run function.
            # Assert mocks were called
            await self.plugin.setup_receptor_and_ligands(receptor_index, ligand_indices)
            self.assertEqual(self.plugin.request_complexes.call_count, 1)
            self.assertEqual(shapes.Shape.upload_multiple.call_count, 1)
            self.assertEqual(self.plugin.create_writing_stream.call_count, 3)
            self.assertEqual(self.plugin.receptor_comp, self.receptor_comp)
            self.assertEqual(self.plugin.ligand_comps, [self.ligand_comp])
        run_awaitable(validate_setup_receptor_and_ligands, self)

    def test_score_ligands(self):
        async def validate_score_ligands(self):
            self.plugin.receptor_comp = self.receptor_comp
            self.plugin.ligand_comps = [self.ligand_comp]
            self.plugin.color_stream = MagicMock()
            self.plugin.size_stream = MagicMock()
            self.plugin.label_stream = MagicMock()
            # Run function
            self.plugin.settings.update_labels = True
            await self.plugin.score_ligands()
            # Assert stream updates were called
            self.assertEqual(self.plugin.color_stream.update.call_count, 1)
            self.assertEqual(self.plugin.size_stream.update.call_count, 1)
            self.assertEqual(self.plugin.label_stream.update.call_count, 1)
            # Call again with labels disabled, and make sure update was not called
            self.plugin.update_content = MagicMock()
            self.plugin.settings.update_labels = False
            await self.plugin.score_ligands()
            self.assertEqual(self.plugin.color_stream.update.call_count, 2)
            self.assertEqual(self.plugin.size_stream.update.call_count, 2)
            self.assertEqual(self.plugin.label_stream.update.call_count, 1)
        run_awaitable(validate_score_ligands, self)

    def test_dsx_parse_output(self):
        results_file = os.path.join(assets_dir, 'dsx_output.txt')
        with open(results_file, 'r') as f:
            dsx_output = f.read()
        expected_atoms_with_scores = 29
        atom_scores = dsx_parse_output(dsx_output, self.ligand_comp)
        self.assertEqual(len(atom_scores), expected_atoms_with_scores)
