import asyncio
import itertools
import nanome
import os
import unittest
from unittest.mock import MagicMock

from nanome.api import structure, PluginInstance, shapes
from dsx.scoring_algo import parse_output
from plugin.RealtimeScoring import RealtimeScoring
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
        cls.receptor_pdb = os.path.join(assets_dir, '5ceo_protein.pdb')
        cls.receptor_comp = structure.Complex.io.from_pdb(path=cls.receptor_pdb)
        cls.ligand_pdb = os.path.join(assets_dir, '50D.pdb')
        cls.ligand_comp = structure.Complex.io.from_pdb(path=cls.ligand_pdb)
        cls.generate_random_indices(cls.receptor_comp)
        cls.generate_random_indices(cls.ligand_comp)
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

    def test_setup_receptor_and_ligands(self):
        async def validate_setup_receptor_and_ligands(self):
            receptor_index = self.receptor_comp.index
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
            ligand_residues = list(self.ligand_comp.residues)
            # Run function.
            await self.plugin.setup_receptor_and_ligands(receptor_index, ligand_residues)
            # Assert mocks were called
            self.assertEqual(self.plugin.request_complexes.call_count, 1)
            self.assertEqual(shapes.Shape.upload_multiple.call_count, 1)
            self.assertEqual(self.plugin.create_writing_stream.call_count, 3)
            # Assert receptor and ligand were set
            self.assertEqual(self.plugin.receptor_comp, self.receptor_comp)
            self.assertEqual(self.plugin.ligand_residues, ligand_residues)
        run_awaitable(validate_setup_receptor_and_ligands, self)

    def test_score_ligands(self):
        async def validate_score_ligands(self):
            self.plugin.complex_cache = [self.receptor_comp, self.ligand_comp]
            self.plugin.receptor_index = self.receptor_comp.index
            self.plugin.ligand_residues = list(self.ligand_comp.residues)
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
            # Call again with labels disabled, and make sure label update was not called
            self.plugin.update_content = MagicMock()
            self.plugin.settings.update_labels = False
            await self.plugin.score_ligands()
            self.assertEqual(self.plugin.color_stream.update.call_count, 2)
            self.assertEqual(self.plugin.size_stream.update.call_count, 2)
            self.assertEqual(self.plugin.label_stream.update.call_count, 1)
        run_awaitable(validate_score_ligands, self)
    
    def test_score_ligand_one_complex(self):
        """Validate score ligand when ligand and receptor are same Complex."""
        async def validate_score_ligands_one_complex(self):
            pdb_path = os.path.join(assets_dir, '5ceo.pdb')
            comp = structure.Complex.io.from_pdb(path=pdb_path)
            self.generate_random_indices(comp)
            ligand_residues = [res for res in comp.residues if res.name == '50D']
            self.assertEqual(len(ligand_residues), 1)

            # Set up plugin state.
            self.plugin.complex_cache = [self.receptor_comp]
            self.plugin.receptor_index = self.receptor_comp.index
            self.plugin.ligand_residues = ligand_residues
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

            # Call again with labels disabled, and make sure label update was not called
            self.plugin.update_content = MagicMock()
            self.plugin.settings.update_labels = False
            await self.plugin.score_ligands()
            self.assertEqual(self.plugin.color_stream.update.call_count, 2)
            self.assertEqual(self.plugin.size_stream.update.call_count, 2)
            self.assertEqual(self.plugin.label_stream.update.call_count, 1)
        run_awaitable(validate_score_ligands_one_complex, self)

    @staticmethod
    def generate_random_indices(comp):
        min_index = 1000000000
        max_index = 9999999999
        comp.index = randint(min_index, max_index)
        for residue in comp.residues:
            residue.index = randint(min_index, max_index)
            for atom in residue.atoms:
                atom.index = randint(min_index, max_index)