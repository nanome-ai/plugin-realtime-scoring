import asyncio
import nanome
import os
import unittest
from nanome.api import structure, PluginInstance
from plugin.RealtimeScoring import RealtimeScoring
from random import randint
from unittest.mock import MagicMock


assets_dir = os.path.join(os.path.dirname(__file__), 'assets')


def run_awaitable(awaitable, *args, **kwargs):
    loop = asyncio.get_event_loop()
    if loop.is_running:
        loop = asyncio.new_event_loop()
    result = loop.run_until_complete(awaitable(*args, **kwargs))
    loop.close()
    return result


class RealtimeScoringTestCase(unittest.TestCase):

    def setUp(self):
        self.plugin = RealtimeScoring()
        PluginInstance._instance = self.plugin
        self.plugin._network = MagicMock()
        nanome._internal.network.plugin_network.PluginNetwork._instance = MagicMock()

        self.receptor_pdb = os.path.join(assets_dir, '5ceo.pdb')
        self.receptor_comp = structure.Complex.io.from_pdb(path=self.receptor_pdb)
        self.ligand_pdb = os.path.join(assets_dir, '50D.pdb')
        self.ligand_comp = structure.Complex.io.from_pdb(path=self.ligand_pdb)

    def test_start_scoring(self):
        async def validate_start_scoring(self):
            receptor_index = self.receptor_comp.index
            ligand_indices = [self.ligand_comp.index]
            request_complexes_fut = asyncio.Future()
            request_complexes_fut.set_result([self.receptor_comp, self.ligand_comp])
            self.plugin.request_complexes = MagicMock(return_value=request_complexes_fut)
            await self.plugin.score_ligand(receptor_index, ligand_indices)
        run_awaitable(validate_start_scoring, self)