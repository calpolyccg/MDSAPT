import pytest

import os

import MDAnalysis as mda
from MDAnalysis.converters.RDKit import atomgroup_to_mol
from MDAnalysis.topology.guessers import guess_types

from rdkit.Chem.rdmolops import AddHs
import rdkit.Chem as chem

from ..reader import InputReader
from ..sapt import TrajectorySAPT
from ..optimizer import Optimizer


class TestSAPT(object):

    def setup(self):
        self.settings = InputReader(os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        self.opt = Optimizer(self.settings)
        self.Unv = mda.Universe(self.settings.top_path, self.settings.trj_path)
        elements = guess_types(self.Unv.atoms.names)
        self.Unv.add_TopologyAttr('elements', elements)

    def test_run_sapt(self):
        SAPT12 = TrajectorySAPT(self.settings, self.opt)
        SAPT12.run(self.settings.trj_settings['start'], self.settings.trj_settings['stop'], self.settings.trj_settings['step'])
        cols_act = SAPT12.results.columns
        cols_exp = ['residues', 'time', 'total', 'electrostatic', 'exchange', 'induction', 'dispersion']
        assert (len(cols_act) == len(cols_exp))
        assert all([cols_act[i] == cols_exp[i] for i in range(len(cols_act))])
        SAPT12.results.to_csv('sapt_test.csv')
