import pytest

import os

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_types

from ..config import Config, load_from_yaml_file
from ..sapt import TrajectorySAPT


class TestSAPT:

    def setup(self):
        self.settings: Config = load_from_yaml_file(
            os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        self.unv = self.settings.analysis.create_universe()
        elements = guess_types(self.unv.atoms.names)
        self.unv.add_TopologyAttr('elements', elements)

    def test_run_sapt(self):
        SAPT12 = TrajectorySAPT(self.settings)
        SAPT12.run(1, 2)
        cols_act = SAPT12.results.columns
        cols_exp = ['residues', 'time', 'total', 'electrostatic', 'exchange', 'induction', 'dispersion']
        assert (len(cols_act) == len(cols_exp))
        assert all([cols_act[i] == cols_exp[i] for i in range(len(cols_act))])
        SAPT12.results.to_csv('sapt_test.csv')
