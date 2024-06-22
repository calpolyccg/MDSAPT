"""
Tests for the SAPT analysis objects
"""
import os

from MDAnalysis.topology.guessers import guess_types

from ..config import load_from_yaml_file
from ..sapt import TrajectorySAPT, DockingSAPT

traj_settings = load_from_yaml_file(
    os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))

dock_settings = load_from_yaml_file(
    os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'docking_in.yaml'))

unv = traj_settings.analysis.create_universe()
elements = guess_types(unv.atoms.names)
unv.add_TopologyAttr('elements', elements)


class TestSAPT:
    """
    Test object for SAPT analysis
    """

    def test_run_traj_sapt(self) -> None:
        """
        Test running trajectory SAPT
        """
        sapt12 = TrajectorySAPT(traj_settings)
        sapt12.run(1, 2)
        cols_act = sapt12.results.columns
        cols_exp = ['residues', 'time', 'total', 'electrostatic',
                    'exchange', 'induction', 'dispersion']
        assert (len(cols_act) == len(cols_exp))
        assert all((cols_act[i] == cols_exp[i] for i in range(len(cols_act))))
        sapt12.results.to_csv('sapt_test.csv')

    def test_run_dock_sapt(self) -> None:
        """
        Test running docking SAPT
        """
        docking = DockingSAPT(dock_settings)
        docking.run()
        cols_act = docking.results.columns
        cols_exp = ['structure', 'pair', 'total', 'electrostatic',
                    'exchange', 'induction', 'dispersion']
        assert (all((cols_act[i] == cols_exp[i] for i in range(len(cols_act)))))
