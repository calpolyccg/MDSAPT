"""
Tests for Ensemble framework
"""
from pathlib import Path

from ..utils.ensemble import Ensemble

resources_dir = Path(__file__).parent / 'testing_resources'

working_docking_settings = dict(
    type='docking',
    protein=f'{resources_dir}/docking_sep_test/2hnt.pdb',
    ligands=[f'{resources_dir}/docking_sep_test/ligands/15U0.pdb',
             f'{resources_dir}/docking_sep_test/ligands/15U1.pdb',
             f'{resources_dir}/docking_sep_test/ligands/98P_2.pdb'],
    pairs=[(3, 14)]
)


def test_list_of_merged() -> None:
    """
    Test loading merged
    """
    merged_list = [f'{resources_dir}/docking_merged_test/2hnt_15U0.pdb',
                   f'{resources_dir}/docking_merged_test/2hnt_98P.pdb']

    ens: Ensemble = Ensemble.build_from_files(merged_list)
    assert all(k in merged_list for k in ens.keys())


def test_dir_of_merged() -> None:
    """
    Test loading a directory of merged
    """
    ens: Ensemble = Ensemble.build_from_dir(Path(f'{resources_dir}/docking_merged_test'))

    assert len(ens) == 2
    assert f'{ens}' == "<Ensemble Containing 2 System>"
