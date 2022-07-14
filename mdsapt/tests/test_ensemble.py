from pathlib import Path
from typing import List, Literal

import pytest
from pydantic import FilePath

from ..config import DockingAnalysisConfig
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
    merged_list = [f'{resources_dir}/docking_merged_test/2hnt_15U0.pdb',
                   f'{resources_dir}/docking_merged_test/2hnt_98P.pdb']

    ens: Ensemble = Ensemble.build_from_files(merged_list)
    assert all(k in merged_list for k in ens.keys())


def test_dir_of_merged() -> None:
    ens: Ensemble = Ensemble.build_from_dir(f'{resources_dir}/docking_merged_test')

    assert len(ens) == 2

