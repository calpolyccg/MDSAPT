from pathlib import Path
from typing import List, Literal

import pytest
from pydantic import FilePath

from ..config import DockingStructureMode, DockingAnalysisConfig
from ..utils.ensemble import Ensemble

resources_dir = Path(__file__).parent / 'testing_resources'


working_docking_settings = dict(
    type='docking',
    mode=DockingStructureMode.SeparateLigand,
    protein=f'{resources_dir}/docking_sep_test/2hnt.pdb',
    ligands=[f'{resources_dir}/docking_sep_test/ligands/15U0.pdb',
             f'{resources_dir}/docking_sep_test/ligands/15U1.pdb',
             f'{resources_dir}/docking_sep_test/ligands/98P_2.pdb'],
    combined_topologies=[f'{resources_dir}/docking_merged_test/2hnt_15U0.pdb',
                         f'{resources_dir}/docking_merged_test/2hnt_98P.pdb'],
    pairs=[(3, 14)]
)

docking_config: DockingAnalysisConfig = DockingAnalysisConfig(**working_docking_settings)


def test_list_of_merged() -> None:
    merged_list = [f'{resources_dir}/docking_merged_test/2hnt_15U0.pdb',
                   f'{resources_dir}/docking_merged_test/2hnt_98P.pdb']

    ens: Ensemble = Ensemble(DockingStructureMode.MergedLigand,
                             systems_dir=docking_config.combined_topologies)



