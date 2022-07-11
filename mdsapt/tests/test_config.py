from pathlib import Path
from typing import Dict, Any, Tuple, List

import pydantic
import pytest

from pydantic import ValidationError, FilePath

from mdsapt.config import load_from_yaml_file, RangeFrameSelection, \
    TrajectoryAnalysisConfig, DockingStructureMode, DockingAnalysisConfig, TopologySelection

resources_dir = Path(__file__).parent / 'testing_resources'

traj_setting_list: List[Tuple[str, Any]] = [
    ('trajectories', [f'{resources_dir}/test_read_error.dcd']),
    ('frames', {'start': 1, 'stop': 120}),
    ('frames', [1, 4, 6, 120]),
    ('pairs', [(250, 251)])
]

working_docking_settings = dict(
    mode=DockingStructureMode.SeparateLigand,
    protein=f'{resources_dir}/docking_sep_test/2hnt.pdb',
    ligands=[f'{resources_dir}/docking_sep_test/15U0.pdb',
             f'{resources_dir}/docking_sep_test/15U1.pdb',
             f'{resources_dir}/docking_sep_test/98P_2.pdb'],
    combined_topologies=[f'{resources_dir}/docking_merged_test/2hnt_15U0.pdb',
                         f'{resources_dir}/docking_merged_test/2hnt_98P.pdb'],
    pairs=[('L', 1)]
)

working_docking_settings_list: List[Dict[str, Any]] = [

    # Testing missing protein
    dict(
        mode=DockingStructureMode.SeparateLigand,
        ligands=[f'{resources_dir}/docking_sep_test/15U0.pdb',
                 f'{resources_dir}/docking_sep_test/15U1.pdb',
                 f'{resources_dir}/docking_sep_test/98P_2.pdb'],
        combined_topologies=[f'{resources_dir}/docking_merged_test/2hnt_15U0.pdb',
                             f'{resources_dir}/docking_merged_test/2hnt_98P.pdb'],
        pairs=[('L', 1)]
    ),

    # Testing missing ligand
    dict(
        mode=DockingStructureMode.SeparateLigand,
        protein=f'{resources_dir}/docking_sep_test/2hnt.pdb',
        combined_topologies=[f'{resources_dir}/docking_merged_test/2hnt_15U0.pdb',
                             f'{resources_dir}/docking_merged_test/2hnt_98P.pdb'],
        pairs=[('L', 1)]
    ),

    # Testing missing combined topologies
    dict(
        mode=DockingStructureMode.SeparateLigand,
        ligands=[f'{resources_dir}/docking_sep_test/15U0.pdb',
                 f'{resources_dir}/docking_sep_test/15U1.pdb',
                 f'{resources_dir}/docking_sep_test/98P_2.pdb'],
        pairs=[('L', 1)]
    ),

    # Testing separate with protein set to None
    dict(
        mode=DockingStructureMode.SeparateLigand,
        protein=None,
        ligands=[f'{resources_dir}/docking_sep_test/15U0.pdb',
                 f'{resources_dir}/docking_sep_test/15U1.pdb',
                 f'{resources_dir}/docking_sep_test/98P_2.pdb'],
        combined_topologies=[f'{resources_dir}/docking_merged_test/2hnt_15U0.pdb',
                             f'{resources_dir}/docking_merged_test/2hnt_98P.pdb'],
        pairs=[('L', 1)]
    ),

    # Testing separate with ligand set to None
    dict(
        mode=DockingStructureMode.SeparateLigand,
        protein=f'{resources_dir}/docking_sep_test/2hnt.pdb',
        ligands=None,
        combined_topologies=[f'{resources_dir}/docking_merged_test/2hnt_15U0.pdb',
                             f'{resources_dir}/docking_merged_test/2hnt_98P.pdb'],
        pairs=[('L', 1)]
    ),

    # Testing separate with protein and ligand set to None
    dict(
        mode=DockingStructureMode.SeparateLigand,
        protein=None,
        ligands=None,
        combined_topologies=[f'{resources_dir}/docking_merged_test/2hnt_15U0.pdb',
                             f'{resources_dir}/docking_merged_test/2hnt_98P.pdb'],
        pairs=[('L', 1)]
    ),

    # Testing combine with top set to none
    dict(
        mode=DockingStructureMode.MergedLigand,
        protein=f'{resources_dir}/docking_sep_test/2hnt.pdb',
        ligands=[f'{resources_dir}/docking_sep_test/15U0.pdb',
                 f'{resources_dir}/docking_sep_test/15U1.pdb',
                 f'{resources_dir}/docking_sep_test/98P_2.pdb'],
        combined_topologies=None,
        pairs=[('L', 1)]
    ),


]


def test_frame_range_selection():
    frame_range: Dict[str, int] = {'start': 2,
                                   'step': 1,
                                   'stop': 0}

    with pytest.raises(ValidationError):
        RangeFrameSelection(**frame_range)


@pytest.mark.parametrize('setting', traj_setting_list)
def test_traj_analysis_config(setting: Tuple[str, Any]) -> None:
    traj_analysis_dict: Dict[str, Any] = dict(
        topology=f'{resources_dir}/testtop.psf',
        trajectories=[f'{resources_dir}/testtraj.dcd'],
        pairs=[(132, 152), (34, 152)],
        frames=[1, 4, 6],
        output=True)

    traj_analysis_dict[setting[0]] = setting[1]

    with pytest.raises(ValidationError):
        TrajectoryAnalysisConfig(**traj_analysis_dict)


@pytest.mark.parametrize('setting', working_docking_settings_list)
def test_dock_analysis_config(setting: Dict[str, Any]) -> None:

    with pytest.raises(ValidationError):
        DockingAnalysisConfig(**setting)


def test_template_can_be_loaded():
    load_from_yaml_file(resources_dir / "test_input.yaml")


def test_load_from_yaml_fail():
    with pytest.raises(ValidationError):
        load_from_yaml_file(resources_dir / 'test_broken_input.yaml')


def test_topology_selection_parses_from_string():
    data: str = str(resources_dir / 'test_input.yaml')
    result = pydantic.parse_obj_as(TopologySelection, data)

    assert result.path == resources_dir / 'test_input.yaml'
    assert result.charge_overrides == {}


def test_topology_selection_parses_from_obj_with_overrides():
    data = {'path': str(resources_dir / 'test_input.yaml'), 'charge_overrides': {'13': 3}}
    result = pydantic.parse_obj_as(TopologySelection, data)

    assert result.path == resources_dir / 'test_input.yaml'
    assert result.charge_overrides == {13: 3}


def test_topology_selection_parses_from_obj_without_overrides():
    data = {'path': str(resources_dir / 'test_input.yaml')}
    result = pydantic.parse_obj_as(TopologySelection, data)

    assert result.path == resources_dir / 'test_input.yaml'
    assert result.charge_overrides == {}
