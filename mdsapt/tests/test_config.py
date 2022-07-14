from pathlib import Path
from typing import Dict, Any, Tuple, List

import pydantic
import pytest

from pydantic import ValidationError

from mdsapt.config import load_from_yaml_file, RangeFrameSelection, \
    TrajectoryAnalysisConfig, DockingAnalysisConfig, TopologySelection
from mdsapt.utils.ensemble import Ensemble, EnsembleAtomGroup

resources_dir = Path(__file__).parent / 'testing_resources'


def test_frame_range_selection():
    frame_range: Dict[str, int] = {'start': 2,
                                   'step': 1,
                                   'stop': 0}

    with pytest.raises(ValidationError):
        RangeFrameSelection(**frame_range)


@pytest.mark.parametrize('k,v', [
    ('trajectories', [f'{resources_dir}/test_read_error.dcd']),
    ('frames', {'start': 1, 'stop': 120}),
    ('frames', [1, 4, 6, 120]),
    ('pairs', [(250, 251)])
])
def test_traj_analysis_config(k: str, v: Any) -> None:
    traj_analysis_dict: Dict[str, Any] = dict(
        type='trajectory',
        topology=f'{resources_dir}/testtop.psf',
        trajectories=[f'{resources_dir}/testtraj.dcd'],
        pairs=[(132, 152), (34, 152)],
        frames=[1, 4, 6],
        output=True
    )

    traj_analysis_dict[k] = v

    with pytest.raises(ValidationError):
        TrajectoryAnalysisConfig(**traj_analysis_dict)


def test_traj_sel() -> None:
    traj_analysis_dict: Dict[str, Any] = dict(
        type='trajectory',
        topology=f'{resources_dir}/testtop.psf',
        trajectories=[f'{resources_dir}/testtraj.dcd'],
        pairs=[(132, 152), (34, 152)],
        frames=[1, 4, 6],
        output=True)

    cfg: TrajectoryAnalysisConfig = TrajectoryAnalysisConfig(**traj_analysis_dict)
    assert {34, 132, 152} == cfg.get_selections()


@pytest.mark.parametrize('setting', [
    # Testing missing ligand
    dict(
        type='docking',
        protein=f'{resources_dir}/docking_sep_test/2hnt.pdb',
        pairs=[('L', 1)]
    ),

    # Testing missing combined topologies
    dict(
        type='docking',
        pairs=[('L', 1)]
    ),

    # Testing separate with protein set to None
    dict(
        type='docking',
        ligands=[f'{resources_dir}/docking_sep_test/15U0.pdb',
                 f'{resources_dir}/docking_sep_test/15U1.pdb',
                 f'{resources_dir}/docking_sep_test/98P_2.pdb'],
        pairs=[('L', 1)]
    ),

    # Testing separate with ligand set to None
    dict(
        type='docking',
        protein=f'{resources_dir}/docking_sep_test/2hnt.pdb',
        pairs=[('L', 1)]
    ),

    # Testing combine with top set to none
    dict(
        type='docking',
        pairs=[('L', 1)]
    ),

    dict(
        type='docking',
        protein=f'{resources_dir}/docking_sep_test/2hnt.pdb',
        ligands=[f'{resources_dir}/docking_sep_test/15U0.pdb',
                 f'{resources_dir}/docking_sep_test/15U1.pdb',
                 f'{resources_dir}/docking_sep_test/98P_2.pdb'],
        pairs=[(400, 410)]
    ),
])
def test_dock_analysis_config(setting: Dict[str, Any]) -> None:
    with pytest.raises(ValidationError):
        DockingAnalysisConfig(**setting)


def test_template_can_be_loaded() -> None:
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


def test_separated_docking_list() -> None:
    config_dict: Dict[str, Any] = dict(
        type='docking',
        protein=resources_dir / 'docking_sep_test/2hnt.pdb',
        ligands=[resources_dir / 'docking_sep_test/ligands/15U0.pdb',
                 resources_dir / 'docking_sep_test/ligands/15U1.pdb',
                 resources_dir / 'docking_sep_test/ligands/98P_1.pdb'],
        pairs=[(13, 18)]
    )

    cfg: DockingAnalysisConfig = DockingAnalysisConfig(**config_dict)
    ens: Ensemble = cfg.build_ensemble()
    assert len(ens) == 3
    ligands: EnsembleAtomGroup = ens.select_atoms('resid -1')
    assert all([len(ligands[k]) != 0 for k in ligands.keys()])


def test_seperated_docking_dir() -> None:
    config_dict: Dict[str, Any] = dict(
        type='docking',
        protein=resources_dir / 'docking_sep_test/2hnt.pdb',
        ligands=resources_dir / 'docking_sep_test/ligands',
        pairs=[(13, 18)]
    )

    cfg: DockingAnalysisConfig = DockingAnalysisConfig(**config_dict)
    ens: Ensemble = cfg.build_ensemble()
    assert len(ens) == 7


def test_combined_docking_list() -> None:
    config_dict: Dict[str, Any] = dict(
        type='docking',
        combined_topologies=[resources_dir / 'docking_merged_test/2hnt_15U0.pdb',
                             resources_dir / 'docking_merged_test/2hnt_98P.pdb'],
        pairs=[(13, 18)]
    )

    cfg: DockingAnalysisConfig = DockingAnalysisConfig(**config_dict)
    ens: Ensemble = cfg.build_ensemble()
    assert len(ens) == 2


def test_combined_docking_dir() -> None:
    config_dict: Dict[str, Any] = dict(
        type='docking',
        combined_topologies=resources_dir / 'docking_merged_test',
        pairs=[(13, 18)]
    )

    cfg: DockingAnalysisConfig = DockingAnalysisConfig(**config_dict)
    ens: Ensemble = cfg.build_ensemble()
    assert len(ens) == 2
