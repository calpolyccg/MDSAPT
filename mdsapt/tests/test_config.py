"""
Tests for the config objects
"""
from pathlib import Path
from tempfile import mktemp
from typing import Dict, Any

import pydantic
import pytest

from pydantic import ValidationError

from mdsapt.config import load_from_yaml_file, RangeFrameSelection, \
    TrajectoryAnalysisConfig, DockingAnalysisConfig, TopologySelection
from mdsapt.utils.ensemble import Ensemble, EnsembleAtomGroup

resources_dir = Path(__file__).parent / 'testing_resources'


def test_frame_range_selection() -> None:
    """
    Test frame selection object validator.
    """
    frame_range: Dict[str, int] = {'start': 2,
                                   'step': 1,
                                   'stop': 0}

    with pytest.raises(ValidationError):
        RangeFrameSelection(**frame_range)


def test_traj_sel() -> None:
    """
    Test getting set of selections
    """
    traj_analysis_dict: Dict[str, Any] = dict(
        type='trajectory',
        topology=f'{resources_dir}/testtop.psf',
        trajectories=[f'{resources_dir}/testtraj.dcd'],
        pairs=[(132, 152), (34, 152)],
        frames={'start': 1, 'stop': 4},
        output='something.csv')

    cfg: TrajectoryAnalysisConfig = TrajectoryAnalysisConfig.model_validate(traj_analysis_dict)
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
    """
    Test docking analysis config validation
    """
    with pytest.raises(ValidationError):
        DockingAnalysisConfig(**setting)


def test_template_can_be_loaded() -> None:
    """
    Tests that a trajectory analysis can be loaded from a yaml file
    """
    load_from_yaml_file(resources_dir / "test_input.yaml")


def test_load_from_yaml_fail() -> None:
    """
    Tests that errors are raised with an incorrect yaml is loaded
    """
    with pytest.raises(ValidationError):
        load_from_yaml_file(resources_dir / 'test_broken_input.yaml')


def test_topology_selection_parses_from_string() -> None:
    """
    Test topology selection parses
    """
    data: str = str(resources_dir / 'test_input.yaml')
    result = TopologySelection.model_validate(data)

    assert result.path == resources_dir / 'test_input.yaml'
    assert result.charge_overrides == {}


def test_topology_selection_parses_from_obj_with_overrides() -> None:
    """
    Test alternative topology selection
    """
    data = {'path': str(resources_dir / 'test_input.yaml'), 'charge_overrides': {'13': 3}}
    result = TopologySelection.model_validate(data)

    assert result.path == resources_dir / 'test_input.yaml'
    assert result.charge_overrides == {13: 3}


def test_topology_selection_parses_from_obj_without_overrides() -> None:
    """
    Tests no charge overrides topology selection
    """
    data = {'path': str(resources_dir / 'test_input.yaml')}
    result = TopologySelection.model_validate(data)

    assert result.path == resources_dir / 'test_input.yaml'
    assert result.charge_overrides == {}


def test_separated_docking_list() -> None:
    """
    Test seperated docking list setup for correct ensemble
    """
    config_dict: Dict[str, Any] = dict(
        type='docking',
        protein=resources_dir / 'docking_sep_test/2hnt.pdb',
        ligands=[resources_dir / 'docking_sep_test/ligands/15U0.pdb',
                 resources_dir / 'docking_sep_test/ligands/15U1.pdb',
                 resources_dir / 'docking_sep_test/ligands/98P_1.pdb'],
        pairs=[(13, 18)],
        output=mktemp(),
    )

    cfg: DockingAnalysisConfig = DockingAnalysisConfig(**config_dict)
    ens: Ensemble = cfg.build_ensemble()
    assert len(ens) == 3
    ligands: EnsembleAtomGroup = ens.select_atoms('resid -1')
    assert all((len(ag) != 0 for ag in ligands.values()))


def test_seperated_docking_dir() -> None:
    """
    Test seperated docking directory of ligands for correct ensemble
    """
    config_dict: Dict[str, Any] = dict(
        type='docking',
        protein=resources_dir / 'docking_sep_test/2hnt.pdb',
        ligands=resources_dir / 'docking_sep_test/ligands',
        pairs=[(13, 18)],
        output=mktemp(),
    )

    cfg: DockingAnalysisConfig = DockingAnalysisConfig.model_validate(config_dict)
    ens: Ensemble = cfg.build_ensemble()
    assert len(ens) == 7


def test_combined_docking_list() -> None:
    """
    Test combined topologies list ensemble
    """
    config_dict: Dict[str, Any] = dict(
        type='docking',
        combined_topologies=[resources_dir / 'docking_merged_test/2hnt_15U0.pdb',
                             resources_dir / 'docking_merged_test/2hnt_98P.pdb'],
        pairs=[(13, 18)],
        output=mktemp(),
    )

    cfg: DockingAnalysisConfig = DockingAnalysisConfig(**config_dict)
    ens: Ensemble = cfg.build_ensemble()
    assert len(ens) == 2


def test_combined_docking_dir() -> None:
    """
    Test combined topologies given in a directory
    """
    config_dict: Dict[str, Any] = dict(
        type='docking',
        combined_topologies=resources_dir / 'docking_merged_test',
        pairs=[(13, 18)],
        output=mktemp(),
    )

    cfg: DockingAnalysisConfig = DockingAnalysisConfig(**config_dict)
    ens: Ensemble = cfg.build_ensemble()
    assert len(ens) == 2
