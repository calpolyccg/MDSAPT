from pathlib import Path
from typing import Dict, Any, Tuple, List

import pytest

from pydantic import ValidationError, FilePath

from mdsapt.config import load_from_yaml_file, RangeFrameSelection, \
    DetailedTopologySelection, topology_selection_path, \
    TrajectoryAnalysisConfig

resources_dir = Path(__file__).parent / 'testing_resources'

setting_list: List[Tuple[str, Any]] = [
    ('trajectories', [f'{resources_dir}/test_read_error.dcd']),
    ('frames', {'start': 1, 'stop': 120}),
    ('frames', [1, 4, 6, 120]),
    ('pairs', [(250, 251)])
]


def test_frame_range_selection():
    frame_range: Dict[str, int] = {'start': 2,
                                   'step': 1,
                                   'stop': 0}

    with pytest.raises(ValidationError):
        RangeFrameSelection(**frame_range)


@pytest.mark.parametrize('setting', setting_list)
def test_traj_analysis_config(setting: Tuple[str, Any]):
    traj_analysis_dict: Dict[str, Any] = {
        'topology': f'{resources_dir}/testtop.psf',
        'trajectories': [f'{resources_dir}/testtraj.dcd'],
        'pairs': [(132, 152), (34, 152)],
        'frames': [1, 4, 6],
        'output': True
    }

    traj_analysis_dict[setting[0]] = setting[1]

    with pytest.raises(ValidationError):
        TrajectoryAnalysisConfig(**traj_analysis_dict)


def test_template_can_be_loaded():
    load_from_yaml_file(resources_dir / "test_input.yaml")

def test_load_from_yaml_fail():
    with pytest.raises(ValidationError):
        load_from_yaml_file(resources_dir / 'test_broken_input.yaml')
