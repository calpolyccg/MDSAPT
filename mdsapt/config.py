r"""
:mod:`mdsapt.reader` -- Reads input file and saves configuration
================================================================

MDSAPT uses an yaml file to get user's configurations for SAPT calculations
class`mdsapt.reader.InputReader` is responsible for reading the yaml file and
returning the information from it. If a yaml file is needed it can be generated
using the included *mdsapt_get_runinput* script.

.. autoexception:: ConfigurationError

.. autoclass:: InputReader
    :members:
    :inherited-members:
"""

import os
from enum import Enum
from pathlib import Path
from typing import List, Dict, Tuple, Literal, Optional, \
    Union, Any, Set

import yaml

import MDAnalysis as mda

import logging

from pydantic import BaseModel, conint, Field, root_validator, \
    FilePath, ValidationError, DirectoryPath

logger = logging.getLogger('mdsapt.config')


class Psi4Config(BaseModel):
    """Psi4 configuration details

    The SAPT method to use.

    NOTE: You can use any valid Psi4 method, but it might fail if you don't use a SAPT method.

    The basis to use in Psi4.

    NOTE: We do not verify if this is a valid basis set or not.
    """

    method: str
    basis: str
    save_output: bool  # whether to save the raw output of Psi4. May be useful for debugging.
    settings: Dict[str, str]  # Other Psi4 settings you would like to provide.


class SysLimitsConfig(BaseModel):
    """Resource limits for your system."""
    ncpus: conint(ge=1)
    memory: str


class ChargeGuesser(Enum):
    Standard = 'standard'
    RDKit = 'rdkit'


class SimulationConfig(BaseModel):
    ph: float
    charge_guesser: ChargeGuesser


class DetailedTopologySelection(BaseModel):
    path: FilePath
    charge_overrides: Dict[int, int] = Field(default_factory=dict)


TopologySelection = Union[FilePath, DetailedTopologySelection]


def topology_selection_path(sel: TopologySelection) -> FilePath:
    if isinstance(sel, DetailedTopologySelection):
        return sel.path
    return sel


class RangeFrameSelection(BaseModel):
    start: Optional[conint(ge=0)]
    stop: Optional[conint(ge=0)]
    step: Optional[conint(ge=1)]

    @root_validator()
    def check_start_before_stop(cls, values: Dict[str, int]) -> Dict[str, int]:
        assert values['start'] <= values['stop'], "start must be before stop"
        return values


class TrajectoryAnalysisConfig(BaseModel):
    """
    A selection of the frames used in this analysis.

    Serialization behavior
    ----------------------
    If this value is a range, it will be serialized using start/stop/step.
    Otherwise, it will be serialized into a List[int].
    """
    type: Literal['trajectory']
    topology: TopologySelection
    trajectories: List[FilePath]
    pairs: List[Tuple[conint(ge=0), conint(ge=0)]]
    frames: Union[List[int], RangeFrameSelection]
    output: str

    @root_validator()
    def check_valid_md_system(cls, values: Dict[str, Any]) -> Dict[str, Any]:
        errors: List[str] = []

        top_path: TopologySelection = values['topology']
        trj_path: List[FilePath] = values['trajectories']
        ag_pair: List[Tuple[conint(ge=0), conint(ge=0)]] = values['pairs']
        frames: Union[List[int], RangeFrameSelection] = values['frames']

        try:
            unv = mda.Universe(str(topology_selection_path(top_path)),
                               [str(p) for p in trj_path])
        except mda.exceptions.NoDataError or ValueError:
            raise ValueError('Error while creating universe using provided topology and trajectories')

        # Ensure that
        items: Set[int] = {i for pair in ag_pair for i in pair}

        for sel in items:
            try:
                unv.select_atoms(f'resid {sel} and protein')
            except mda.SelectionError:
                errors.append('Error in selection: {}'.format(sel))

        trajlen: int = len(unv.trajectory)
        if isinstance(frames, RangeFrameSelection):
            if trajlen <= frames.stop:
                errors.append(f'Stop {frames.stop} exceeds trajectory length {trajlen}.')
        else:
            frames: List[int]
            for frame in frames:
                if frame >= trajlen:
                    errors.append(f'Frame {frame} exceeds trajectory length {trajlen}')

        if len(errors) > 0:
            raise ValidationError(errors)
        return values


class DockingAnalysisConfig(BaseModel):
    type: Literal['docking']
    topologies: Union[List[TopologySelection], DirectoryPath]
    pairs: List[Tuple[int, int]]


class Config(BaseModel):
    psi4: Psi4Config
    simulation: SimulationConfig
    system_limits: SysLimitsConfig
    analysis: Union[TrajectoryAnalysisConfig, DockingAnalysisConfig] = Field(..., discriminator='type')


def load_from_yaml_file(path: Union[str, Path]) -> Config:
    if isinstance(path, str):
        path = Path(path)

    with path.open() as f:
        try:
            return Config(**yaml.safe_load(f))
        except ValidationError as err:
            logger.exception(f"Error while loading {path}")
            raise err
