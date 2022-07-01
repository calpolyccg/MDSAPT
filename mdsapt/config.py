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
from typing import List, Dict, Tuple, Iterable, Generic, TypeVar, Self, Literal, Optional, Union
from dataclasses import dataclass
from abc import ABC, abstractmethod

import yaml
import psi4.driver.procrouting.proc_table as proc_table

import MDAnalysis as mda

import logging

from pydantic import BaseModel, validator, conint, Field, root_validator, FilePath, ValidationError

logger = logging.getLogger('mdsapt.reader')


class Psi4Config(BaseModel):
    """Psi4 configuration details"""

    method: str
    """
    The SAPT method to use.
    
    NOTE: You can use any valid Psi4 method, but it might fail if you don't use a SAPT method.
    """

    basis: str
    """
    The basis to use in Psi4.
    
    NOTE: We do not verify if this is a valid basis set or not.
    """
    save_output: bool  # whether to save the raw output of Psi4. May be useful for debugging.

    settings: Dict[str, str]  # Other Psi4 settings you would like to provide.

    @validator('method')
    def check_valid_method(self, v):
        if v not in proc_table.procedures['energy'].keys():
            raise ValueError(f"method {v} not supported by Psi4! See Psi4 docs for list of valid methods.")
        return v


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


class TopologySelection(BaseModel):
    path: str
    charge_overrides: Dict[int, int]

    @validator('path')
    def check_valid_path(self, v):
        path = os.path.join(os.curdir, v)
        if not os.path.exists(path):
            raise FileNotFoundError(f"specified topology file: {path} not found")
        return v


class RangeFrameSelection(BaseModel):
    start: Optional[conint(ge=0)] = None
    stop: Optional[conint(ge=0)] = None
    step: Optional[conint(ge=1)] = None

    @root_validator()
    def check_start_before_stop(self, values):
        assert values['start'] <= values['stop'], "start must be before stop"
        return values


class TrajectoryAnalysisConfig(BaseModel):
    type: Literal['trajectory']
    topology: TopologySelection
    trajectories: List[FilePath]
    pairs: List[Tuple[int, int]]
    frames: Union[List[int], RangeFrameSelection]
    """
    A selection of the frames used in this analysis.
    
    Serialization behavior
    ----------------------
    If this value is a range, it will be serialized using start/stop/step.
    Otherwise, it will be serialized into a List[int].
    """
    output: str

    @root_validator()
    def check_valid_md_system(self, values):
        errors = []

        top_path = values['topology']
        trj_path = values['trajectories']
        ag_pair = values['pairs']
        frames = values['frames']

        try:
            unv = mda.Universe(top_path, trj_path)
        except mda.exceptions.NoDataError or ValueError:
            raise ValueError('Error while creating universe using provided topology and trajectories')

        # Ensure that
        items = {i for pair in ag_pair for i in pair}
        for sel in items:
            try:
                unv.select_atoms(f'resid {sel} and protein')
            except mda.SelectionError:
                errors.append('Error in selection: {}'.format(sel))

        trajlen = len(unv.trajectory)
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
    topologies: List[str]
    pairs: List[Tuple[int, int]]
    output: str


class Config(BaseModel):
    psi4: Psi4Config
    simulation: SimulationConfig
    system_limits: SysLimitsConfig
    analysis: Union[TrajectoryAnalysisConfig, DockingAnalysisConfig] = Field(..., discriminator='type')


class InputReader:
    """Reader for yaml inputs"""

    top_path: str
    trj_path: str
    ag_sel: List[int]
    ag_pair: List[List[int]]
    trj_settings: dict
    sys_settings: dict
    opt_settings: dict
    sapt_settings: dict
    sapt_method: str
    sapt_basis: str
    sapt_out: bool
    start: int
    stop: int
    step: int
    pH: float
    ncpus: int
    memory: int
    walltime: str

    def __init__(self, path) -> None:
        """Reads input file, checks it for validity,
         and saves its data as instance variables.

         Errors in run input will result in
         :class`mdsapt.reader.ConfigurationError` being
         raised."""
        self.load_input(path)

    def load_input(self, path: str) -> None:
        """Loads input file from path and records settings.
         If an error is found :class:`mdsapt.reader.ConfigurationError`
         is raised.

         :Arguments:
            *path*
                Path to yaml input file."""
        try:
            in_cfg = yaml.safe_load(open(path))
            self._check_inputs(in_cfg)
            self._save_params(in_cfg)
        except IOError or ConfigurationError:
            logger.fatal(f'error loading file {path}')
            raise ConfigurationError

    def _save_params(self, yaml_dict: dict) -> None:
        self.top_path = yaml_dict['topology_path']
        self.trj_path = yaml_dict['trajectory_paths']
        self.ag_sel = yaml_dict['selection_resid_num']
        self.ag_pair = yaml_dict['int_pairs']
        self.trj_settings = yaml_dict['trajectory_settings']
        self.sys_settings = yaml_dict['system_settings']
        self.opt_settings = yaml_dict['opt_settings']
        self.sapt_settings = yaml_dict['sapt_settings']
        self.sapt_method = yaml_dict['sapt_settings']['method']
        self.sapt_basis = yaml_dict['sapt_settings']['basis']
        self.sapt_out = yaml_dict['sapt_settings']['save_psi4_output']
        self.start = yaml_dict['trajectory_settings']['start']
        self.stop = yaml_dict['trajectory_settings']['stop']
        self.step = yaml_dict['trajectory_settings']['step']
        self.pH = yaml_dict['opt_settings']['pH']
        self.ncpus = yaml_dict['system_settings']['ncpus']
        self.memory = yaml_dict['system_settings']['memory']
        self.walltime = yaml_dict['system_settings']['time']

    @staticmethod
    def _check_inputs(yaml_dict: dict) -> None:
        # Checking inputs of yaml file
        try:
            top_path = yaml_dict['topology_path']
            trj_path = yaml_dict['trajectory_paths']
            ag_sel = yaml_dict['selection_resid_num']
            ag_pair = yaml_dict['int_pairs']
            trj_settings = yaml_dict['trajectory_settings']
            sys_settings = yaml_dict['system_settings']
            opt_settings = yaml_dict['opt_settings']
            sapt_settings = yaml_dict['sapt_settings']
        except KeyError as err:
            logger.fatal(f'{err}: missing from YAML file')
            raise ConfigurationError

        try:
            if not os.path.exists(os.path.join(os.getcwd(), top_path)):
                raise ConfigurationError
            for f in trj_path:
                if not os.path.exists(os.path.join(os.getcwd(), f)):
                    raise ConfigurationError
            unv = mda.Universe(os.path.join(os.getcwd(), top_path), [os.path.join(os.getcwd(), x) for x in trj_path])
        except mda.exceptions.NoDataError or ConfigurationError or ValueError:
            logger.fatal('MD file error')
            raise ConfigurationError

        # Testing names and selections
        for sel in ag_sel:
            try:
                ag = unv.select_atoms(f'resid {sel} and protein')
            except mda.SelectionError:
                raise ConfigurationError('Error in selection: {}'.format(sel))

        try:
            start = trj_settings['start']
            step = trj_settings['step']
            stop = trj_settings['stop']
            cpu = sys_settings['ncpus']
            mem = sys_settings['memory']
            time = sys_settings['time']
            pH = opt_settings['pH']
            method = sapt_settings['method']
            basis = sapt_settings['basis']
            settings = sapt_settings['settings']
            save_sapt_out = sapt_settings['save_psi4_output']

        except KeyError as err:

            logger.fatal(f'{err}: missing data in input file')
            raise ConfigurationError

        for pair in ag_pair:
            if len(pair) != 2:
                logger.fatal('Pairs must be a python list of integers with 2 items')
                raise ConfigurationError
            found0 = False
            found1 = False
            for name in ag_sel:
                if pair[0] == name:
                    found0 = True
                if pair[1] == name:
                    found1 = True
            if found0 is False:
                logger.fatal(f'{pair[0]} in {pair} group_pair_selections is not in defined in atom_group_names')
                raise ConfigurationError
            if found1 is False:
                logger.fatal(f'{pair[1]} in {pair} group_pair_selections is not in defined in atom_group_names')
                raise ConfigurationError

            if start >= stop:
                logger.fatal('Start is greater than or equal to stop')
                raise ConfigurationError
            if step >= stop:
                logger.fatal('Step is greater than or equal to stop')
                raise ConfigurationError
            if step == 0:
                logger.fatal('Step cannot be 0')
                raise ConfigurationError

            if len(unv.trajectory) < stop:
                logger.fatal('Stop exceeds length of trajectory.')
                raise ConfigurationError

        logger.info('Input Parameters Accepted')
