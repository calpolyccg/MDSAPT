r"""
:mod:`mdsapt.reader` -- Reads input file and saves configuration
================================================================

MDSAPT uses an yaml file to get user's configurations for SAPT calculations
class`mdsapt.reader.InputReader` is responsible for reading the yaml file and
returning the information from it. If a yaml file is needed it can be generated
using the included *mdsapt_get_runinput* script.

"""

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

    @classmethod
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
        except (mda.exceptions.NoDataError, OSError, ValueError):
            errors.append('Error while creating universe using provided topology and trajectories')
            raise ValidationError(errors)  # If Universe doesn't load need to stop

        # Ensure that
        items: Set[int] = {i for pair in ag_pair for i in pair}

        for sel in items:
            ag: mda.AtomGroup = unv.select_atoms(f'resid {sel}')
            if len(ag) == 0:
                errors.append(f"Selection {sel} returns an empty AtomGroup.")



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

    def get_universe(self, **universe_kwargs) -> mda.Universe:
        return mda.Universe(str(topology_selection_path(self.topology)),
                            [str(path) for path in self.trajectories],
                            **universe_kwargs)

    def get_selections(self) -> Set[int]:
        return {i for pair in self.pairs for i in pair}


class DockingStructureMode(Enum):
    Protein_Ligand = 'protein-ligand'
    Separate_Ligand = 'separate-ligand'


DockingElement = Union[Literal['L'], int]


class DockingAnalysisConfig(BaseModel):
    type: Literal['docking']
    mode: DockingStructureMode
    protein: Optional[TopologySelection]
    ligands: Optional[Union[List[TopologySelection], DirectoryPath]]
    combined_topologies: Union[List[TopologySelection], DirectoryPath]
    pairs: List[Tuple[DockingElement, DockingElement]]

    @root_validator()
    def check_valid_config(self, values: Dict[str, Any]) -> Dict[str, Any]:
        mode: DockingStructureMode = values['mode']
        protein: Optional[TopologySelection] = None
        ligands: Optional[Union[List[TopologySelection], DirectoryPath]] = None
        combined_topologies: Optional[Union[List[TopologySelection], DirectoryPath]] = None
        errors: List[str] = []

        # Trying to get values some will not be given
        try:
            protein = values['protein']
        except KeyError:
            pass

        try:
            ligands = values['ligands']
        except KeyError:
            pass
        try:
            combined_topologies = values['combined_topologies']
        except KeyError:
            pass

        # Checking if necessary values given
        if mode == DockingStructureMode.Protein_Ligand and \
                combined_topologies is None:
            errors.append("protein and ligands must be specified when using 'protein-ligand' mode")
        elif mode == DockingStructureMode.Separate_Ligand and \
                (protein is None or ligands is None):
            errors.append("topologies must be specified with using 'combined-topologies mode")

        self.replace_ligand_alias()

        pairs: List[Tuple[DockingElement, DockingElement]] = values['pairs']
        selections: Set[DockingElement] = {
            i for pair in self.pairs for i in pair
        }

        if mode == DockingStructureMode.Protein_Ligand:
            sys_dict: Dict[TopologySelection, mda.Universe] = {}

            for top in combined_topologies:
                try:
                    sys_dict[top] = mda.Universe(str(topology_selection_path(top)))
                except (mda.exceptions.NoDataError, OSError, ValueError):
                    errors.append('Error while creating universe using provided topology and trajectories')
                    raise ValidationError(errors)  # If Universe doesn't load need to stop

            for selection in selections:
                for k in sys_dict:
                    ag: mda.AtomGroup = sys_dict[k].select_atoms(f'resid {selection}')
                    if len(ag) == 0:
                        errors.append(f"Selection {selection} returns an empty AtomGroup.")
        elif mode == DockingStructureMode.Separate_Ligand:
            try:
                protein_sys: mda.Universe = mda.Universe(str(protein))
            except (mda.exceptions.NoDataError, OSError, ValueError):
                errors.append('Error while creating universe using provided topology and trajectories')
                raise ValidationError(errors)  # If Universe doesn't load need to stop

            for selection in selections:
                ag: mda.AtomGroup = protein_sys.select_atoms(f'resid {selection}')
                if len(ag) == 0:
                    errors.append(f"Selection {selection} returns an empty AtomGroup.")

        return values

    def replace_ligand_alias(self) -> None:
        for ind, pair in enumerate(self.pairs):
            if pair[0] == Literal['L']:
                self.pairs[ind] = (-1, self.pairs[ind][1])
            if pair[1] == Literal['L']:
                self.pairs[ind] = (self.pairs[ind][0], -1)

    def get_selections(self) -> Set[Union[Literal['L'], int]]:
        return {i for pair in self.pairs for i in pair}


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
