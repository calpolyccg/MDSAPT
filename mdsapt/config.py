r"""
:mod:`mdsapt.config` -- Reads input file and saves configuration
================================================================

"""
import dataclasses
from dataclasses import dataclass
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import List, Dict, Tuple, Literal, Optional, \
    Union, Any, Set, Iterable

import logging

import pydantic
from pydantic import BaseModel, conint, Field, root_validator, \
    FilePath, ValidationError, DirectoryPath
import yaml

import MDAnalysis as mda

from mdsapt.utils.ensemble import Ensemble

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
    STANDARD = 'standard'
    RDKIT = 'rdkit'


class SimulationConfig(BaseModel):
    ph: float
    charge_guesser: ChargeGuesser


@dataclass
class TopologySelection:
    class _TopologySelection(BaseModel):
        path: FilePath
        topology_format: Optional[str]
        charge_overrides: Dict[int, int] = Field(default_factory=dict)

    path: Path
    topology_format: Optional[str] = None
    charge_overrides: Dict[int, int] = dataclasses.field(default_factory=dict)

    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, values):
        result = pydantic.parse_obj_as(Union[FilePath, cls._TopologySelection], values)
        if isinstance(result, PathLike):
            return TopologySelection(path=Path(result))
        return TopologySelection(path=result.path, topology_format=result.topology_format,
                                 charge_overrides=result.charge_overrides)

    def create_universe(self, *coordinates: Any, **kwargs) -> mda.Universe:
        """Create a universe based on this topology and the given arguments.."""
        return mda.Universe(str(self.path), *coordinates,
                            topology_format=self.topology_format, **kwargs)


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

    # noinspection PyMethodParameters
    @root_validator
    def check_valid_md_system(cls, values: Dict[str, Any]) -> Dict[str, Any]:
        errors: List[str] = []

        topology: TopologySelection = values['topology']
        trajectories: List[FilePath] = values['trajectories']
        ag_pair: List[Tuple[conint(ge=0), conint(ge=0)]] = values['pairs']
        frames: Union[List[int], RangeFrameSelection] = values['frames']

        try:
            unv = topology.create_universe([str(p) for p in trajectories])
        except OSError as e:
            raise ValueError("Error while creating the universe") from e

        missing_selections = get_invalid_residue_selections({r for p in ag_pair for r in p}, unv)
        if len(missing_selections) > 0:
            errors.append(f'Selected residues are missing from topology: {missing_selections}')

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
            raise ValidationError([errors], cls)
        return values

    def create_universe(self, **universe_kwargs) -> mda.Universe:
        return self.topology.create_universe([str(p) for p in self.trajectories], **universe_kwargs)

    def get_selections(self) -> Set[int]:
        return {i for pair in self.pairs for i in pair}


def get_invalid_residue_selections(residues: Iterable[int], unv: mda.Universe) -> Iterable[int]:
    """Helper function to find selected residues that aren't in the universe."""
    return [
        i for i in residues
        if len(unv.select_atoms(f'resid {i}')) == 0
    ]


DockingElement = Union[Literal['L'], conint(ge=-1)]
"""
A single element to analyze in docking.

The literal 'L' specifies the ligand, whereas an integer specifies the protein residue number.
"""


class TopologyGroupSelection(BaseModel):
    __root__: Union[DirectoryPath, List[TopologySelection]]

    def get_individual_topologies(self) -> List[TopologySelection]:
        if isinstance(self.__root__, list):
            return self.__root__

        return [
            TopologySelection(path=f)
            for f in self.__root__.iterdir()
            if f.is_file()
        ]


# noinspection PyMethodParameters
class DockingAnalysisConfig(BaseModel):
    type: Literal['docking']
    pairs: List[Tuple[DockingElement, DockingElement]]
    combined_topologies: Optional[TopologyGroupSelection]
    protein: Optional[TopologySelection]
    ligands: Optional[TopologyGroupSelection]
    output: str

    @root_validator
    def check_valid_config(cls, values: Dict[str, Any]) -> Dict[str, Any]:
        errors: List[str] = []

        pairs: List[Tuple[DockingElement, DockingElement]] = values['pairs']
        protein_selections: Set[int] = {
            i for pair in pairs for i in pair if i != 'L'
        }
        ens: Ensemble = cls._build_ensemble(combined_topologies=values.get('combined_topologies'),
                                            protein=values.get('protein'),
                                            ligands=values.get('ligands'))
        missing_selections: List[int] = []

        for k in ens.keys():
            missing_selections += get_invalid_residue_selections(protein_selections, ens[k])

        if len(missing_selections) > 0:
            errors.append(f'Selected residues are missing from topology: {missing_selections}')

        return values

    def build_ensemble(self):
        return self._build_ensemble(combined_topologies=self.combined_topologies, protein=self.protein,
                                    ligands=self.ligands)

    @classmethod
    def _build_ensemble(
            cls,
            *,
            combined_topologies: Optional[TopologyGroupSelection],
            protein: Optional[TopologySelection],
            ligands: Optional[TopologyGroupSelection],
    ) -> Ensemble:
        """Fails if the wrong types of arguments are provided."""
        if combined_topologies is not None and (protein, ligands) == (None, None):
            return Ensemble.build_from_files(
                [top.path for top in combined_topologies.get_individual_topologies()]
            )

        if combined_topologies is None and None not in (protein, ligands):
            ens: Ensemble = Ensemble.build_from_files([top.path for top
                                                       in ligands.get_individual_topologies()])
            protein_sys: mda.Universe = mda.Universe(str(protein.path))
            protein_mol: mda.AtomGroup = protein_sys.select_atoms("protein")
            ens = ens.merge(protein_mol)
            return ens

        raise ValueError('Must provide `protein` and `ligands` keys, or only `combined_topologies`')

    def get_selections(self) -> Set[int]:
        return {i for pair in self.pairs for i in pair}


class Config(BaseModel):
    """

    """
    psi4: Psi4Config
    simulation: SimulationConfig
    system_limits: SysLimitsConfig
    analysis: Union[TrajectoryAnalysisConfig, DockingAnalysisConfig] = Field(..., discriminator='type')


def load_from_yaml_file(path: Union[str, Path]) -> Config:
    """

    """
    path = Path(path)

    with path.open() as file:
        try:
            return Config(**yaml.safe_load(file))
        except ValidationError as err:
            logger.exception(f"Error while loading {path}")
            raise err
