r"""
:mod:`mdsapt.config` -- Reads input file and saves configuration
================================================================

"""

# There's lots of implicit class methods because pydantic decorators are stupid.
# Thus, we will disable this lint for this file.
# pylint: disable=no-self-argument

import dataclasses
from dataclasses import dataclass
from enum import Enum
from os import PathLike
import os
from pathlib import Path
from typing import List, Dict, Tuple, Literal, Optional, \
    Union, Any, Set, Iterable

import logging

import pydantic
from pydantic import BaseModel, conint, Field, root_validator, \
    FilePath, ValidationError, DirectoryPath
import yaml

import MDAnalysis as mda
from mdsapt.repair import ChargeStrategy, StandardChargeStrategy

from mdsapt.utils.ensemble import Ensemble

logger = logging.getLogger(__name__)


class Psi4Config(BaseModel):
    """Psi4 configuration details

    Attributes:
        method:
            The SAPT method to use.
            NOTE: You can use any valid Psi4 method, but it might fail if you don't
            use a valid one.
        basis:
            The basis to use.
            NOTE: We do not verify if this is a valid basis set or not.
        save_output:
            Whether to save the raw output of Psi4. May be useful for debugging.
        settings:
            Other Psi4 settings you would like to provide. These will be passed into
            `psi4.set_options <https://psicode.org/psi4manual/master/api/psi4.driver.set_options.html>`_.
    """

    method: str
    basis: str
    save_output: bool
    settings: Dict[str, str]


class SysLimitsConfig(BaseModel):
    """
    Resource limits for your system.
    """
    ncpus: conint(ge=1)
    memory: str


class ChargeGuesser(Enum):
    """
    Specifies the charge guesser used in the analysis. The standard one is faster but will not
    work on more complex sulfur-based ligands.
    """
    STANDARD = 'standard'
    RDKIT = 'rdkit'

    @property
    def charge_strategy(self) -> ChargeStrategy:
        if self == ChargeGuesser.STANDARD:
            return StandardChargeStrategy()
        if self == ChargeGuesser.RDKIT:
            raise NotImplemented('RDKit guesser is not implemented yet.')
        raise ValueError(f'Unknown charge guesser {self}')


class SimulationConfig(BaseModel):
    """
    Configuration options for the simulation.

    Attributes:
        ph (float): The pH of the simulation.
        charge_guesser: The default charge guesser to use for atoms.
    """
    ph: float
    charge_guesser: ChargeGuesser


@dataclass
class TopologySelection:
    """
    A configuration item for selecting a single topology. To successfully
    import a topology, it must be supported by MDAnalysis.

    Attributes:
        path: Where the topology file is located.
        topology_format: If specified, overrides the format to import with.
            charge_overrides: An optional dictionary, where keys are atom numbers and
            values are their charges.

    .. seealso::
        `List of topology formats that MDAnalysis supports <https://docs.mdanalysis.org/1.1.1/documentation_pages/topology/init.html>`_
    """
    class _TopologySelection(BaseModel):
        path: FilePath
        topology_format: Optional[str]
        charge_overrides: Dict[int, int] = Field(default_factory=dict)

    path: Path
    topology_format: Optional[str] = None
    charge_overrides: Dict[int, int] = dataclasses.field(default_factory=dict)

    @classmethod
    def __get_validators__(cls):
        yield cls._validate

    @classmethod
    def _validate(cls, values):
        """
        Validates the topology. You should not call this directly.
        """
        result = pydantic.parse_obj_as(Union[FilePath, cls._TopologySelection], values)
        try:
            path = Path(result)
        except TypeError:
            return TopologySelection(path=result.path, topology_format=result.topology_format,
                                     charge_overrides=result.charge_overrides)
        return TopologySelection(path=path)

    def create_universe(self, *coordinates: Any, **kwargs) -> mda.Universe:
        """Create a universe based on this topology and the given arguments.."""
        return mda.Universe(str(self.path), *coordinates,
                            topology_format=self.topology_format, **kwargs)


class RangeFrameSelection(BaseModel):
    """
    A selection of contiguous frames.

    Attributes:
        start: the first frame to use, inclusive.
        stop: the last frame to use, inclusive.
        step: step between frames.
    """
    start: Optional[conint(ge=0)]
    stop: Optional[conint(ge=0)]
    step: Optional[conint(ge=1)] = 1

    @root_validator()
    def _check_start_before_stop(cls, values: Dict[str, int]) -> Dict[str, int]:
        """
        Ensures that a valid range is selected for frame iteration.
        """
        assert values['start'] <= values['stop'], "start must be before stop"
        return values


class TrajectoryAnalysisConfig(BaseModel):
    """
    Config for performing a Trajectory SAPT analysis.

    .. seealso::
        :obj:`mdsapt.sapt.TrajectorySAPT`

    Attributes:
        topology: The topology to analyze.
        trajectories: A list of trajectories to analyze.
        pairs: Interaction pairs to study
        frames:
            A selection of frames to analyze.

            This must be a :obj:`RangeFrameSelection` with start/stop/step.
        output: A file to write an output CSV to.
    """
    type: Literal['trajectory']
    topology: TopologySelection
    trajectories: List[FilePath]
    pairs: List[Tuple[conint(ge=0), conint(ge=0)]]
    frames: RangeFrameSelection
    output: str

    # noinspection PyMethodParameters
    @root_validator
    def check_valid_md_system(cls, values: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validates that setting work with the selected MD system
        """
        errors: List[str] = []

        topology: TopologySelection = values['topology']
        trajectories: List[FilePath] = values['trajectories']
        ag_pair: List[Tuple[conint(ge=0), conint(ge=0)]] = values['pairs']
        frames: RangeFrameSelection = values['frames']

        try:
            unv = topology.create_universe([str(p) for p in trajectories])
        except OSError as err:
            raise ValueError("Error while creating the universe") from err

        missing_selections = get_invalid_residue_selections({r for p in ag_pair for r in p}, unv)
        if len(missing_selections) > 0:
            errors.append(f'Selected residues are missing from topology: {missing_selections}')

        trajlen: int = len(unv.trajectory)
        if trajlen <= frames.stop:
            errors.append(f'Stop {frames.stop} exceeds trajectory length {trajlen}.')

        if len(errors) > 0:
            raise ValidationError([errors], cls)
        return values

    def create_universe(self, **universe_kwargs) -> mda.Universe:
        """
        Loads a universe from the given topology and trajectory
        """
        return self.topology.create_universe([str(p) for p in self.trajectories], **universe_kwargs)

    def get_selections(self) -> Set[int]:
        """
        Returns a set of ints for the selected residues for trajectory analysis
        """
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
    """
    A selection of a group of topologies.

    In a YAML config, this may either be a path to a flat directory full of topologies
    or a list of :obj:`TopologySelection`s.
    """
    __root__: Union[DirectoryPath, List[TopologySelection]]

    def get_individual_topologies(self) -> List[TopologySelection]:
        """
        It
        """
        if isinstance(self.__root__, list):
            return self.__root__

        return [
            TopologySelection(path=f)
            for f in self.__root__.iterdir()
            if f.is_file()
        ]


# noinspection PyMethodParameters
class DockingAnalysisConfig(BaseModel):
    """
    Config for performing a docking SAPT analysis.

    There are two valid modes of selecting the system to analyze:
        - Only specifying <combined_topologies>
        - Specifying <protein> and <ligands> together

    You must choose one mode or the other, you cannot mix the two (i.e. specify
    <combined_topologies> and <protein>).

    .. seealso::
        :obj:`mdsapt.sapt.DockingSAPT`

    Attributes:
        pairs: Interaction pairs to study.
        combined_topologies:
            A selection of topologies such that every topology contains a ligand
            and a protein.

            NOTE: These must all be numbered consistently.
        protein:
            The topology of the protein to study.
        ligands:
            A selection of topologies of a ligand such that every ligand is offset relative to
            the protein specified in <protein>.
        output: A file to write an output CSV to.
    """
    type: Literal['docking']
    pairs: List[Tuple[DockingElement, DockingElement]]
    combined_topologies: Optional[TopologyGroupSelection]
    protein: Optional[TopologySelection]
    ligands: Optional[TopologyGroupSelection]
    output: str

    @root_validator
    def _check_valid_config(cls, values: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validates that the provided settings are valid.
        """
        errors: List[str] = []

        pairs: List[Tuple[DockingElement, DockingElement]] = values['pairs']
        protein_selections: Set[int] = {
            i for pair in pairs for i in pair if i != 'L'
        }
        ens: Ensemble = cls._build_ensemble(combined_topologies=values.get('combined_topologies'),
                                            protein=values.get('protein'),
                                            ligands=values.get('ligands'))
        missing_selections: List[int] = []

        for val in ens.values():
            missing_selections += get_invalid_residue_selections(protein_selections, val)

        if len(missing_selections) > 0:
            errors.append(f'Selected residues are missing from topology: {missing_selections}')

        return values

    def build_ensemble(self) -> Ensemble:
        """
        Builds an ensemble from this configuration data.
        """
        return self._build_ensemble(combined_topologies=self.combined_topologies,
                                    protein=self.protein,
                                    ligands=self.ligands)

    def build_overrides_dict(self) -> Dict[str, Dict[int, int]]:
        if self.combined_topologies is not None and (self.protein, self.ligands) == (None, None):
            sels = self.combined_topologies.get_individual_topologies()
            return {
                str(top.path): top.charge_overrides
                for top in sels
            }

        if self.combined_topologies is None and None not in (self.protein, self.ligands):
            ligand_sels = self.ligands.get_individual_topologies()
            result = {
                str(top.path): top.charge_overrides
                for top in ligand_sels
            }
            result[str(self.protein.path)] = self.protein.charge_overrides
            return result

        raise ValueError('Must provide `protein` and `ligands` keys, or only `combined_topologies`')

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
        """
        Returns selected residues in a set
        """
        return {i for pair in self.pairs for i in pair}


class Config(BaseModel):
    """
    The root configuration object for MDSAPT.
    """
    psi4: Psi4Config
    simulation: SimulationConfig
    system_limits: SysLimitsConfig
    analysis: Union[TrajectoryAnalysisConfig, DockingAnalysisConfig] = \
        Field(..., discriminator='type')


def load_from_yaml_file(path: Union[str, PathLike]) -> Config:
    """
    Loads a config from a YAML file.
    """
    with Path(path).open('r', encoding='utf8') as file:
        try:
            return Config(**yaml.safe_load(file))
        except ValidationError as err:
            logger.exception("Error while loading config from %r", path)
            raise err
