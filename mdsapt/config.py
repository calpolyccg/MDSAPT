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
from pydantic import BaseModel, TypeAdapter, conint, Field, model_validator, root_validator, \
    FilePath, ValidationError, DirectoryPath
import yaml

import MDAnalysis as mda

from mdsapt.utils.ensemble import Ensemble

logger = logging.getLogger(__name__)


class Psi4Config(BaseModel):
    """Psi4 configuration details

    Attributes:
        method:
            The SAPT method to use.
            NOTE: You can use any valid Psi4 method, but it might fail if you don't use a SAPT method.
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
    A configuration item for selecting a single topology. To successfully import a topology,
    it must be supported by MDAnalysis.

    Attributes:
        path: Where the topology file is located.
        topology_format: If specified, overrides the format to import with.
        charge_overrides: An optional dictionary, where keys are atom numbers and values are their charges.

    .. seealso::
        `List of topology formats that MDAnalysis supports <https://docs.mdanalysis.org/1.1.1/documentation_pages/topology/init.html>`_
    """

    path: Path
    topology_format: Optional[str] = None
    charge_overrides: Dict[int, int] = dataclasses.field(default_factory=dict)

    @model_validator(mode='before')
    @classmethod
    def _accept_bare_string(cls, data: Any) -> Any:
        try:
            return Path(data)
        except TypeError:
            return data

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

    @model_validator(mode='after')
    def _check_start_before_stop(self) -> 'RangeFrameSelection':
        """
        Ensures that a valid range is selected for frame iteration.
        """
        if self.start is not None and self.stop is not None and self.start > self.stop:
            raise ValueError('Start must be before stop')
        return self


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


TopologyGroupSelection = Union[DirectoryPath, List[TopologySelection]]
"""
A selection of a group of topologies.

In a YAML config, this may either be a path to a flat directory full of topologies
or a list of :obj:`TopologySelection`s.
"""


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

    def build_ensemble(self) -> Ensemble:
        return self._build_ensemble(combined_topologies=self.combined_topologies,
                                    protein=self.protein,
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
