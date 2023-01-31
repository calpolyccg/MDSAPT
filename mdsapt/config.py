r"""
:mod:`mdsapt.config` -- Reads input file and saves configuration
================================================================

In MD-SAPT, :obj:`mdsapt.config.Config` objects store user settings, which are needed for running analyses.
The config system was designed with a type driven development approach with the main object,
:obj:`mdsapt.config.Config`, containing a set objects each holding the settings for a portion of the overall MD-SAPT
process.
This design, along with `Pydantic <https://docs.pydantic.dev>`_ checks on each object, ensures that
any instance of :obj:`mdsapt.config.Config` is valid.
If any issues are detected, such as an error in the file paths, will raise an Error explaining that issue.

Summary of Config
_________________

:obj:`mdsapt.config.Config` contains:

- :obj:`mdsapt.config.Psi4Config`
- :obj:`mdsapt.config.SimulationConfig`
- :obj:`mdsapt.config.SystemLimitsConfig`
- :obj:`mdsapt.config.TrajectoryAnalysisConfig` or :obj:`mdsapt.config.DockingAnalysisConfig`

In most cases user configurations are provided with a yaml input file.
The input file is read using the :obj:`mdsapt.config.load_yaml_file` function, which returns a
:obj:`mdsapt.config.Config`.

.. autofunction:: load_from_yaml_file

Writing Input Files
___________________

Input files allow MD-SAPT to be used without writing code.
They contain the information needed for running analyses, and are written in a yaml format.
There are two similar but distinct input formats, one each of MD-SAPT's two analysis methods.
An example of each is provided below, and blank versions can be generated
#TODO: figure out how to make input

Example Trajectory Input
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    psi4:
      method: 'sapt0'
      basis: 'jun-cc-pvdz'
      save_output: true
      settings:
        reference: 'rhf'
    simulation:
      ph: 7.0
      charge_guesser: 'standard'
      # charge_guesser: 'rdkit'  # to use rdkit. Make sure it is installed first.
    system_limits:
      ncpus: 5
      memory: '1GB'
    analysis:
      ### This section is for running TrajectorySAPT. To run other types of analyses, see below.
      type: 'trajectory'

      topology: 'mdsapt/tests/testing_resources/testtop.psf'
      # If you want to override the charges of specific, you may specify it this way.
      # topology:
      #    path: 'your/file/path.pdb'
      #    charge_overrides:
      #      132: -2

      trajectories:
        - 'mdsapt/tests/testing_resources/testtraj.dcd'
      pairs:
        - [11, 119]
      frames:
        start: 1
        stop: 2
        #step: 3

      output: 'output.csv'

Example Docking Input
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    psi4:
      method: 'sapt0'
      basis: 'jun-cc-pvdz'
      save_output: true
      settings:
        reference: 'rhf'
    simulation:
      ph: 7.0
      charge_guesser: 'standard'
      # charge_guesser: 'rdkit'  # to use rdkit. Make sure it is installed first.
    system_limits:
      ncpus: 5
      memory: '1GB'
    analysis:
      ### This section is for running TrajectorySAPT. To run other types of analyses, see below.
      type: 'docking'
      combined_topologies: 'mdsapt/tests/testing_resources/docking_merged_test'
      pairs:
        - [11, 119]
      output: 'dock_out.csv'

Documentation for Configuration Objects
_______________________________________

.. autoclass:: Config
    :members:
    :inherited-members:

.. autoclass:: Psi4Config
    :members:
    :inherited-members:

.. autoclass:: SimulationConfig
    :members:
    :inherited-members:

.. autoclass:: SystemLimitsConfig
    :members:
    :inherited-members:

.. autoclass:: TrajectoryAnalysisConfig
    :members:
    :inherited-members:

.. autoclass:: DockingAnalysisConfig
    :members:
    :inherited-members:

"""

# There's lots of implicit class methods because pydantic decorators are stupid.
# Thus, we will disable this lint for this file.
# pylint: disable=no-self-argument

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


DockingElement: type = Union[Literal['L'], conint(ge=-1)]
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

        for v in ens.values():
            missing_selections += get_invalid_residue_selections(protein_selections, v)

        if len(missing_selections) > 0:
            errors.append(f'Selected residues are missing from topology: {missing_selections}')

        return values

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
    It can read both docking and trajectory inputs.
    """
    with Path(path).open('r', encoding='utf8') as file:
        try:
            return Config(**yaml.safe_load(file))
        except ValidationError as err:
            logger.exception("Error while loading config from %r", path)
            raise err
