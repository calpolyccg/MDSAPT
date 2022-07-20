r"""
:mod:`mdsapt.sapt` -- Tools for calculating SAPT energy from MD data
====================================================================

Sets up and runs `SAPT <https://psicode.org/psi4manual/master/sapt.html>`_
calculations between the residues selected in the input file.

.. autofunction:: build_psi4_imput_str

.. autofunction:: calc_sapt

.. autoclass:: TrajectorySAPT
    :members:
    :inherited-members:

.. autoclass:: DockingSAPT
    :members:
    :inherited-members:

"""
from typing import Dict, List, Set, Tuple, Optional, Final, Union, Any

import logging

import pandas as pd

import MDAnalysis as mda

from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.topology.guessers import guess_types
from MDAnalysis.converters.RDKit import atomgroup_to_mol
from MDAnalysis.lib.log import ProgressBar

import psi4

from pydantic import ValidationError

from rdkit import Chem

from .config import Config, TrajectoryAnalysisConfig, DockingAnalysisConfig, \
    Psi4Config, SysLimitsConfig
from .repair import rebuild_resid, get_spin_multiplicity
from .utils.ensemble import Ensemble, EnsembleAtomGroup

logger = logging.getLogger('mdsapt.sapt')

MHT_TO_KCALMOL: Final[float] = 627.509


def build_psi4_input_str(resid: int, residue: mda.AtomGroup) -> str:
    """
    Generates Psi4 input file the specified residue. Prepares amino acids for SAPT using
    :class:`mdsapt.optimizer.Optimizer`. Adds charge and spin multiplicity to top of cooridnates.
    """
    repaired_resid: mda.AtomGroup = rebuild_resid(resid, residue)
    rd_mol = atomgroup_to_mol(repaired_resid)

    coords: str = f'{Chem.GetFormalCharge(rd_mol)} {get_spin_multiplicity(rd_mol)}'
    for atom in repaired_resid.atoms:
        coords += f'\n{atom.element} {atom.position[0]} {atom.position[1]} {atom.position[2]}'
    return coords


def calc_sapt(psi4_input: str, psi4_cfg: Psi4Config, sys_cfg: SysLimitsConfig,
              outfile: Optional[str]) -> Dict[str, float]:
    """
    Runs SAPT on the molecules given in the input string. If `save_psi4_output`
    is set to true the output will be saved as the given filename.

    Results are returned in a dictionary with the SAPT energy broken down by
    type with the following keys.

        1. SAPT TOTAL ENERGY
        2. SAPT ELST ENERGY
        3. SAPT EXCH ENERGY
        4. SAPT IND ENERGY
        5. SAPT DISP ENERGY
    """
    dimer = psi4.geometry(psi4_input)
    psi4.set_options(psi4_cfg.settings)
    psi4.set_memory(sys_cfg.memory)
    psi4.set_num_threads(sys_cfg.ncpus)

    if outfile is not None:
        psi4.set_output_file(outfile)  # Saves output file

    # Calculating SAPT
    psi4.energy(f'{psi4_cfg.method}/{psi4_cfg.basis}', molecule=dimer)

    result: Dict[str, float] = {
        'SAPT TOTAL ENERGY': psi4.variable('SAPT TOTAL ENERGY') * MHT_TO_KCALMOL,
        'SAPT ELST ENERGY': psi4.variable('SAPT ELST ENERGY') * MHT_TO_KCALMOL,
        'SAPT EXCH ENERGY': psi4.variable('SAPT EXCH ENERGY') * MHT_TO_KCALMOL,
        'SAPT IND ENERGY': psi4.variable('SAPT IND ENERGY') * MHT_TO_KCALMOL,
        'SAPT DISP ENERGY': psi4.variable('SAPT DISP ENERGY') * MHT_TO_KCALMOL
    }

    # Getting results
    return result


class TrajectorySAPT(AnalysisBase):
    """
    Handles iterating over MD trajectory frames,
    setting up SAPT calculations, and processing results.

    Results are stored in a Pandas :class:`DataFrame` following the
    `"tidy dataframe" <https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html>`_
    convention.
    """

    _cfg: Config
    _unv: mda.Universe
    _sel: Dict[int, mda.AtomGroup]
    _sel_pairs: List[Tuple[int, int]]
    _res_dict: Dict[str, List[Any]]
    results: pd.DataFrame
    _COL: Final[List[str]] = ['residues', 'time', 'total', 'electrostatic',
                              'exchange', 'induction', 'dispersion']
    _SAPT_KEYS: Final[List[str]] = ['SAPT TOTAL ENERGY', 'SAPT ELST ENERGY',
                                    'SAPT EXCH ENERGY', 'SAPT IND ENERGY',
                                    'SAPT DISP ENERGY']

    def __init__(self, config: Config, **universe_kwargs) -> None:
        """
        Sets up Trajectory and residue selections.

        :Arguments:
            *config*
                :class:`mdsapt.config.Config` containing data for running calculations
            *universe_arguments*
                keyword arguments for loading the trajectory into a MDAnalysis
                :class:`Universe <MDAnalysis.core.groups.universe.Universe>`
        """
        try:
            # Ensuring config type is correct
            if not isinstance(config.analysis, TrajectoryAnalysisConfig):
                raise ValidationError("config.analysis.type is not trajectory")
        except ValidationError as err:
            logger.exception(err)
            raise err
        self._cfg = config
        self._unv = config.analysis.create_universe(**universe_kwargs)
        elements = guess_types(self._unv.atoms.names)
        self._unv.add_TopologyAttr('elements', elements)
        ag_sel: Set[int] = config.analysis.get_selections()
        self._sel = {
            x: self._unv.select_atoms(f'resid {x} and not (name OH2 or name H1 or name H2)')
            for x in ag_sel
        }

        self._sel_pairs = config.analysis.pairs
        AnalysisBase.__init__(self, self._unv.trajectory)

    def _prepare(self) -> None:
        self.results = pd.DataFrame(columns=self._COL)
        self._res_dict = {x: [] for x in self._COL}

    def _single_frame(self) -> None:
        outfile: Optional[str] = None
        xyz_dict = {k: build_psi4_input_str(k, self._sel[k]) for k in self._sel.keys()}
        for pair in self._sel_pairs:
            coords = xyz_dict[pair[0]] + '\n--\n' + xyz_dict[pair[1]] + '\nunits angstrom'

            logger.info(f'Starting SAPT for {pair}')

            if self._cfg.psi4.save_output:
                outfile = f'sapt_{pair[0]}-{pair[1]}_{self._ts.time}.out'

            sapt: Dict[str, float] = calc_sapt(coords, self._cfg.psi4,
                                               self._cfg.system_limits,
                                               outfile)
            result = [f'{pair[0]}-{pair[1]}', self._ts.time] + [sapt[x] for x in self._SAPT_KEYS]

            for r in range(len(result)):
                self._res_dict[self._COL[r]].append(result[r])

    def _conclude(self) -> None:
        for k in self._COL:
            self.results[k] = self._res_dict[k]


class DockingSAPT:
    """
    Handles running SAPT calculations on a collection of protein-ligand topologies.

    Results are stored in a Pandas :class:`DataFrame` following the
    `"tidy dataframe" <https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html>`_
    convention.
    """

    _cfg: Config
    _ens: Ensemble
    _sel: Dict[int, EnsembleAtomGroup]
    _sel_pairs: List[Tuple[int, int]]
    _COL: Final[List[str]] = ['structure', 'pair', 'total', 'electrostatic',
                              'exchange', 'induction', 'dispersion']
    _key: str
    key_name: Dict[str, str]
    _pair_names: Dict[Tuple[int, int], str]
    _SAPT_KEYS: Final[List[str]] = [
        'SAPT TOTAL ENERGY',
        'SAPT ELST ENERGY',
        'SAPT EXCH ENERGY',
        'SAPT IND ENERGY',
        'SAPT DISP ENERGY'
    ]

    def __init__(self, config: Config) -> None:
        """
        Sets up ligand topology systems.

        :Arguments:
            *config*
                A :class:`mdsapt.config.Config` for a docking analysis.
        """
        self._cfg = config

        try:
            # Ensuring config type is correct
            if not isinstance(config.analysis, DockingAnalysisConfig):
                raise ValidationError("config.analysis.type is not docking")
        except ValidationError as err:
            logger.exception(err)
            raise err

        self._ens = self._cfg.analysis.build_ensemble()

        self._sel_pairs = self._cfg.analysis.pairs

        self._sel = {
            k: self._ens.select_atoms(f'resid {k} and not (name OH2 or name H1 or name H2)')
            for k in self._cfg.analysis.get_selections()
        }

    def _prepare(self) -> None:
        self.results = pd.DataFrame(columns=self._COL)
        self._res_dict = {x: [] for x in self._COL}
        self._key_names = {k: k.split("/")[-1].split(".")[0] for k in self._ens.keys()}
        self._pair_names = {pair: f'{pair[0]}-{pair[1]}' for pair in self._sel_pairs}

    def _single_system(self) -> None:
        xyz_dict = {k: build_psi4_input_str(k, self._sel[k][self._key]) for k in self._sel.keys()}
        outfile: Optional[str] = None

        for pair in self._sel_pairs:

            coords = xyz_dict[pair[0]] + '\n--\n' + xyz_dict[pair[1]] + '\nunits angstrom'

            logger.info(f'Starting SAPT for {pair}')

            if self._cfg.psi4.save_output:
                outfile = f'sapt_{self._key_names[self._key]}_{self._pair_names[pair]}.out'

            sapt_result = calc_sapt(coords, self._cfg.psi4, self._cfg.system_limits, outfile)
            result: List[Union[str, float]] = \
                [self._key_names[self._key], self._pair_names[pair]] + \
                [sapt_result[k] for k in self._SAPT_KEYS]

            for i, res in enumerate(result):
                self._res_dict[self._COL[i]].append(res)

    def _conclude(self) -> None:
        for k in self._COL:
            self.results[k] = self._res_dict[k]

    def run(self) -> None:
        """
        Runs _single_universe on each system and _single_frame
        on each frame in the system.
        First iterates through keys of ensemble, then runs _setup_system
        which defines the system and trajectory. Then iterates over
        trajectory frames.
        """
        logger.info("Setting up systems")
        for self._key in ProgressBar(self._ens.keys(), verbose=True):
            self._prepare()
            self._single_system()
            logger.info("Moving to next universe")
        logger.info("Finishing up")
        self._conclude()
