r"""
:mod:`mdsapt.sapt` -- Tools for calculating SAPT energy from MD data
====================================================================

Sets up and runs `SAPT <https://psicode.org/psi4manual/master/sapt.html>`_
calculations between the residues selected in the input file.

 autoclass:: SAPT
    :members:
    :inherited-members:

.. autoclass:: TrajectorySAPT
    :members:
    :inherited-members:

"""

from typing import Dict, List, Set, Tuple, Optional, Final

import pandas as pd

import MDAnalysis as mda

from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.topology.guessers import guess_types
from MDAnalysis.converters.RDKit import atomgroup_to_mol
from MDAnalysis.lib.log import ProgressBar

import psi4

from pydantic import ValidationError

from rdkit import Chem

from .config import Config, TrajectoryAnalysisConfig, DockingAnalysisConfig, Psi4Config, SysLimitsConfig
from .repair import is_amino, rebuild_resid, get_spin_multiplicity
from .utils.ensemble import Ensemble, EnsembleAtomGroup

import logging

logger = logging.getLogger('mdsapt.sapt')


MHT_TO_KCALMOL: Final[float] = 627.509


def build_psi4_input_str(resid: int, residue: mda.AtomGroup) -> str:
    """Generates Psi4 input file the specified residue. Prepares amino acids for SAPT using
    :class:`mdsapt.optimizer.Optimizer`. Adds charge and spin multiplicity to top of cooridnates. """
    repaired_resid: mda.AtomGroup = rebuild_resid(resid, residue)
    rd_mol = atomgroup_to_mol(repaired_resid)

    coords: str = f'{Chem.GetFormalCharge(rd_mol)} {get_spin_multiplicity(rd_mol)}'
    for atom in repaired_resid.atoms:
        coords += f'\n{atom.element} {atom.position[0]} {atom.position[1]} {atom.position[2]}'
    return coords


def calc_sapt(psi4_input: str, psi4_cfg: Psi4Config, sys_cfg: SysLimitsConfig,
              outfile: Optional[str]) -> Dict[str, float]:
    """Runs SAPT on the molecules given in the input string. If `save_psi4_output` is set to true the output will
    be saved as the given filename.

        Results are returned in a dictionary with the SAPT energy broken down by type with the following keys.

        1. SAPT TOTAL ENERGY
        2. SAPT ELST ENERGY
        3. SAPT EXCH ENERGY
        4. SAPT IND ENERGY
        5. SAPT DISP ENERGY"""
    dimer = psi4.geometry(psi4_input)
    psi4.set_options(psi4_cfg.settings)
    psi4.set_memory(sys_cfg.memory)
    psi4.set_num_threads(sys_cfg.ncpus)

    if outfile is not None:
        psi4.set_output_file(outfile)  # Saves output file

    # Calculating SAPT
    psi4.energy(f'{psi4_cfg.method}/{psi4_cfg.basis}', molecule=dimer)

    result: Dict[str, float] = {'SAPT TOTAL ENERGY': psi4.variable('SAPT TOTAL ENERGY') * MHT_TO_KCALMOL,
                                'SAPT ELST ENERGY': psi4.variable('SAPT ELST ENERGY') * MHT_TO_KCALMOL,
                                'SAPT EXCH ENERGY': psi4.variable('SAPT EXCH ENERGY') * MHT_TO_KCALMOL,
                                'SAPT IND ENERGY': psi4.variable('SAPT IND ENERGY') * MHT_TO_KCALMOL,
                                'SAPT DISP ENERGY': psi4.variable('SAPT DISP ENERGY') * MHT_TO_KCALMOL}

    # Getting results

    return result


class TrajectorySAPT(AnalysisBase):
    """Handles iterating over MD trajectory frames,
    setting up SAPT calculations, and processing results.

    Results are stored in a Pandas :class:`DataFrame` following the
    `"tidy dataframe" <https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html>`_
    convention.
    """

    _cfg: Config
    _unv: mda.Universe
    _sel: Dict[int, mda.AtomGroup]
    _sel_pairs: List[Tuple[int, int]]
    results: pd.DataFrame

    def __init__(self, config: Config, **universe_kwargs) -> None:
        """Sets up Trajectory and residue selections.

        :Arguments:
            *config*
                :class:`mdsapt.reader.InputReader containing data for running calculations
            *optimizer*
                :class:`mdsapt.optimizer.Optimizer` for preparing residues by replacing missing protons
                and providing a balanced spin state.
            *universe_arguments*
                keyword arguments for loading the trajectory into a MDAnalysis :class:`Universe <MDAnalysis.core.groups.universe.Universe>`
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
        self._sel = {x: self._unv.select_atoms(f'resid {x} and not (name OH2 or name H1 or name H2)')
                     for x in ag_sel}
        self._sel_pairs = config.analysis.pairs
        AnalysisBase.__init__(self, self._unv.trajectory)

    def _prepare(self) -> None:
        self._col = ['residues', 'time', 'total', 'electrostatic',
                     'exchange', 'induction', 'dispersion']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {x: [] for x in self._col}

    def _single_frame(self) -> None:
        outfile: Optional[str] = None
        xyz_dict = {k: build_psi4_input_str(k, self._sel[k]) for k in self._sel.keys()}
        for pair in self._sel_pairs:
            coords = xyz_dict[pair[0]] + '\n--\n' + xyz_dict[pair[1]] + '\nunits angstrom'

            logger.info(f'Starting SAPT for {pair}')
            if self._cfg.analysis.output:
                outfile = f'sapt_{pair[0]}-{pair[1]}_{self._ts.time}.out'

            sapt: Dict[str, float] = calc_sapt(coords, self._cfg.psi4,
                                               self._cfg.system_limits,
                                               outfile)
            result = [f'{pair[0]}-{pair[1]}', self._ts.time] + [sapt[x] for x in
                                                                ['SAPT TOTAL ENERGY', 'SAPT ELST ENERGY',
                                                                 'SAPT EXCH ENERGY', 'SAPT IND ENERGY',
                                                                 'SAPT DISP ENERGY']]

            for r in range(len(result)):
                self._res_dict[self._col[r]].append(result[r])

    def _conclude(self) -> None:
        for k in self._col:
            self.results[k] = self._res_dict[k]


class DockingSAPT:
    """"""

    _cfg: Config
    _ens: Ensemble
    _sel: Dict[int, EnsembleAtomGroup]
    _sel_pairs: List[Tuple[int, int]]
    _key: int
    _sapt_keys: List[str] = [
        'SAPT TOTAL ENERGY',
        'SAPT ELST ENERGY',
        'SAPT EXCH ENERGY',
        'SAPT IND ENERGY',
        'SAPT DISP ENERGY'
    ]

    def __init__(self, config: Config) -> None:
        self._cfg = config

        try:
            # Ensuring config type is correct
            if not isinstance(config.analysis, DockingAnalysisConfig):
                raise ValidationError("config.analysis.type is not docking")
        except ValidationError as err:
            logger.exception(err)
            raise err

        self._sel_pairs = self._cfg.analysis.pairs

        if config.analysis.mode == DockingStructureMode.MergedLigand:
            self._ens = Ensemble(systems_dir=config.analysis.combined_topologies)
        elif config.analysis.mode == DockingStructureMode.SeparateLigand:
            self._ens = Ensemble(protein_dir=config.analysis.protein,
                                 ligands_dir=config.analysis.ligands)
        else:
            err = ValidationError("Must specify either 'protein-ligand' or 'separate-ligand' mode in config")
            logger.exception(err)
            raise err

        self._sel = {k: self._ens.select_atoms(f'resid {k} and not (name OH2 or name H1 or name H2') for k in
                     self._cfg.ag_sel}

    def _prepare(self) -> None:
        self._col = ['structure', 'time', 'total', 'electrostatic',
                     'exchange', 'induction', 'dispersion']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {x: [] for x in self._col}

    def _single_system(self) -> None:
        xyz_dict = {k: build_psi4_input_str(k, self._sel[k][self._key]) for k in self._sel.keys()}
        outfile: Optional[str] = None

        for pair in self._sel_pairs:
            coords = xyz_dict[pair[0]] + '\n--\n' + xyz_dict[pair[1]] + '\nunits angstrom'

            logger.info(f'Starting SAPT for {pair}')

            if self._cfg.analysis.output:
                outfile = f'sapt_{pair[0]}-{pair[1]}_{self._key}.out'

            sapt: Dict[str, float] = calc_sapt(coords, self._cfg.psi4, self._cfg.system_limits, outfile)
            result = [f'{pair[0]}-{pair[1]}', self._key] + [sapt[x] for x in self._sel_pairs]

            for r in range(len(result)):
                self._res_dict[self._col[r]].append(result[r])

    def _conclude(self) -> None:
        for k in self._col:
            self.results[k] = self._res_dict[k]

    def run(self):
        """Runs _single_universe on each system and _single_frame
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
        return self
