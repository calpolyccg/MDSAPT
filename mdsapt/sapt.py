r"""
:mod:`mdsapt.sapt` -- Calculate SAPT energy between selections
==============================================================

Sets up and runs `SAPT <https://psicode.org/psi4manual/master/sapt.html>`_
calculations between the residues selected in the input file.

Required Input:

- :class:`-mdsapt.reader.InputReader`
- :class:`-mdsapt.optimizer.Optimizer`

.. autoclass:: TrajectorySAPT
    :members:
    :inherited-members:
"""

from typing import Dict, List

import pandas as pd

import MDAnalysis as mda

from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.topology.guessers import guess_types
from MDAnalysis.converters.RDKit import atomgroup_to_mol

import psi4

from rdkit import Chem

from .reader import InputReader
from .optimizer import Optimizer, get_spin_multiplicity

import logging

logger = logging.getLogger('mdsapt.sapt')


class TrajectorySAPT(AnalysisBase):
    """Handles iterating over MD trajectory frames,
    setting up SAPT calculations, and processing results.

    Results are stored in a Pandas :class:`DataFrame` following the
    `"tidy dataframe" <https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html>`_
    convention.
    """

    _unv: mda.Universe
    _sel: Dict[int, mda.AtomGroup]
    _sel_pairs: List[List[int]]
    _mem: str
    _cfg: InputReader
    _opt: Optimizer
    _save_psi_out: bool
    _method: str
    _basis: str
    _settings: Dict[str, str]
    results: pd.DataFrame
    _mht_to_kcalmol: float = 627529

    def __init__(self, config: InputReader, optimizer: Optimizer, **universe_kwargs) -> None:
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

        self._unv = mda.Universe(config.top_path, config.trj_path, **universe_kwargs)
        elements = guess_types(self._unv.atoms.names)
        self._unv.add_TopologyAttr('elements', elements)
        self._sel = {x: self._unv.select_atoms(f'resid {x}') for x in config.ag_sel}
        self._sel_pairs = config.ag_pair
        self._mem = config.sys_settings['memory']
        self._cfg = config
        self._opt = optimizer
        self._save_psi_out = config.sapt_out
        self._method = config.sapt_method
        self._basis = config.sapt_basis
        self._settings = config.sapt_settings['settings']

        super(TrajectorySAPT, self).__init__(self._unv.trajectory)

    def _prepare(self) -> None:
        self._col = ['residues', 'time', 'energy']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {x: [] for x in self._col}

    def _get_psi_mol(self, key: int):
        resid: mda.AtomGroup = self._opt.rebuild_resid(key, self._sel[key])
        rd_mol = atomgroup_to_mol(resid)

        coords: str = f'{Chem.GetFormalCharge(rd_mol)} {get_spin_multiplicity(rd_mol)}'
        for atom in resid.atoms:
            coords += f'\n{atom.name[0]} {atom.position[0]} {atom.position[1]} {atom.position[2]}'
        return coords

    def _single_frame(self) -> None:
        xyz_dict = {k: self._get_psi_mol(k) for k in self._sel.keys()}
        for pair in self._sel_pairs:
            coords = xyz_dict[pair[0]] + '\n--\n' + xyz_dict[pair[1]] + '\nunits angstrom'
            dimer = psi4.geometry(coords)
            psi4.set_options(self._settings)
            psi4.set_memory(self._mem)

            logger.info(f'Starting SAPT for {pair}')

            if self._save_psi_out:
                psi4.set_output_file(f'sapt_{pair[0]}-{pair[1]}_{self._ts.time}.out')  # Saves output file

            psi4.energy(f'{self._method}/{self._basis}', molecule=dimer)
            sapt = psi4.variable('SAPT TOTAL ENERGY')
            result = [f'{pair[0]}-{pair[1]}', self._ts.time, sapt*self._mht_to_kcalmol]
            for r in range(len(result)):
                self._res_dict[self._col[r]].append(result[r])

    def _conclude(self) -> None:
        for k in self._col:
            self.results[k] = self._res_dict[k]
