r"""
:mod:`mdsapt.sapt` -- Calculate SAPT energy between selections
==============================================================

Sets up and runs `SAPT <https://psicode.org/psi4manual/master/sapt.html>`_
calculations between the residues selected in the input file.

Required Input:

- :class:`-mdsapt.reader.InputReader
- :class:`-mdsapt.optimizer.Optimizer



.. autoclass:: TrajectorySAPT
    :members:
    :inherited-members:
"""

from typing import Dict, List

import pandas as pd

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.topology.guessers import guess_types

import psi4

from .reader import InputReader
from .optimizer import Optimizer


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
    results: pd.DataFrame

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
        super(TrajectorySAPT, self).__init__(self._unv.trajectory)

    def _prepare(self) -> None:
        self._col = ['residues', 'time', 'energy']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {x: [] for x in self._col}

    def _get_psi_mol(self, key: int):
        resid: mda.AtomGroup = self._opt.rebuild_resid(key, self._sel[key])
        coords: str = ''
        for atom in resid.atoms:
            coords += f'\n{atom.name[0]} {atom.position[0]} {atom.position[1]} {atom.position[2]}'
        return coords

    def _single_frame(self) -> None:
        xyz_dict = {k: self._get_psi_mol(k) for k in self._sel.keys()}
        with open('test.xyz', 'w+') as r:
            r.write(xyz_dict[2])
        for pair in self._sel_pairs:
            coords = xyz_dict[pair[0]] + '\n--\n' + xyz_dict[pair[1]] + '\nunits angstrom'
            dimer = psi4.geometry(coords)
            psi4.set_options({'scf_type': 'df',
                              'freeze_core': 'true'})
            psi4.set_memory(self._mem)

            sapt = psi4.energy('sapt0/jun-cc-pvdz', molecule=dimer)
            result = [f'{pair[0]}-{pair[1]}', self._ts.time, sapt]
            for r in range(len(result)):
                self._res_dict[self._col[r]].append(result[r])

    def _conclude(self) -> None:
        for k in self._col:
            self.results[k] = self._res_dict[k]
