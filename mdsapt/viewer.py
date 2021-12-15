"""Viewer for inputted selections"""

from typing import Union

import nglview as nv
import numpy as np

from .sapt import TrajectorySAPT

import MDAnalysis as mda

from .reader import InputReader
from .optimizer import Optimizer


class Viewer(object):

    _unv: mda.Universe
    _opt: Optimizer

    def __init__(self, settings: InputReader) -> None:
        self._unv = mda.Universe(settings.top_path, settings.trj_path)
        self._opt = Optimizer(settings)

    @staticmethod
    def _launch_viewer(system: Union[mda.Universe, mda.AtomGroup], **nglview_kwargs) -> None:
        nv.show_mdanalysis(system, **nglview_kwargs)

    def view_system(self, **nglview_kwargs) -> None:
        self._launch_viewer(self._unv, **nglview_kwargs)

    def view_residue(self, resid: int, **nglview_kwargs) -> None:
        resid: mda.AtomGroup = self._unv.select_atoms(f'resid {resid}')
        self._launch_viewer(resid, **nglview_kwargs)

    def view_interaction_pair(self, resid1: int, resid2: int, **nglview_kwargs) -> None:
        r1: mda.AtomGroup = self._unv.select_atoms(f'resid {resid1}')
        r2: mda.AtomGroup = self._unv.select_atoms(f'resid {resid2}')
        r_pair: mda.AtomGroup = r1 + r2
        self._launch_viewer(r_pair, **nglview_kwargs)

    def view_optimized_residue(self, resid: int, **nglview_kwargs) -> None:
        resid = self._opt[resid]
        self._launch_viewer(resid, **nglview_kwargs)

    def view_optimized_interaction_pair(self, resid1: int, resid2: int, **nglview_kwargs) -> None:
        r1 = self._opt[resid1]
        r2 = self._opt[resid2]
        r_pair: mda.Universe = mda.Universe.empty(n_atoms=(r1.n_atoms + r2.n_atoms), trajectory=True)
        r_pair.add_TopologyAttr('masses', [x for x in r1.mases] + [x for x in r2.masses])
        r_pair.add_TopologyAttr('name', [x for x in r1.names] + [x for x in r2.names])
        r_pair.atoms.positions = np.row_stack((r1.positions + r2.positions))
        self._launch_viewer(r_pair, **nglview_kwargs)
