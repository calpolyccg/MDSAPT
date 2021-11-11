"""Viewer for inputted selections"""

import nglview as nv

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

    def view_system(self, **nglview_kwargs) -> None:
        nv.show_mdanalysis(self._unv, **nglview_kwargs)

    def view_residue(self, resid: int, **nglview_kwargs) -> None:
        resid: mda.AtomGroup = self._opt[resid]
        nv.show_mdanalysis(resid, **nglview_kwargs)

    def view_interaction_pair(self, resid1: int, resid2: int, **nglview_kwargs) -> None:
        r1: mda.AtomGroup = self._opt[resid1]
        r2: mda.AtomGroup = self._opt[resid2]
        r_pair: mda.AtomGroup = r1 + r2
        nv.show_mdanalysis(r_pair, **nglview_kwargs)

