"""Provide the primary functions."""

from typing import List, Optional
import os

from MDAnalysis.analysis.base import AnalysisBase
import MDAnalysis as mda
import psi4

from .reader import CoordsReader


class TrajectorySAPT(AnalysisBase):
    """"""
    def __init__(self, universe: mda.Universe, selections: List[mda.AtomGroup]):
        self._unv: mda.Universe = universe
        self._sel: List[mda.AtomGroup] = selections
        super(TrajectorySAPT, self).__init__(universe.trajectory)

    def _prepare(self):
        self._crd: CoordsReader = CoordsReader(self._unv, self._sel)
        self._eng: List = []

    def _single_frame(self):



