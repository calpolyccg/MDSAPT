r"""
:mod:`mdsapt.viewer` -- Visualize residues and interaction pairs
================================================================

Allows for visualization of trajectories using `NGLView <http://nglviewer.org>`_
in a Jupyter Notebook.

Required Input:

- :class:`-mdsapt.reader.InputReader
- :class:`-mdsapt.optimizer.Optimizer

.. autoclass:: Viewer
    :members:
    :inherited-members:
"""

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
        """Sets up visualizations for selected residues"""
        self._unv = mda.Universe(settings.top_path, settings.trj_path)
        self._opt = Optimizer(settings)

    @staticmethod
    def _launch_viewer(system: Union[mda.Universe, mda.AtomGroup], **nglview_kwargs) -> None:
        return nv.show_mdanalysis(system, **nglview_kwargs)

    def view_system(self, **nglview_kwargs) -> nv.NGLWidget:
        self._launch_viewer(self._unv, **nglview_kwargs)

    def view_residue(self, resid: int, **nglview_kwargs) -> nv.NGLWidget:
        """Shows selected residue.

        :Arguments:
            *resid*
                number of selected residue in polypeptide chain
            *nglview_kwargs*
                arguments passed to the viewer"""
        resid: mda.AtomGroup = self._unv.select_atoms(f'resid {resid}')
        return self._launch_viewer(resid, **nglview_kwargs)

    def view_interaction_pair(self, resid1: int, resid2: int, **nglview_kwargs) -> nv.NGLWidget:
        """Shows selected pair of residues.

        :Arguments:
            *resid1*
                number of first selected residue in polypeptide chain
            *resid2*
                number of second selected residue in polypeptide chain
            *nglview_kwargs*
                arguments passed to the viewer"""
        r1: mda.AtomGroup = self._unv.select_atoms(f'resid {resid1}')
        r2: mda.AtomGroup = self._unv.select_atoms(f'resid {resid2}')
        r_pair: mda.AtomGroup = r1 + r2
        return self._launch_viewer(r_pair, **nglview_kwargs)

    def view_optimized_residue(self, resid: int, **nglview_kwargs) -> nv.NGLWidget:
        """Shows selected residue after prepared for SAPT by
        :class:`mdsapt.optimizer.

        :Arguments:
            *resid*
                number of selected residue in polypeptide chain
            *nglview_kwargs*
                arguments passed to the viewer"""
        resid = self._opt.rebuild_resid(resid, self._unv.select_atoms(f'resid {resid}'))
        return self._launch_viewer(resid, **nglview_kwargs)

    def view_optimized_interaction_pair(self, resid1: int, resid2: int, **nglview_kwargs) -> nv.NGLWidget:
        """Shows selected pair of residues after prepared for SAPT by
        :class:`mdsapt.optimizer.

        :Arguments:
            *resid*
                number of selected residue in polypeptide chain
            *nglview_kwargs*
                arguments passed to the viewer"""
        r1 = self._opt.rebuild_resid(resid1, self._unv.select_atoms(f'resid {resid1}'))
        r2 = self._opt.rebuild_resid(resid2, self._unv.select_atoms(f'resid {resid2}'))
        r_pair: mda.Universe = mda.Universe.empty(n_atoms=(r1.n_atoms + r2.n_atoms), trajectory=True)
        r_pair.add_TopologyAttr('masses', [x for x in r1.mases] + [x for x in r2.masses])
        r_pair.add_TopologyAttr('name', [x for x in r1.names] + [x for x in r2.names])
        r_pair.atoms.positions = np.row_stack((r1.positions + r2.positions))
        return self._launch_viewer(r_pair, **nglview_kwargs)
