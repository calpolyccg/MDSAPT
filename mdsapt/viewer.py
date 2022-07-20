r"""
:mod:`mdsapt.viewer` -- Visualize residues and interaction pairs
================================================================

Allows for visualization of trajectories using `NGLView <http://nglviewer.org>`_
in a Jupyter Notebook.

.. note::
    This module only works if NGLView is installed. It is likely not
    automatically installed by your package manager because it is an optional
    dependency.

Required Input:

- :class:`-mdsapt.reader.InputReader
- :class:`-mdsapt.optimizer.Optimizer

.. autoclass:: Viewer
    :members:
    :inherited-members:
"""

from typing import Union

try:
    import nglview as nv
except ImportError:
    raise ImportError(
        "nglview is not installed! Please install it to use the viewer module."
    )
    
import numpy as np

from .sapt import TrajectorySAPT

import MDAnalysis as mda

from .config import Config
from .repair import rebuild_resid


class Viewer:

    _unv: mda.Universe

    def __init__(self, settings: Config) -> None:
        """Sets up visualizations for selected residues"""
        self._unv = mda.Universe(settings.analysis.topology, settings.analysis.trajectories)

    @staticmethod
    def _launch_viewer(system: Union[mda.Universe, mda.AtomGroup], **nglview_kwargs)\
            -> nv.NGLWidget:
        return nv.show_mdanalysis(system, **nglview_kwargs)

    def view_system(self, **nglview_kwargs) -> nv.NGLWidget:
        """Shows MD system."""
        return self._launch_viewer(self._unv, **nglview_kwargs)

    def view_residue(self, resid: int, **nglview_kwargs) -> nv.NGLWidget:
        """Shows selected residue.

        :Arguments:
            *resid*
                number of selected residue in polypeptide chain
            *nglview_kwargs*
                arguments passed to the viewer"""
        resid: mda.AtomGroup = self._unv.select_atoms(f'resid {resid} and protein')
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
        r1: mda.AtomGroup = self._unv.select_atoms(f'resid {resid1} and protein')
        r2: mda.AtomGroup = self._unv.select_atoms(f'resid {resid2} and protein')
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
        resid = rebuild_resid(resid, self._unv.select_atoms(f'resid {resid} and protein'))
        return self._launch_viewer(resid, **nglview_kwargs)

    def view_optimized_interaction_pair(self, resid1: int, resid2: int, **nglview_kwargs)\
            -> nv.NGLWidget:
        """Shows selected pair of residues after prepared for SAPT by
        :class:`mdsapt.optimizer.

        :Arguments:
            *resid*
                number of selected residue in polypeptide chain
            *nglview_kwargs*
                arguments passed to the viewer"""
        r1 = rebuild_resid(resid1, self._unv.select_atoms(f'resid {resid1} and protein'))
        r2 = rebuild_resid(resid2, self._unv.select_atoms(f'resid {resid2} and protein'))
        r_pair: mda.Universe = mda.Universe.empty(n_atoms=(r1.n_atoms + r2.n_atoms),
                                                  trajectory=True)
        r_pair.add_TopologyAttr('name', [x for x in r1.names] + [x for x in r2.names])
        r_pair.atoms.positions = np.row_stack((r1.positions, r2.positions))
        return self._launch_viewer(r_pair, **nglview_kwargs)
