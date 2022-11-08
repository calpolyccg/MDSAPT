# Copied from MDPOW
# Written 2021 by Alia Lescoulie for MDPOW by Oliver Beckstein et al.
r"""
:mod"`mdsapt.utils.ensemble` -- Managing collections of MD systems
==================================================================

A set of objects for manging and analyzing collections of similar MD simulations. Based on work by Alia Lescoulie during SPIDAL REU summer 2021 for `MDPOW <mdpow.readthedocs.io>`_

.. autoclass:: Ensemble
    :members:
    :inherited-members:

.. autoclass:: EnsembleAtomGroup
    :members:
    :inherited-members:

"""
import os
from os import PathLike
from pathlib import Path
from typing import Optional, List, Tuple, Union, Dict, Iterable, Set

import MDAnalysis as mda
from MDAnalysis.core.universe import Merge
from MDAnalysis.exceptions import FileFormatWarning, NoDataError, MissingDataWarning, SelectionError
from pydantic import DirectoryPath

import numpy as np

from .utils import in_dir

import logging

logger = logging.getLogger('mdsapt.utils.ensemble')


class Ensemble:
    """Collection of related :class:`MDAnalysis.Universe <MDAnalysis.core.groups.universe.Universe>`
    objects.
    Stores systems produced by running mdpow-fep organized
    by solvent, interaction, and lambda.
    Given a mdpow simulation directory will load the MD
    simulation files with the directory structure as keys.
    :Keywords:
    *dirname*
        Molecule Simulation directory. Loads simulation files present in
        lambda directories into the new instance. With this method for
        generating an :class:`~mdpow.analysis.ensemble.Ensemble` the
        lambda directories are explored and
        :meth:`~mdpow.analysis.ensemble.Ensemble._load_universe_from_dir`
        searches for .gro, .gro.b2z, .gro.gz, and .tpr files for topology,
        and .xtc files for trajectory. It will default to using the tpr file
        available.
    *solvents*
        Solvents from directory given to the new instance. Default
        :code:`solvents=('water', 'octanol')`
    *topology_paths*
        Specifies topologies used in loading simulated systems. Given
        with a dictionary with keys-value pair for each solvent and
        its respective topology path.
    *interactions*
        Interactions from directory given to the instance. Default
        :code:`interactions=('Coulomb', 'VDW')`
    *universe_kwargs*
        `Keywords arguments <https://docs.mdanalysis.org/stable/documentation_pages/core/universe>`_
        for loading :class:`MDAnalysis.Universe <MDAnalysis.core.groups.universe.Universe>`
        objects from MDPOW files in :code:`dirname` argument directory when creating an
        :class:`~mdpow.analysis.ensemble.Ensemble` .
    .. rubric:: Examples
    Typical workflow for MDPOW directory::
        ens = Ensemble(dirname='molecule')
    Typical workflow for adding universes individually::
        ens = Ensemble()
        u = mda.Universe(md.gro', 'md.xtc')
        ens.add_system(u)
    Topology paths can be specified when defining the _ensemble
    by giving the paths to each solvent topology in a dictionary
    with the topology_paths argument::
        ens = Ensemble(dirname='molecule', topology_paths={'water': water_path,
                                                           'octanol': octanol_path}
    Interactions can also be specified when initializing the with
    a list using the interactions argument::
        ens = Ensemble(dirname='molecule', interactions=['Coulomb']
    .. versionadded:: 2.0.0
    """

    _top_types: List[str] = ['psf', 'crd', 'pdb', 'ent', 'pqr', 'pdbqt',
                             'gro', 'top', 'prmtop', 'parm7', 'dms', 'tpr', 'itp',
                             'mol2', 'data', 'lammpsdump', 'xyz', 'txyz', 'arc',
                             'gms', 'log', 'config', 'history', 'xml', 'gsd', 'mmtf',
                             'in']
    _ensemble: Dict[str, mda.Universe]

    def __init__(self, universes: Optional[Dict[str, mda.Universe]] = None):
        self._ensemble = {} if universes is None else universes

    def __repr__(self) -> str:
        return f"<Ensemble Containing {len(self)} System>"

    def __len__(self) -> int:
        return len(self._ensemble)

    def __getitem__(self, index: str) -> mda.Universe:
        """Allows dictionary like indexing"""
        return self._ensemble[index]

    def keys(self) -> Iterable[str]:
        """Returns list of system keys"""
        return self._ensemble.keys()

    def items(self) -> Iterable[Tuple[str, mda.Universe]]:
        """Returns an iterable of key/value pairs"""
        return self._ensemble.items()

    def values(self) -> Iterable[mda.Universe]:
        """Returns an iterable of values"""
        return self._ensemble.values()

    def merge(self, *args, ligand_id: int = -1) -> 'Ensemble':
        """
        Merge a list of atom group into the existing ensemble returning a new merged ensemble,
        the existing item is set to resid id -1 by default, intended for adding proteins to a ligand
        """
        _ens: Dict[str, mda.Universe] = {}

        for k, univ in self.items():
            self[k].universe.residues.resids = [ligand_id]
            _ens[k] = Merge(univ.select_atoms('name *'), *args)

        return Ensemble(_ens)

    @classmethod
    def build_from_dir(cls, ensemble_dir: DirectoryPath, **universe_kwargs) -> 'Ensemble':
        """
        Finds simulation files genderated by MDPOW and attempts to build
        :class:`MDAnalysis.Universe <MDAnalysis.core.groups.universe.Universe>`
        in the lambda directories.
        Run if :code:`dirname` argument is given when initializing the class.
        First enters FEP directory, then traverses solvent and interaction
        directories to search lambda directories for system files.
        """

        with in_dir(str(ensemble_dir), create=False):
            cur_dir = os.listdir(os.curdir)
            top = []

            for file in cur_dir:
                if any([file.endswith(x) for x in cls._top_types]):
                    # Saving topology directories
                    top.append(file)

            if len(top) == 0:
                logger.warning('No MD files detected in %s', os.curdir)
                return Ensemble()

            _ens: Dict[str, mda.Universe] = {}

            for f in top:
                try:
                    u = mda.Universe(os.path.abspath(f), **universe_kwargs)
                    _ens[f.split('.')[0]] = u
                except (ValueError, FileFormatWarning, NoDataError, MissingDataWarning, OSError) as err:
                    logger.exception('Error while loading topology %r in dir %r', top[0], cur_dir)
                    raise NoDataError from err
            return Ensemble(_ens)

    @classmethod
    def build_from_files(cls, topologies: List[Union[str, Path]], **universe_kwargs) -> 'Ensemble':
        """Constructs an ensemble from a list of files."""
        _ens: Dict[str, mda.Universe] = {}
        for path in topologies:
            name: str = str(path)
            try:
                unv = mda.Universe(name, **universe_kwargs)
                _ens[name] = unv
            except (mda.exceptions.NoDataError, OSError, ValueError) as err:
                logger.exception(err)
                raise err
        return Ensemble(_ens)

    def select_atoms(self, *args, **kwargs):
        """Returns :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` containing selections
        from the :class:`~mdpow.analysis.ensemble.Ensemble`
        Uses the same
        `selection commands <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_
        as MDAnalysis, and has the same keys as the :class:`~mdpow.analysis.ensemble.Ensemble`"""
        selections = {}
        for k, unv in self.items():
            try:
                ag = unv.select_atoms(*args, **kwargs)
            except SelectionError as err:
                logger.exception("Failed to select system %r with selection settings %r %r", k, args, kwargs)
                raise err
            else:
                selections[k] = ag
        return EnsembleAtomGroup(selections, ensemble=self)

    def select_systems(self, keys: List[str]) -> 'Ensemble':
        """
        Select specific subset of systems and returns them in an Ensemble.
        This can be accomplished in two ways, by specific keys, or by
        specifying the desired system attributes solvents, interactions and
        lambdas. All arguments are stored in list form.
        :keywords:
        *keys*
            System keys from :class:`~mdpow.analysis.ensemble.Ensemble`
            to be returned.
        """
        _ens: Dict[str, mda.Universe] = {}

        for k in keys:
            logger.info('adding system %r to ensemble', k)
            _ens[k] = self[k]
        return Ensemble(_ens)


class EnsembleAtomGroup:
    """Group for storing selections from :class:`~mdpow.analysis.ensemble.Ensemble`
     objects made using the :meth:`~mdpow.analysis.ensemble.Ensemble.select_atoms` method.
    :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` is not set up for manual initialization,
    they should be obtained by selecting atoms from an existing object.
    """

    _groups: Dict[str, mda.AtomGroup]
    _ensemble: Ensemble

    def __init__(self, group_dict: dict, ensemble: Ensemble):
        self._groups = group_dict
        self._ensemble = ensemble

    def __getitem__(self, index) -> mda.AtomGroup:
        return self._groups[index]

    def __eq__(self, other) -> bool:
        if self.keys() == other.keys():
            return all(self[k] == other[k] for k in self.keys())
        return False

    def __len__(self) -> int:
        return len(self._ensemble)

    def items(self) -> Iterable[Tuple[str, mda.AtomGroup]]:
        """Returns an iterable of key/value pairs"""
        return self._groups.items()

    def values(self) -> Iterable[mda.AtomGroup]:
        """Returns an iterable of values"""
        return self._groups.values()

    def keys(self) -> Iterable[str]:
        """List of keys to specific atom groups in the system"""
        return self._groups.keys()

    def positions(self, keys=None) -> Dict[str, np.ndarray]:
        """Returns the positions of the keys of the selected atoms.
        If no keys are specified positions for all keys are returned"""
        positions = {}
        if keys is not None:
            for k in keys:
                positions[k] = self._groups[k].positions
        else:
            for k, u in self.items():
                positions[k] = u.positions
        return positions

    def select_atoms(self, *args, **kwargs):
        """Returns :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` containing selections
        from the :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`
        Uses the same
        `selection commands <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_
        as MDAnalysis, and has the same keys as :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`"""
        selections = {}
        for key in self.items():
            try:
                ag = self[key].select_atoms(*args, **kwargs)
            except SelectionError as err:
                logger.exception("Failed to select in system %r with selection traj_settings %r %r", key, args, kwargs)
                raise
            else:
                selections[key] = ag
        return EnsembleAtomGroup(selections, ensemble=self._ensemble)

    @property
    def ensemble(self) -> Ensemble:
        """Returns the ensemble of the EnsembleAtomGroup"""
        return self._ensemble
