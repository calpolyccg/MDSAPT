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
import errno
from typing import Optional, List

import numpy as np

import MDAnalysis as mda
from MDAnalysis.lib.log import ProgressBar
from MDAnalysis.exceptions import FileFormatWarning, NoDataError, MissingDataWarning, SelectionError

from utils import in_dir

import logging

logger = logging.getLogger('mdsapt.utils')

class Ensemble(object):
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
    .. versionadded:: 0.8.0
    """

    _top_types: List[str] = ['psf', 'crd', 'pdb', 'ent', 'pqr', 'pdbqt',
                             'gro', 'top', 'prmtop', 'parm7', 'dms', 'tpr', 'itp',
                             'mol2', 'data', 'lammpsdump', 'xyz', 'txyz', 'arc',
                             'gms', 'log', 'config', 'history', 'xml', 'gsd', 'mmtf',
                             'in']

    def __init__(self, dirname=None, **universe_kwargs):
        self._num_systems = 0
        self._ensemble = {}
        self._keys = []
        self.unv_kwargs = universe_kwargs

        if dirname is None:
            return

        if not os.path.exists(dirname):
            logger.error(f"Directory {dirname} does not exist")
            raise FileNotFoundError(errno.ENOENT, 'Directory does not'
                                                  'exist', dirname)

        self._ensemble_dir = dirname
        self._build_ensemble(**universe_kwargs)

    def __repr__(self):
        return f"<Ensemble Containing {self._num_systems} System>"

    def __len__(self):
        return self._num_systems

    def __getitem__(self, index):
        """Allows dictionary like indexing"""
        return self._ensemble[index]

    def keys(self):
        """Returns list of system keys"""
        return self._keys

    def _build_ensemble(self, **universe_kwargs) -> None:
        """Finds simulation files genderated by MDPOW and attempts to build
        :class:`MDAnalysis.Universe <MDAnalysis.core.groups.universe.Universe>`
        in the lambda directories.
        Run if :code:`dirname` argument is given when initializing the class.
        First enters FEP directory, then traverses solvent and interaction
        directories to search lambda directories for system files."""

        with in_dir(self._ensemble_dir, create=False):
            cur_dir = os.listdir(os.curdir)
            top = []

            for file in cur_dir:
                if any([file.endswith(x) for x in self._top_types]):
                    # Saving topology directories
                    top.append(file)

            if len(top) == 0:
                logger.warning('No MD files detected in %s', os.curdir)
                return

            for f in top:
                try:
                    u = mda.Universe(os.path.abspath(f), **universe_kwargs)
                    self.add_system(f.split('.')[0], u)
                except (ValueError, FileFormatWarning, NoDataError, MissingDataWarning, OSError) as err:
                    logger.error(f'{err} raised while loading {top[0]} in dir {cur_dir}')
                    raise NoDataError

    def add_system(self, key, universe: mda.Universe):
        """Adds system from universe object for trajectory and topology files
        Existing mda.Universe object or trajectory and topology path. Ensure
        that paths are set to absolute when creating the universe."""
        self._ensemble[key] = universe
        self._keys.append(key)
        self._num_systems += 1

    def pop(self, key):
        """Removes and returns system at specified key.
        Logs if KeyError is raised."""
        system = self._ensemble.pop(key)
        self._num_systems -= 1
        return system

    def select_atoms(self, *args, **kwargs):
        """Returns :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` containing selections
        from the :class:`~mdpow.analysis.ensemble.Ensemble`
        Uses the same
        `selection commands <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_
        as MDAnalysis, and has the same keys as the :class:`~mdpow.analysis.ensemble.Ensemble`"""
        selections = {}
        for key in self.keys():
            try:
                ag = self[key].select_atoms(*args, **kwargs)
            except SelectionError as err:
                logger.error("%r on system %r with selection settings %r %r", err, key, args, kwargs)
                raise
            else:
                selections[key] = ag
        return EnsembleAtomGroup(selections, ensemble=self)

    def select_systems(self, keys: List[str]):
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
        new_ens = Ensemble()
        for k in keys:
            logger.info('adding system %r to ensemble', k)
            new_ens.add_system(k, universe=self[k])
        new_ens._ensemble_dir = self._ensemble_dir
        return new_ens


class EnsembleAtomGroup(object):
    """Group for storing selections from :class:`~mdpow.analysis.ensemble.Ensemble`
     objects made using the :meth:`~mdpow.analysis.ensemble.Ensemble.select_atoms` method.
    :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` is not set up for manual initialization,
    they should be obtained by selecting atoms from an existing object.
    """

    def __init__(self, group_dict: dict, ensemble: Ensemble):
        self._groups = group_dict
        self._ensemble = ensemble
        self._keys = group_dict.keys()

    def __getitem__(self, index):
        return self._groups[index]

    def __eq__(self, other):
        if self.keys() == other.keys():
            return all(self[k] == other[k] for k in self.keys())
        return False

    def __len__(self):
        return len(self.keys())

    def keys(self):
        """List of keys to specific atom groups in the system"""
        return self._keys

    def positions(self, keys=None):
        """Returns the positions of the keys of the selected atoms.
        If no keys are specified positions for all keys are returned"""
        positions = {}
        if not keys is None:
            for k in keys:
                positions[k] = self._groups[k].positions
        else:
            for k in self.keys():
                positions[k] = self[k].positions
        return positions

    def select_atoms(self, *args, **kwargs):
        """Returns :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` containing selections
        from the :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`
        Uses the same
        `selection commands <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_
        as MDAnalysis, and has the same keys as :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`"""
        selections = {}
        for key in self.keys():
            try:
                ag = self[key].select_atoms(*args, **kwargs)
            except SelectionError as err:
                logger.error("%r on system %r with selection settings %r %r", err, key, args, kwargs)
                raise
            else:
                selections[key] = ag
        return EnsembleAtomGroup(selections, ensemble=self._ensemble)

    @property
    def ensemble(self):
        """Returns the ensemble of the EnsembleAtomGroup"""
        return self._ensemble


class EnsembleAnalysis(object):
    """Base class for running multi-system analyses
    The class is designed based on the `AnalysisBase
    class in MDAnalysis <https://docs.mdanalysis.org/stable/documentation_pages/analysis/base.html>`_
    and is a template for creating multi-universe multi-frame analyses using the
    :class:`~mdpow.analysis.ensemble.Ensemble` object
    :Keywords:
    *ensemble*
        The :class:`~mdpow.analysis.ensemble.Ensemble` being analyzed in the class
    .. rubric:: Example Analysis
    Dihedral Analysis Demonstration::
        class DihedralAnalysis(mdpow.ensemble.EnsembleAnalysis):
            def __init__(self, DihedralEnsembleGroup):
                super(DihedralAnalysis, self).__init__(DihedralEnsembleGroup.ensemble)
                self._sel = DihedralEnsembleGroup
            def _prepare_ensemble(self):
                self.result_dict = {}
                for s in ['water', 'octanol']:
                    self.result_dict[s] = {'Coulomb': {},
                                           'VDW': {}}
                for key in self._sel.group_keys():
                    self.result_dict[key[0]][key[1]][key[2]] = None
            def _prepare_universe(self):
                self.angle_dict = {'angle': None,
                                   'time': None}
                self.angles = []
            def _single_frame(self):
                angle = calc_dihedrals(self._sel[self._key].positions[0], self._sel[self._key].positions[1],
                                       self._sel[self._key].positions[2], self._sel[self._key].positions[3])
                self.angles.append(angle)
            def _conclude_universe(self):
                self.angle_dict['time'] = self.times
                self.angle_dict['angle'] = self.angles
                self.result_dict[self._key[0]][self._key[1]][self._key[2]] = self.angle_dict
            def _conclude_ensemble(self):
                self.results = pd.DataFrame(data=self.result_dict)
        D = DihedralAnalysis.run(start=0 stop=10, step=1)
    """

    def __init__(self, ensemble=None):
        self._ensemble = ensemble

    def _setup_system(self, key, start=None, stop=None, step=None):
        self._system = self._ensemble[key]
        self._key = key
        self.start = start
        self.stop = stop
        self.step = step
        self._setup_frames(self._system.trajectory)

    def _setup_frames(self, trajectory):
        self._trajectory = trajectory
        start, stop, step = trajectory.check_slice_indices(self.start, self.stop, self.step)
        self.n_frames = len(range(start, stop, step))
        self.frames = np.zeros(self.n_frames, dtype=int)
        self.times = np.zeros(self.n_frames)

    def _single_universe(self):
        """Calculations on a single Universe object.
            Run on each universe in the ensemble during when
            self.run in called.
        """
        pass  # pragma: no cover

    def _single_frame(self):
        """Calculate data from a single frame of trajectory
        Called on each frame for universes in the Ensemble.
        """
        pass  # pragma: no cover

    def _prepare_ensemble(self):
        """For establishing data structures used in running
        analysis on the entire ensemble.
        Data structures will not be overwritten upon moving to
        next system in ensemble.
        """
        pass  # pragma: no cover

    def _prepare_universe(self):
        """For establishing data structures used in running
        analysis on each trajectory in ensemble
        Data structures will be overwritten between upon after
        each trajectory has been run
        """
        pass  # pragma: no cover

    def _conclude_universe(self):
        """Run after each trajectory is finished"""
        pass  # pragma: no cover

    def _conclude_ensemble(self):
        """Run after all trajectories in ensemble are finished"""
        pass  # pragma: no cover

    def run(self, start=None, stop=None, step=None):
        """Runs _single_universe on each system and _single_frame
        on each frame in the system.
        First iterates through keys of ensemble, then runs _setup_system
        which defines the system and trajectory. Then iterates over
        trajectory frames.
        """
        logger.info("Setting up systems")
        self._prepare_ensemble()
        for self._key in ProgressBar(self._ensemble.keys(), verbose=True):
            self._setup_system(self._key, start=start, stop=stop, step=step)
            self._prepare_universe()
            self._single_universe()
            for i, ts in enumerate(ProgressBar(self._trajectory[self.start:self.stop:self.step], verbose=True,
                                               postfix=f'running system {self._key}')):
                self._frame_index = i
                self._ts = ts
                self.frames[i] = ts.frame
                self.times[i] = ts.time
                self._single_frame()
            self._conclude_universe()
            logger.info("Moving to next universe")
        logger.info("Finishing up")
        self._conclude_ensemble()
        return self

    @staticmethod
    def check_groups_from_common_ensemble(groups: List[EnsembleAtomGroup]):
        """Checks if inputted list of
        :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` originate from
        the same :class:`~mdpow.analysis.ensemble.Ensemble`
        Checks every :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup` in
        list to determine if their :meth:`ensemble` references the same object
        in memory. If two :class:`~mdpow.analysis.ensemble.EnsembleAtomGroup`
        object don't have a common :class:`~mdpow.analysis.ensemble.Ensemble`
        :class:`ValueError` is raised."""
        for i in range((len(groups) - 1) // 2):
            # Checking if EnsembleAtomGroup.ensemble references same object in memory
            if groups[i].ensemble is not groups[-1 - i]:
                msg = '''Dihedral selections from different Ensembles,
                          ensure that all EnsembleAtomGroups are created
                          from the same Ensemble.'''
                logger.error(msg)
                raise ValueError(msg)
