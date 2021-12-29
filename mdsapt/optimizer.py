r"""
:mod:`mdsapt.optimizer` -- Prepare residues for SAPT calculations
=================================================================

Prepares residues for SAPT calculations by adding protons and replacing missing atoms

When pulled out of the peptide backbone residues are missing protons on both the C and N
terminus giving an unbalanced spin multiplicity. This causes SAPT calculations to fail.

Required Input:

- :class:`mdsapt.reader.InputReader`


.. autofunction:: get_spin_multiplicity

.. autoclass:: Optimizer
    :members:
    :inherited-members:

"""

from typing import Dict, Optional, KeysView

import MDAnalysis as mda


from MDAnalysis.converters.RDKit import atomgroup_to_mol
from MDAnalysis.topology.guessers import guess_types, guess_atom_element

import psi4

from rdkit import Chem

from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

import numpy as np

from .reader import InputReader

import logging

logger = logging.getLogger('mdsapt.optimizer')


def get_spin_multiplicity(molecule: Chem.Mol) -> int:
    """Returns the spin multiplicity of a :class:`RDKit.Mol`.
    Based on method in http://www.mayachemtools.org/docs/modules/html/code/RDKitUtil.py.html .

    :Arguments:
        *molecule*
            :class:`RDKit.Mol` object
    """
    radical_electrons: int = 0

    for atom in molecule.GetAtoms():
        radical_electrons += atom.GetNumRadicalElectrons()

    total_spin: int = radical_electrons // 2
    spin_mult: int = total_spin + 1
    return spin_mult


class Optimizer(object):
    """Prepares amino acid residues from the peptide backbone
    for sapt calculations. This requires selecting the residues
    given in :class:`mdsapt.reader.InputReader`, adding protons
    to the N terminus, then adding a proton to the C terminus and
    optimizing the bond length.

    For after the length has been determined once for a residue
    the proton can be added again without needing to preform a
    bond optimization. This is so that while running the SAPT
    calculation getting valid coordinates is faster."""

    _resids: Dict[int, mda.AtomGroup]
    _unv: mda.Universe
    _settings: InputReader
    _bond_lengths: Dict[int, float]
    _opt_set: Dict[str, str]
    _basis: str

    def __init__(self, settings: InputReader) -> None:
        """Prepares selected residues for SAPT calculations
        by adding missing protons. Obtains new bond lengths
        for protons added to C terminus of the amino acids.

        :Arguments:
            *settings*
                :class:`mdsapt.reader.InputReader`
        """
        self._settings = settings
        self._unv = mda.Universe(self._settings.top_path, self._settings.trj_path)
        self._resids = {x: self._unv.select_atoms(f"resid {x}") for x in self._settings.ag_sel}
        self._bond_lengths = {}
        self._opt_set = settings.opt_settings['settings']
        self._basis = settings.opt_settings['basis']
        self._prepare_resids()

    @property
    def num_failed_residue(self) -> int:
        return len(self._resids.keys()) - len(self._bond_lengths)

    def _prepare_resids(self) -> None:
        logger.info('Running optimiztions of residues')
        logger.info(f'Attempting optimization with {self._basis} basis, \n and {self._settings} settings')
        self._run_opts(self._resids.keys())

    def _run_opts(self, resid_list: KeysView[int]) -> None:
        for k in resid_list:  # List of residue numbers should be obtained from keys
            self._run_opt_steps(k)

    def _run_opt_steps(self, key: int) -> None:
        step0: mda.AtomGroup = self._resids[key]  # Get resid for optimization
        step1: mda.AtomGroup = self._fix_amino(step0)  # Fix amino group
        # Add proton to backbone
        step2: mda.Universe = self._protonate_backbone(step1)
        logger.info(f'Optimizing new bond for residue {key}')
        length: Optional[float] = self._opt_geometry(step2, key)
        if length is not None:  # Check if value for length obtained
            self._bond_lengths[key] = length  # Hash new bond length

    def rebuild_resid(self, key: int, resid: mda.AtomGroup) -> mda.AtomGroup:
        """Rebuilds residue by replacing missing protons and adding a new proton
        with the bond length determined by the optimization. Raises key error if class
        has no value for that optimization."""
        try:
            length = self._bond_lengths[key]
        except KeyError:
            logger.error(f'No bond length for {key}')
            raise KeyError

        step0: mda.AtomGroup = self._fix_amino(resid)
        step1: mda.Universe = self._protonate_backbone(step0, length=length)
        return step1.select_atoms("all")

    def _fix_amino(self, resid: mda.AtomGroup) -> mda.AtomGroup:
        resid.write('resid.pdb', file_format='PDB')  # Saving residue
        fixer = PDBFixer(filename='resid.pdb')
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingHydrogens(self._settings.trj_settings['pH'])  # Adding protons at pH value
        PDBFile.writeFile(fixer.topology, fixer.positions, open('resid_fixed.pdb', 'w'))

        res_fixed = mda.Universe('resid_fixed.pdb')
        resid: mda.AtomGroup = res_fixed.select_atoms("resname *")
        resid.guess_bonds()
        return resid

    @staticmethod
    def _get_new_pos(backbone: mda.AtomGroup, length: float):
        c_pos = backbone.select_atoms('name C').positions[0]
        o_pos = backbone.select_atoms('name O').positions[0]
        a_pos = backbone.select_atoms('name CA').positions[0]
        o_pos = o_pos - c_pos  # Translate coords such that C in at origin
        a_pos = a_pos - c_pos
        o_norm = o_pos/np.linalg.norm(o_pos)
        a_norm = a_pos/np.linalg.norm(a_pos)
        h_pos = -(o_norm + a_norm)
        h_norm = h_pos/np.linalg.norm(h_pos)
        h_norm = (h_norm * length) + c_pos
        return h_norm

    def _protonate_backbone(self, resid: mda.AtomGroup, length: float = 1.128) -> mda.Universe:
        backbone = resid.select_atoms('backbone')

        protonated: mda.Universe = mda.Universe.empty(n_atoms=resid.n_atoms + 1, trajectory=True)
        protonated.add_TopologyAttr('masses', [x for x in resid.masses] + [1])
        protonated.add_TopologyAttr('name', [x for x in resid.names] + ['H*'])
        protonated.add_TopologyAttr('types', guess_types(protonated.atoms.names))
        protonated.add_TopologyAttr('elements', [guess_atom_element(atom) for atom in protonated.atoms.names])
        new_pos = resid.positions
        h_pos = self._get_new_pos(backbone, length)
        protonated.atoms.positions = np.row_stack((new_pos, h_pos))
        return protonated

    def _opt_geometry(self, system: mda.Universe, key: int, basis: str = 'scf/cc-pvdz') -> Optional[float]:
        resid: mda.AtomGroup = system.select_atoms('all')
        rd_mol: Chem.Mol = atomgroup_to_mol(resid)
        coords: str = f'{Chem.GetFormalCharge(rd_mol)} {get_spin_multiplicity(rd_mol)}'
        freeze_list: str = ''
        for n in range(len(system.atoms)):
            atom = system.atoms[n]
            coords += f'\n{atom.name[0]} {atom.position[0]} {atom.position[1]} {atom.position[2]}'
            if atom.name != 'H*':
                freeze_list += f'\n{n + 1} xyz'
            elif atom.name == 'CA':
                ca_ind = n

        mol: psi4.core.Molecule = psi4.geometry(coords)
        psi4.set_memory(self._settings.sys_settings['memory'])
        psi4.set_options(self._opt_set)
        try:
            psi4.optimize(basis, freeze_list=freeze_list, opt_cooridnates='cartesian',  molecule=mol)
        except psi4.PsiException:
            logger.warning(f'Optimization failed on resid {key}')
            return None

        opt_coords = mol.create_psi4_string_from_molecule()

        for line in opt_coords:
            if 'H' in line:
                line = line.split()
                h_coord = [float(line[1]), float(line[2]), float(line[3])]

        c_line = opt_coords[ca_ind].split()
        c_coord = [float(c_line[1]), float(c_line[2]), float(c_line[3])]
        return np.sqrt((c_coord[0] - h_coord[0])**2 + (c_coord[1] - h_coord[1])**2 + (c_coord[2] - h_coord[1])**2)

    def set_opt_setting(self, opt_settings: Dict[str, str]) -> None:
        """Changes `Psi4 <https.psicode.org>`_ optimization settings.
        The optimization must be rerun for this to take effect.

        :Arguments:
            *opt_settings*
                Setting applied before optimization
            """
        self._opt_set = opt_settings

    def set_basis(self, basis: str) -> None:
        """Changes basis used in optimizations see `Psi4 <https.psicode.org>`_
        docs for information of available basis sets. Default is scf/cc-pvdz.

       :Arguments:
            *basis*
                Changes basis set used in optimization
        """
        self._basis = basis

    def rerun_optimizations(self) -> None:
        """Reruns optimization of bond lengths."""
        logger.info('Rerunning optimization')
        self._prepare_resids()

    def rerun_failed_optimizations(self) -> None:
        """Reruns failed optimizations, allowing for new
        for running optimization with new settings."""
        logger.info('Rerunning optimization of failed residues')
        logger.info(f'Attempting optimization with {self._basis} basis, \n and {self._settings} settings')
        self._run_opts(self._bond_lengths.keys())

    def rerun_residue(self, residue_id: int) -> None:
        """Reruns failed optimization on specified residue allowing
        for optimization with new settings.

        :Arguments:
            *residue_id*
                Number of the residue in the polypeptide chain"""
        logger.info(f'Rerunning optimization of failed residues {residue_id}')
        logger.info(f'Attempting optimization with {self._basis} basis, \n and {self._settings} settings')
        self._run_opt_steps(residue_id)
