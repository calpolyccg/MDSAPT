r"""
:mod:`mdsapt.optimizer` -- Prepare residues for SAPT calculations
=================================================================

Prepares residues for SAPT calculations by adding protons and replacing missing atoms

When pulled out of the peptide backbone residues are missing protons on both the C and N
terminus giving an unbalanced spin multiplicity. This causes SAPT calculations to fail.

Required Input:

- :class:`mdsapt.reader.InputReader`


.. autoclass:: Optimizer
    :members:
    :inherited-members:

"""

from typing import Dict

import MDAnalysis as mda

import psi4

from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

import numpy as np

from .reader import InputReader

import logging

logger = logging.getLogger('mdsapt')


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
    _bond_lenghts: Dict[int, np.ndarray]
    _opt_set: Dict[str, str]

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
        self._bond_lenghts = {}
        self._opt_set = {'reference': 'rhf'}
        self._prepare_resids()

    def _prepare_resids(self) -> None:
        for k in self._resids:
            step0: mda.AtomGroup = self._resids[k]  # Get resid for optimization
            step1: mda.AtomGroup = self._fix_amino(step0)  # Fix amino group
            step2: mda.Universe = self._protonate_backbone(step1)  # Add proton to backbone
            logger.info(f'Optimizing new bond for residue {k}')
            length: float = self._opt_geometry(step2)  # Get optimized new C-H bond
            self._bond_lenghts[k] = length  # Hash new bond length

    def rebuild_resid(self, key: int, resid: mda.AtomGroup) -> mda.AtomGroup:
        """Rebuilds residue by replacing missing protons and adding a new proton
        with the bond length determined by the optimization. Raises key error if class
        has no value for that optimization."""
        try:
            length = self._bond_lenghts[key]
        except KeyError:
            logger.error(f'No bond length for {key}')
            raise KeyError

        step0: mda.AtomGroup = self._fix_amino(resid)
        step1: mda.Universe = self._protonate_backbone(step0, length=length, just_backbone=False)
        return step1.select_atoms(f"resid {key}")

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

    def _protonate_backbone(self, resid: mda.AtomGroup, length: float = 1.128, just_backbone=True) -> mda.Universe:
        new_resid: mda.AtomGroup
        backbone = resid.select_atoms('backbone')
        if just_backbone:
            new_resid = backbone
        else:
            new_resid = resid

        protonated: mda.Universe = mda.Universe.empty(n_atoms=new_resid.n_atoms + 1, trajectory=True)
        protonated.add_TopologyAttr('masses', [x for x in new_resid.masses] + [1])
        protonated.add_TopologyAttr('name', [x for x in new_resid.names] + ['H'])
        new_pos = new_resid.positions
        h_pos = self._get_new_pos(backbone, length)
        protonated.atoms.positions = np.row_stack((new_pos, h_pos))
        return protonated

    def _opt_geometry(self, system: mda.Universe, basis: str = 'scf/dz') -> float:
        coords: str = ''
        freeze_list: str = ''
        for n in range(len(system.atoms)):
            atom = system.atoms[n]
            coords += f'\n{atom.name[0]} {atom.position[0]} {atom.position[1]} {atom.position[2]}'
            if atom.name != 'H':
                freeze_list += f'\n{n + 1} xyz'
            elif atom.name == 'CA':
                ca_ind = n

        mol: psi4.core.Molecule = psi4.geometry(coords)
        psi4.set_memory(self._settings.sys_settings['memory'])
        psi4.set_options(self._opt_set)
        psi4.optimize(basis, freeze_list=freeze_list, opt_cooridnates='cartesian',  molecule=mol)
        opt_coords = mol.create_psi4_string_from_molecule()

        for line in opt_coords:
            if 'H' in line:
                line = line.split()
                h_coord = [float(line[1]), float(line[2]), float(line[3])]

        c_line = opt_coords[ca_ind].split()
        c_coord = [float(c_line[1]), float(c_line[2]), float(c_line[3])]
        return np.sqrt((c_coord[0] - h_coord[0])**2 + (c_coord[1] - h_coord[1])**2 + (c_coord[2] - h_coord[1])**2)

    def set_opt_setting(self, opt_settings: Dict[str, str]) -> None:
        """Changes `Psi4 <https.psicode.org>` optimization settings.
        The optimization must be rerun for this to take effect."""
        self._opt_set = opt_settings

    def re_run_optimizations(self) -> None:
        """Reruns optimization of bond lengths."""
        logger.info('Rerunning optimization')
        self._prepare_resids()
