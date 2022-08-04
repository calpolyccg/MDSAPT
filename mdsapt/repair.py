r"""
:mod:`mdsapt.optimizer` -- Prepare residues for SAPT calculations
=================================================================

Prepares residues for SAPT calculations by adding protons and replacing missing atoms

When pulled out of the peptide backbone residues are missing protons on both the C and N
terminus giving an unbalanced spin multiplicity. This causes SAPT calculations to fail.

Required Input:

- :class:`mdsapt.reader.InputReader`


.. autofunction:: get_spin_multiplicity

.. autofunction:: is_amino

.. autofunction:: rebuild_resid

"""

from abc import ABC
from typing import Dict, NamedTuple, Optional, Set, Tuple, Union

import logging

import numpy as np

import MDAnalysis as mda
from MDAnalysis.core.groups import Atom
from MDAnalysis.topology.guessers import guess_types, guess_atom_element

from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

from mdsapt.utils.formal_charge import ElectronInfo, calculate_electron_info, calculate_spin_multiplicity

logger = logging.getLogger('mdsapt.optimizer')


class MoleculeElectronInfo(NamedTuple):
    charge: int
    radical: int


class ChargeStrategy(ABC):
    def calculate(self, ag: mda.AtomGroup, charge_overrides: Optional[Dict[int, int]] = None) -> MoleculeElectronInfo:
        """For the given atom group, returns a pair of (formal charge, number of radical electrons)"""
        raise NotImplemented


class StandardChargeGuesser(ChargeStrategy):
    def calculate(self, ag: mda.AtomGroup, charge_overrides: Optional[Dict[int, int]] = None) -> MoleculeElectronInfo:
        if charge_overrides is None:
            charge_overrides = {}

        results = [StandardChargeGuesser._calc_single_atom(a, charge_overrides.get(a.idx)) for a in ag]

        fc = sum((pair[0] for pair in results))
        total_radicals = sum((pair[1] for pair in results))

        return MoleculeElectronInfo(fc, calculate_spin_multiplicity(total_radicals))

    @staticmethod
    def _calc_single_atom(atom: Atom, charge_override: Optional[int] = None) -> MoleculeElectronInfo:
        bonds = atom.get_connections('bonds', outside=True)
        orders = [b.order for b in bonds]

        # Assume any atom with an aromatic bond has fc=0 and radicals=0
        if 1.5 in bonds or 'ar' in bonds:
            return 0, 0

        bond_count: int = sum(orders)
        info: ElectronInfo = calculate_electron_info(atom.element, bond_count, charge_override)
        return info.fc, info.radical


def is_amino(unv: mda.Universe, resid: int) -> bool:
    std_resids: Set[str] = {
        'ALA',
        'ARG',
        'ASN',
        'ASP',
        'CYS',
        'GLU',
        'GLN',
        'GLY',
        'HIS',
        'ILE',
        'LUE',
        'LYS',
        'MET',
        'PHE',
        'PRO',
        'SER',
        'THR',
        'TRP',
        'TYR',
        'VAL'
    }

    resname_atr = unv._topology.resnames
    return resname_atr.values[resid - 1] in std_resids


def rebuild_resid(resid: int, residue: mda.AtomGroup, sim_ph: float = 7.0) -> mda.AtomGroup:
    """Rebuilds residue by replacing missing protons and adding a new proton
     on the C terminus. Raises key error if class
    has no value for that optimization."""

    def fix_amino(amino: mda.AtomGroup, sys_ph: float = 7.0) -> mda.AtomGroup:
        amino.write('resid.pdb', file_format='PDB')  # Saving residue
        fixer = PDBFixer(filename='resid.pdb')
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingHydrogens(sys_ph)  # Adding protons at pH value

        with open('resid_fixed.pdb', 'w', encoding='utf8') as file:
            PDBFile.writeFile(fixer.topology, fixer.positions, file)

        res_fixed = mda.Universe('resid_fixed.pdb')
        amino: mda.AtomGroup = res_fixed.select_atoms("resname *")
        amino.guess_bonds()
        return amino

    def get_new_pos(backbone: mda.AtomGroup, length: float) -> np.ndarray:
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

    def protonate_backbone(bkbone: mda.AtomGroup, length: float = 1.128) -> \
            Union[mda.AtomGroup, mda.Universe]:
        mol_resid = atomgroup_to_mol(bkbone)
        i: int = 0
        for atom in mol_resid.GetAtoms():
            i += atom.GetNumRadicalElectrons()

        if i > 0:
            backbone = bkbone.select_atoms('backbone')
            protonated: mda.Universe = mda.Universe.empty(n_atoms=bkbone.n_atoms + 1,
                                                          trajectory=True)
            protonated.add_TopologyAttr('masses', np.append(bkbone.masses, [1]))
            protonated.add_TopologyAttr('name', np.append(bkbone.names, ['Hc']))
            protonated.add_TopologyAttr('types', guess_types(protonated.atoms.names))
            protonated.add_TopologyAttr('elements', [guess_atom_element(atom) for
                                                     atom in protonated.atoms.names])
            new_pos = bkbone.positions
            h_pos = get_new_pos(backbone, length)
            protonated.atoms.positions = np.row_stack((new_pos, h_pos))
            return protonated
        return bkbone

    if is_amino(residue.universe, resid):
        step0: mda.AtomGroup = fix_amino(residue, sim_ph)
        step1: mda.Universe = protonate_backbone(step0, length=1.128)
        return step1.select_atoms("all")

    return residue
