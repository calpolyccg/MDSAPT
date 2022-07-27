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

from typing import Set, Union

import logging

import numpy as np

import MDAnalysis as mda
from MDAnalysis.converters.RDKit import atomgroup_to_mol
from MDAnalysis.topology.guessers import guess_types, guess_atom_element

from rdkit import Chem

from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

logger = logging.getLogger('mdsapt.optimizer')


def get_spin_multiplicity(molecule: Chem.Mol) -> int:
    """Returns the spin multiplicity of a :class:`RDKIT.Mol`.
    Based on method in http://www.mayachemtools.org/docs/modules/html/code/RDKitUtil.py.html .

    :Arguments:
        *molecule*
            :class:`RDKIT.Mol` object
    """
    radical_electrons: int = 0

    for atom in molecule.GetAtoms():
        radical_electrons += atom.GetNumRadicalElectrons()

    total_spin: int = radical_electrons // 2
    spin_mult: int = total_spin + 1
    return spin_mult


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
