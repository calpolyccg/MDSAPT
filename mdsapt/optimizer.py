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

from typing import Dict, List

import MDAnalysis as mda


from MDAnalysis.converters.RDKit import atomgroup_to_mol
from MDAnalysis.topology.guessers import guess_types, guess_atom_element

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
    to the N terminus, then adding a proton to the C terminus.

    For the standard 20 amino acids they are prepared via the
    aforementioned process. Selected molecules besides amino
    acids will pass through without modification.
    """
    _resids: Dict[int, mda.AtomGroup]
    _unv: mda.Universe
    _settings: InputReader
    _bond_lengths: Dict[int, float] = {
        'ALA': 1.1,
        'ARG': 1.1,
        'ASN': 1.1,
        'ASP': 1.1,
        'CYS': 1.1,
        'GLU': 1.1,
        'GLN': 1.1,
        'GLY': 1.1,
        'HIS': 1.1,
        'ILE': 1.1,
        'LUE': 1.1,
        'LYS': 1.1,
        'MET': 1.1,
        'PHE': 1.1,
        'PRO': 1.1,
        'SER': 1.1,
        'THR': 1.1,
        'TRP': 1.1,
        'TYR': 1.1,
        'VAL': 1.1
    }

    _std_resids: List[str] = [x for x in _bond_lengths.keys()]

    def __init__(self, settings: InputReader) -> None:
        """Prepares selected residues for SAPT calculations
        by adding missing protons.

        :Arguments:
            *settings*
                :class:`mdsapt.reader.InputReader`
        """
        self._settings = settings
        self._unv = mda.Universe(self._settings.top_path, self._settings.trj_path)
        self._resids = {x: self._unv.select_atoms(f"resid {x}") for x in self._settings.ag_sel}

    def _is_amino(self, key: int) -> bool:
        resname_atr = self._resids[key].universe._topology.resnames
        return resname_atr.values[key - 1] in self._std_resids

    def rebuild_resid(self, key: int, resid: mda.AtomGroup) -> mda.AtomGroup:
        """Rebuilds residue by replacing missing protons and adding a new proton
         on the C terminus. Raises key error if class
        has no value for that optimization."""
        if self._is_amino(key):
            resname_atr = self._resids[key].universe._topology.resnames
            resname = resname_atr.values[key - 1]
            step0: mda.AtomGroup = self._fix_amino(resid)
            step1: mda.Universe = self._protonate_backbone(step0, length=self._bond_lengths[resname])
            return step1.select_atoms("all")
        else:
            return resid

    def _fix_amino(self, resid: mda.AtomGroup) -> mda.AtomGroup:
        resid.write('resid.pdb', file_format='PDB')  # Saving residue
        fixer = PDBFixer(filename='resid.pdb')
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingHydrogens(self._settings.pH)  # Adding protons at pH value
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
        mol_resid = atomgroup_to_mol(resid)
        i: int = 0
        for n in mol_resid.GetAtoms():
            i += n.GetNumRadicalElectrons()

        if i > 0:
            backbone = resid.select_atoms('backbone')
            protonated: mda.Universe = mda.Universe.empty(n_atoms=resid.n_atoms + 1, trajectory=True)
            protonated.add_TopologyAttr('masses', [x for x in resid.masses] + [1])
            protonated.add_TopologyAttr('name', [x for x in resid.names] + ['Hc'])
            protonated.add_TopologyAttr('types', guess_types(protonated.atoms.names))
            protonated.add_TopologyAttr('elements', [guess_atom_element(atom) for atom in protonated.atoms.names])
            new_pos = resid.positions
            h_pos = self._get_new_pos(backbone, length)
            protonated.atoms.positions = np.row_stack((new_pos, h_pos))
            return protonated
        else:
            return resid
