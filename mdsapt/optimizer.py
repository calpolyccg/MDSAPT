from .reader import InputReader

from typing import Dict

import MDAnalysis as mda

import psi4

from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile


class Optimizer(object):
    """Prepares residues for SAPT"""

    _resids: Dict[mda.AtomGroup]
    _unv: mda.Universe
    _settings: InputReader

    def __init__(self, settings: InputReader) -> None:
        self._settings = settings
        self._unv = mda.Universe(self._settings.top_path, self._settings.trj_path)
        self._resids = {x: self._unv.select_atoms(f"resid {x}") for x in self._settings.ag_sel}


    @staticmethod
    def get_psi_mol(molecule: mda.AtomGroup, ph: float = 7.0) -> str:
        # Based on instructions in https://linuxtut.com/en/30aa73dd6bb949d4858b/

        molecule.write('resid.pdb', file_format='PDB')
        fixer = PDBFixer(filename='resid.pdb')
        fixer.findMissingResidues()

        fixer.addMissingHydrogens(ph)
        PDBFile.writeFile(fixer.topology, fixer.positions, open('resid_fixed.pdb', 'w'))

        mol_fixed = mda.Universe('resid_fixed.pdf')
        mol = mol_fixed.select_atoms('resid 1')

        conf = mol.GetConformer()
        mol_coords = ''
        for atom in mol.GetAtoms():
            mol_coords += (f'\n {atom.GetSymbol()} '
                           f'{conf.GetAtomPosition(atom.GetIdx()).x} '
                           f'{conf.GetAtomPosition(atom.GetIdx()).y} '
                           f'{conf.GetAtomPosition(atom.GetIdx()).z}')

        psi_mol: psi4.core.Molecule = psi4.geometry(mol_coords)
        coords = psi_mol.save_string_xyz()
        coord_header = f'{mol.getFormalCharge(), psi_mol.get_fragment_multiplicities()}'
        coords = coords.split('\n')
        coords = coord_header.join(['\n' + coord for coord in coords[1:]])
        return psi_mol.save_string_xyz()

