"""Provide the primary functions."""

from typing import Dict

import pandas as pd

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.topology.guessers import guess_types

import psi4

from rdkit import Chem

from pdbfixer import PDBFixer

from openmm.app import PDBFile

from .reader import InputReader


class TrajectorySAPT(AnalysisBase):
    """Calculates the <SAPT>`` for individual frames of a trajectory.
    """
    def __init__(self, config: InputReader, **universe_kwargs):
        self._unv: mda.Universe = mda.Universe(config.top_path, config.trj_path, **universe_kwargs)
        elements = guess_types(self._unv.atoms.names)
        self._unv.add_TopologyAttr('elements', elements)
        self._sel: Dict[mda.AtomGroup] = {x: self._unv.select_atoms(f'resid {x}') for x in config.ag_sel}
        self._sel_pairs = config.ag_pair
        self._mem = config.sys_settings['memory']
        self._cfg = config
        super(TrajectorySAPT, self).__init__(self._unv.trajectory)

    def _prepare(self):
        self._col = ['residues', 'time', 'energy']
        self.results = pd.DataFrame(columns=self._col)
        self._res_dict = {x: [] for x in self._col}

    @staticmethod
    def get_psi_mol(molecule: mda.AtomGroup, ph: float=7.0):
        # Based on instructions in https://linuxtut.com/en/30aa73dd6bb949d4858b/

        molecule.write('resid.pdb', file_format='PDB')
        fixer = PDBFixer(filename='resid.pbd')
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)
        PDBFile.writeFile(fixer.topology, fixer.positions, open('resid_fixed.pdb', 'w'))

        mol = Chem.MolFromMolFile('resid_fixed.pdb')
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

    def _single_frame(self):
        xyz_dict = {k: self.get_psi_mol(self._sel[k]) for k in self._sel.keys()}
        with open('test.xyz', 'w+') as r:
            r.write(xyz_dict[2])
        for pair in self._sel_pairs:
            coords = xyz_dict[pair[0]] + '--\n' + xyz_dict[pair[1]] + 'units angstrom'
            dimer = psi4.geometry(coords)
            psi4.set_options({'scf_type': 'df',
                              'freeze_core': 'true'})
            psi4.set_memory(self._mem)

            sapt = psi4.energy('sapt0/jun-cc-pvdz', molecule=dimer)
            result = [f'{pair[0]}-{pair[1]}', self._ts.time, sapt]
            for r in range(len(result)):
                self._res_dict[self._col[r]].append(result[r])

    def _conclude(self):
        for k in self._col:
            self.results[k] = self._res_dict[k]
