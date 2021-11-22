"""Prepares AtomGroups for SAPT calculations by adding protons and replacing missing atoms"""

from typing import Dict

import MDAnalysis as mda

import psi4

from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

from .reader import InputReader


class Optimizer(object):
    """Prepares residues for SAPT"""

    _resids: Dict[mda.AtomGroup]
    _unv: mda.Universe
    _settings: InputReader

    def __init__(self, settings: InputReader) -> None:
        self._settings = settings
        self._unv = mda.Universe(self._settings.top_path, self._settings.trj_path)
        self._resids = {x: self._unv.select_atoms(f"resid {x}") for x in self._settings.ag_sel}
        self.fix_protons(settings.trj_settings['pH'])

    def fix_protons(self, ph: float) -> None:
        for k in self._resids:
            res = self._resids[k]
            res.write('resid.pdb', file_format='PDB')
            fixer = PDBFixer(filename='resid.pdb')
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingHydrogens(ph)
            PDBFile.writeFile(fixer.topology, fixer.positions, open('resid_fixed.pdb', 'w'))
            res_fixed = mda.Universe('resid_fixed.pdb')
            resid = mda.Universe.select_atoms('resid *')
            resid.guess_bonds()
            self._resids[k] = resid

    def fix_peptide_carbon(self) -> None:
        res: mda.AtomGroup
        for k in self._resids:
            res = self._resids[k]


    def __getitem__(self, item: int) -> mda.AtomGroup:
        return self._resids[item]
