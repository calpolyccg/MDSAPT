"""Prepares AtomGroups for SAPT calculations by adding protons and replacing missing atoms"""

from typing import Dict

import MDAnalysis as mda

import psi4

from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

import numpy as np

from .reader import InputReader


class Optimizer(object):
    """Prepares residues for SAPT"""

    _resids: Dict[int: mda.AtomGroup]
    _unv: mda.Universe
    _settings: InputReader
    _bond_lenghts: Dict[int: np.ndarray]
    _opt_set: Dict[str, str]

    def __init__(self, settings: InputReader) -> None:
        self._settings = settings
        self._unv = mda.Universe(self._settings.top_path, self._settings.trj_path)
        self._resids = {x: self._unv.select_atoms(f"resid {x}") for x in self._settings.ag_sel}
        self._bond_lenghts = {}
        self._prepare_resids()

    def _prepare_resids(self) -> None:
        for k in self._resids:
            step0: mda.AtomGroup = self._resids[k]
            step1: mda.AtomGroup = self._fix_amino(step0)
            step2: mda.Universe = self._protonate_backbone(step1)
            lenght: float = self._opt_geometry(step2)
            self._bond_lenghts[k] = lenght

    def _fix_amino(self, resid: mda.AtomGroup, ph: float) -> mda.AtomGroup:
        resid.write('resid.pdb', file_format='PDB')  # Saving residue
        fixer = PDBFixer(filename='resid.pdb')
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingHydrogens(ph)  # Adding protons at pH value
        PDBFile.writeFile(fixer.topology, fixer.positions, open('resid_fixed.pdb', 'w'))

        res_fixed = mda.Universe('resid_fixed.pdb')
        resid: mda.AtomGroup = res_fixed.select_atoms("resname *")
        resid.guess_bonds()
        return resid

    def _protonate_backbone(self, resid: mda.AtomGroup, lenght: float = 1.128) -> mda.Universe:
        backbone = resid.select_atoms('backbone')
        carbon = backbone.select_atoms('name C')
        c_pos = carbon.positions
        h_backbone = mda.Universe.empty(n_atoms=backbone + 1, trajectory=True)
        h_backbone.add_TopologyAttr('masses', [x for x in backbone.masses] + [1])
        h_backbone.add_TopologyAttr('name', [x for x in backbone.names] + ['H'])
        back_bone_pos = backbone.positions
        h_backbone.positions = np.row_stack(back_bone_pos, np.array((c_pos[0] - lenght*np.cos(np.pi / 6),
                                                                     c_pos[1] - lenght*np.sin(np.pi / 6),
                                                                     c_pos[2])))
        return h_backbone

    def _opt_geometry(self, system: mda.Universe, basis: str = 'scf/cc-pvdz') -> float:
        tmp_settings: Dict[str, str] = self._opt_set.copy()
        coords: str = ''
        for n in range(len(system.atoms)):
            atom: mda.Atom = system.atoms[n]
            coords += f'\n{atom.name[0]} {atom.position[0]} {atom.position[1]} {atom.position[2]}'
            if not atom.name == 'H':
                tmp_settings['freeze_list'] += f'\n{n + 1} xyz'

        mol = psi4.geometry(coords)
        psi4.set_memory(self._settings.sys_settings['memory'])
        psi4.optimize(basis, molecule=mol)


    def __getitem__(self, item: int) -> mda.AtomGroup:
        return self._resids[item]
