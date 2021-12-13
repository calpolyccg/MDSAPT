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
            step0: mda.AtomGroup = self._resids[k]  # Get resid for optimization
            step1: mda.AtomGroup = self._fix_amino(step0)  # Fix amino group
            step2: mda.Universe = self._protonate_backbone(step1)  # Add proton to backbone
            lenght: float = self._opt_geometry(step2)  # Get optimized new C-H bond
            self._bond_lenghts[k] = lenght  # Hash new bond length

    def _rebuild_resid(self, key: int) -> None:
        pass

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
    def _protonate_backbone(resid: mda.AtomGroup, lenght: float = 1.128) -> mda.Universe:
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

    def __getitem__(self, item: int) -> mda.AtomGroup:
        return self._resids[item]

    def set_opt_setting(self, opt_settings: Dict[str, str]) -> None:
        self._opt_set = opt_settings
