import os
from typing import List, Optional

import MDAnalysis as mda


class CoordsReader(object):
    """"""
    def __init__(self, universe: Optional[mda.Universe] = None, selections: Optional[List[mda.AtomGroup]] = None, **kwargs):
        if not universe is None:
            self._unv = universe
        if not selections is None:
            self._sel = selections

    def write_xyz(self, selection: mda.AtomGroup, time: int, pathname: str):
        pathname += '.xyz'
        with mda.Writer(pathname, selection[time].n_atoms) as coords:
            coords.write(selection)

    @staticmethod
    def read_xyz(xyz_path: str) -> List[str]:
        with open(xyz_path, 'r') as coord_file:
            xyz_data = []
            coord_data = coord_file.readlines()[2:]

            for line in coord_data:
                if '.' in line:
                    xyz_data.append(line)
            return xyz_data

    @staticmethod
    def remove_xyz(pathname: str):
        try:
            os.remove(pathname)
        except IOError:
            pass

    def save_sapt_in(self, path: str, coords0: List[str], coords1: List[str], molecule_name: str,
                     memory: str, char0: str, char1: str):

        coord_data = 'molecule %s {\n' % molecule_name
        coord_data += f'{char0}\n'

        for line0 in coords0:
            items = line0.split()
            line0 = items[0][0]
            for item in items[1:]:
                line0 += ('\t' + item)
            coord_data += (line0 + '\n')

        coord_data += '--\n'
        coord_data += f'{char1}\n'

        for line1 in coords1:
            items = line1.split()
            line1 = items[0][0]
            for item in items[1:]:
                line1 += ('\t' + item)
            coord_data += (line1 + '\n')

        coord_data += '\nunits angstrom\n' \
                      '\n' \
                      '}\n' \
                      '\nset {\n' \
                      'basis jun-cc-pVDZ\n' \
                      'scf_type DF\n' \
                      'freeze_core True\n' \
                      '}\n'

        coord_data += '\n' + 'memory ' + str(memory) + ' GB\n'

        coord_data += "\nenergy('sapt0')\n"

        with open(path, 'w+') as input_file:
            for line in coord_data:
                input_file.write(line)

