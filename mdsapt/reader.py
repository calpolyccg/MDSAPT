import os
from typing import List, Optional

import yaml

import MDAnalysis as mda

import logging
logger = logging.getLogger('mdsapt.config')


class InputError(Exception):
    pass


class CoordsReader(object):
    """"""
    def __init__(self, universe: Optional[mda.Universe] = None, selections: Optional[List[mda.AtomGroup]] = None, **kwargs) -> None:
        if not universe is None:
            self._unv = universe
        if not selections is None:
            self._sel = selections

    def write_xyz(self, selection: mda.AtomGroup, time: int, pathname: str) -> None:
        pathname += '.xyz'
        with mda.Writer(pathname, selection[time].n_atoms) as coords:
            coords.write(selection)

    def get_int_pair(self, sel1: mda.AtomGroup, sel2: mda.AtomGroup, time_step: int) -> List[str]:
        self.write_xyz(sel1, time_step, 'res1')
        self.write_xyz(sel2, time_step, 'res2')

        resid_list: str = ''

        resid_list += f'{sel1.total_charge()} 1\n'


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
    def remove_xyz(pathname: str) -> None:
        try:
            os.remove(pathname)
        except IOError:
            pass

    def save_sapt_in(self, path: str, coords0: List[str], coords1: List[str], molecule_name: str,
                     memory: str, char0: str, char1: str) -> None:

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


class InputReader(object):
    """Reader for yaml inputs"""
    def __init__(self, path):
        self.load_input(path)
        self.top_path: Optional[str]
        self.trj_path: Optional[List[str]]
        self.ag_sel: Optional[List[int]]
        self.ag_names: Optional[List[str]]
        self.ag_pair: Optional[List[List[str]]]
        self.trj_settings: Optional[dict]
        self.sys_settings: Optional[dict]

    def load_input(self, path: str) -> None:
        """Reads yaml file"""
        try:
            in_cfg = yaml.safe_load(path)
            self.check_inputs(in_cfg)
            self.save_params()
        except IOError:
            logger.warning('error loading file')

    def save_params(self, yaml_dict: dict) -> None:
        self.top_path = yaml_dict['topology_path']
        self.trj_path = yaml_dict['trajectory_path']
        self.ag_sel = yaml_dict['selection_resid_nums']
        self.ag_names = yaml_dict['selection_names']
        self.ag_pair = yaml_dict['int_pairs']
        self.trj_settings = yaml_dict['trajectory_settings']
        self.sys_settings = yaml_dict['system_settings']

    @staticmethod
    def check_inputs(yaml_dict: dict) -> None:
        try:
            top_path = yaml_dict['topology_path']
            trj_path = yaml_dict['trajectory_path']
            ag_sel = yaml_dict['selection_resid_nums']
            ag_names = yaml_dict['selection_names']
            ag_pair = yaml_dict['int_pairs']
            trj_settings = yaml_dict['trajectory_settings']
            sys_settings = yaml_dict['system_settings']
        except KeyError:
            logger.fatal('Invalid YAML file')

        try:
            unv = mda.Universe(top_path, trj_path)
        except mda.exceptions.NoDataError:
            logger.fatal('MD file error')

        # Testing names and selections
        if len(ag_sel) > len(ag_names):
            raise InputError('Not all selections are named')
        elif len(ag_sel) < len(ag_names):
            raise InputError('Too many selection names for number of selections')

        for sel in ag_sel:
            try:
                ag = unv.select_atoms(sel)
            except mda.SelectionError:
                raise InputError('Error in selection: {}'.format(sel))

        try:
            start = trj_settings['start']
            step = trj_settings['step']
            stop = trj_settings['stop']
            cpu = sys_settings['ncpu']
            mem = sys_settings['memory']
            time = sys_settings['time']

            for pair in ag_pair:
                if len(pair) != 4:
                    raise InputError('Pairs must be a python list of string with 4 items')
                found0 = False
                found1 = False
                for name in ag_names:
                    if pair[0] == name:
                        found0 = True
                    if pair[1] == name:
                        found1 = True
                if found0 is False:
                    raise InputError(f'{pair[0]} in {pair} group_pair_selections is not in defined in atom_group_names')
                if found1 is False:
                    raise InputError(f'{pair[1]} in {pair} group_pair_selections is not in defined in atom_group_names')

                if start >= stop:
                    raise InputError('Start is greater than or equal to stop')
                if step >= stop:
                    raise InputError('Step is greater than or equal to stop')
                if step == 0:
                    raise InputError('Step cannot be 0')

                if len(unv.trajectory) < stop:
                    raise InputError(f'Stop exceeds length of trajectory, trajectory is '
                                     f'{len(unv.trajectory)} frames')

        except KeyError or InputError as err:
            logger.fatal(f'{err} Error in trajectory settings')

        logger.info('Input Parameters Accepted')
