import copy

import pytest

import os

from ..config import InputReader, InputError


class TestInputReader(object):
    reference = {
        'topology_path': 'mdsapt/tests/testing_resources/testtop.psf',
        'trajectory_paths': ['mdsapt/tests/testing_resources/testtraj.dcd'],
        'selection_resid_num': [11, 199],
        'int_pairs': [[11, 199]],
        'trajectory_settings': {
            'start': 1,
            'stop': 2,
            'step': 1
        },
        'system_settings': {
            'ncpus': 16,
            'memory': '48GB',
            'time': '24:00:00'
        },
        'opt_settings': {
            'pH': 7
        },
        'sapt_settings': {
            'method': 'sapt0',
            'basis': 'jun-cc-pvdz',
            'settings': {'reference': 'rhf'},
            'save_psi4_output': True
        }
    }

    def test_reader(self):
        settings = InputReader(os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        assert settings.top_path == self.reference['topology_path']
        assert settings.trj_path == self.reference['trajectory_paths']
        assert settings.ag_sel == self.reference['selection_resid_num']
        assert settings.ag_pair == self.reference['int_pairs']
        assert settings.trj_settings == self.reference['trajectory_settings']
        assert settings.sys_settings == self.reference['system_settings']

    @pytest.mark.parametrize('key_pair', [['topology_path', 'error'], ['trajectory_paths', ['error']],
                                          ['selection_resid_num', [0]], ['int_pairs', [[1]]],
                                          ['int_pairs', [[1, 199]]], ['int_pairs', [[11, 1]]]])
    def test_check_runinputs(self, key_pair):
        alt_set = self.reference.copy()
        alt_set[key_pair[0]] = key_pair[1]
        with pytest.raises(InputError):
            InputReader._check_inputs(alt_set)

    def test_fail_load_input(self):
        with pytest.raises(InputError):
            settings = InputReader(os.path.join('mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
            settings.load_input('error')

    def test_missing_yaml_data0(self):
        alt_set = copy.deepcopy(self.reference)
        alt_set.pop('topology_path')
        with pytest.raises(InputError):
            InputReader._check_inputs(alt_set)

    def test_missing_yaml_data1(self):
        alt_set = copy.deepcopy(self.reference)
        alt_set['opt_settings'].pop('pH')
        with pytest.raises(InputError):
            InputReader._check_inputs(alt_set)

    @pytest.mark.parametrize('key_pair', [['start', 10], ['step', 10],
                                          ['step', 0], ['stop', 1000]])
    def test_error_in_trj_settings(self, key_pair):
        alt_set = copy.deepcopy(self.reference)
        alt_set['trajectory_settings'][key_pair[0]] = key_pair[1]
        with pytest.raises(InputError):
            InputReader._check_inputs(alt_set)

    def test_check_resids(self):
        alt_set = copy.deepcopy(self.reference)
        alt_set['selection_resid_num'] = [11, 199, None]
        with pytest.raises(InputError):
            InputReader._check_inputs(alt_set)
