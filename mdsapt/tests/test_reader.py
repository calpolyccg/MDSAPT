import pytest

import os

from ..reader import InputReader, InputError


class TestInputReader(object):
    reference = {
        'topology_path': 'mdsapt/tests/testing_resources/testtop.psf',
        'trajectory_path': ['mdsapt/tests/testing_resources/testtraj.dcd'],
        'selection_resid_num': [11, 199],
        'int_pairs': [[11, 199]],
        'trajectory_settings': {
            'start': 1,
            'stop': 2,
            'step': 1
        },
        'system_settings': {
            'ncpus': 4,
            'memory': '12GB',
            'time': '24:00:00'
        }
    }

    def test_reader(self):
        settings = InputReader(os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        assert settings.top_path == self.reference['topology_path']
        assert settings.trj_path == self.reference['trajectory_path']
        assert settings.ag_sel == self.reference['selection_resid_num']
        assert settings.ag_pair == self.reference['int_pairs']
        assert settings.trj_settings == self.reference['trajectory_settings']
        assert settings.sys_settings == self.reference['system_settings']

    @pytest.mark.parametrize('key_pair', [['topology_path', 'error'], ['trajectory_path', ['error']],
                                          ['selection_resid_num', [0]], ['int_pairs', [1]],
                                          ['trajectory_settings', []], ['system_settings', []]])
    def test_check_runinputs(self, key_pair):
        settings = InputReader(os.path.join('mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        alt_set = self.reference
        alt_set[key_pair[0]] = key_pair[1]
        with pytest.raises(InputError):
            settings._check_inputs(alt_set)
