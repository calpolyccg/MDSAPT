import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest


import os

import MDAnalysis as mda

from ..repair import rebuild_resid
from ..config import Config, load_from_yaml_file


class TestOptimizer(object):

    def test_prepare_resids(self):
        settings: Config = load_from_yaml_file(
            os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        settings.analysis.pairs = [(11, 214)]
        Unv: mda.Universe = settings.analysis.create_universe()

        r11: mda.AtomGroup = Unv.select_atoms('resid 11')

        r11_fixed: mda.AtomGroup = rebuild_resid(11, r11)
        assert_array_almost_equal(r11.select_atoms(f'name CA').positions, r11_fixed.select_atoms(f'name CA').positions, decimal=3)
        assert_array_almost_equal(r11.select_atoms(f'name N').positions, r11_fixed.select_atoms(f'name N').positions, decimal=3)
        assert_array_almost_equal(r11.select_atoms(f'name O').positions, r11_fixed.select_atoms(f'name O').positions, decimal=3)
        assert_array_almost_equal(r11.select_atoms(f'name C').positions, r11_fixed.select_atoms(f'name C').positions, decimal=3)

    def test_prepare_end(self):
        settings: Config = load_from_yaml_file(
            os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        Unv: mda.Universe = settings.analysis.create_universe()
        settings.analysis.pairs = [(1, 214)]

        r215 = Unv.select_atoms('resid 214')
        r215_fixed = rebuild_resid(214, r215)
        assert len(r215_fixed.select_atoms('name Hc')) == 0

    def test_non_amino(self):
        settings: Config = load_from_yaml_file(
            os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        settings.analysis.pairs = [(126, 214)]
        Unv: mda.Universe = settings.analysis.create_universe()

        r126 = Unv.select_atoms('resid 126')
        r126_fixed = rebuild_resid(126, r126)
        assert len(r126_fixed.select_atoms('name Hc')) == 0
