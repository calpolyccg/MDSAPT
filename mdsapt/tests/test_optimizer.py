import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest


import os

import MDAnalysis as mda

from ..optimizer import Optimizer
from ..reader import InputReader


class TestOptimizer(object):

    def test_prepare_resids(self):
        settings = InputReader(os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        U = mda.Universe(settings.top_path, settings.trj_path)
        Opt: Optimizer = Optimizer(settings)

        r11: mda.AtomGroup = U.select_atoms('resid 11')

        r11_fixed: mda.AtomGroup = Opt.rebuild_resid(11, r11)
        assert_array_almost_equal(r11.select_atoms(f'name CA').positions, r11_fixed.select_atoms(f'name CA').positions, decimal=3)
        assert_array_almost_equal(r11.select_atoms(f'name N').positions, r11_fixed.select_atoms(f'name N').positions, decimal=3)
        assert_array_almost_equal(r11.select_atoms(f'name O').positions, r11_fixed.select_atoms(f'name O').positions, decimal=3)
        assert_array_almost_equal(r11.select_atoms(f'name C').positions, r11_fixed.select_atoms(f'name C').positions, decimal=3)

    def test_prepare_end(self):
        settings = InputReader(os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        U = mda.Universe(settings.top_path, settings.trj_path)
        Opt: Optimizer = Optimizer(settings)
        Opt._resids[214] = U.select_atoms('resid 214')

        r215 = U.select_atoms('resid 214')
        r215_fixed = Opt.rebuild_resid(214, r215)
        assert len(r215_fixed.select_atoms('name Hc')) == 0


