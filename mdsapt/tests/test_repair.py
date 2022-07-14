from pathlib import Path

from numpy.testing import assert_array_almost_equal

import MDAnalysis as mda

from ..repair import rebuild_resid

resources_dir = Path(__file__).parent / 'testing_resources'

unv: mda.Universe = mda.Universe(f'{resources_dir}/testtop.psf', f'{resources_dir}/testtraj.dcd')


class TestOptimizer(object):

    def test_prepare_resids(self):
        r11: mda.AtomGroup = unv.select_atoms('resid 11')
        r11_fixed: mda.AtomGroup = rebuild_resid(11, r11)
        assert_array_almost_equal(r11.select_atoms(f'name CA').positions, r11_fixed.select_atoms(f'name CA').positions, decimal=3)
        assert_array_almost_equal(r11.select_atoms(f'name N').positions, r11_fixed.select_atoms(f'name N').positions, decimal=3)
        assert_array_almost_equal(r11.select_atoms(f'name O').positions, r11_fixed.select_atoms(f'name O').positions, decimal=3)
        assert_array_almost_equal(r11.select_atoms(f'name C').positions, r11_fixed.select_atoms(f'name C').positions, decimal=3)

    def test_prepare_end(self):
        r215 = unv.select_atoms('resid 214')
        r215_fixed = rebuild_resid(214, r215)
        assert len(r215_fixed.select_atoms('name Hc')) == 0

    def test_non_amino(self):
        r126 = unv.select_atoms('resid 126')
        r126_fixed = rebuild_resid(126, r126)
        assert len(r126_fixed.select_atoms('name Hc')) == 0
