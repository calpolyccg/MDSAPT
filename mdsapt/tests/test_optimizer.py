import numpy as np
import pytest


import os

import MDAnalysis as mda

from ..optimizer import Optimizer
from ..reader import InputReader


class TestOptimizer(object):

    def test_prepare_resids(self):
        settings = InputReader(os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        opt: Optimizer = Optimizer(settings)
        assert len(opt._bond_lengths) == len(settings.ag_sel)
