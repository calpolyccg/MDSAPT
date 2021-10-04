import pytest

import os

import MDAnalysis as mda
from MDAnalysis.converters.RDKit import atomgroup_to_mol
from MDAnalysis.topology.guessers import guess_types

from rdkit.Chem.rdmolops import AddHs
import rdkit.Chem as chem

from ..reader import InputReader
from ..sapt import TrajectorySAPT


class TestSAPT(object):

    met1 = """
1 1
 C    0.640221269280   -2.365527384856   -0.220701583902
 O   -0.375247328132   -2.845901721099   -0.601017364542
 C    1.138295800835   -1.101050608733   -0.829104789774
 N    1.209717424065   -0.993726008513   -2.291155227701
 H    2.146285684258   -0.934034579375   -0.526638397257
 C    0.213311822564    0.075404889009   -0.436624893229
 H    1.838792474419   -1.654587500670   -2.680715927164
 H    1.565622003228   -0.052988284209   -2.570485481302
 H    0.305630357415   -1.185909503081   -2.809805282633
 H   -0.887536375373   -0.224038355925   -0.813826927225
 H    0.568287522942    1.024360424897   -0.987870582621
 C    0.030244500786    0.439716107270    1.041762939413
 H    0.922888429314    0.828102833650    1.460975280722
 H   -0.249183027595   -0.515948527434    1.556871047934
 S   -1.276656477302    1.610320813081    1.244025340994
 C   -0.082423536628    2.997824437043    0.992187610586
 H    0.946120889336    2.799901730439    0.493072143515
 H    0.038891465813    3.573343045137    1.902494064291
 H   -0.643233625740    3.677866703889    0.364156833609
"""

    def setup(self):
        self.settings = InputReader(os.path.join(os.getcwd(), 'mdsapt', 'tests', 'testing_resources', 'test_input.yaml'))
        self.Unv = mda.Universe(self.settings.top_path, self.settings.trj_path)
        elements = guess_types(self.Unv.atoms.names)
        self.Unv.add_TopologyAttr('elements', elements)

    def test_psi4_in(self):
        res1 = self.Unv.select_atoms('resid 1')
        mol1 = res1.convert_to('RDKIT')
        mol1 = chem.AddHs(mol1)
        mol1_in = TrajectorySAPT.get_psi_mol(mol1)
        assert mol1_in.split('\n') == self.met1.split('\n')

    def test_sapt_run(self):
        SAPT1 = TrajectorySAPT(self.settings).run(stop=2)
