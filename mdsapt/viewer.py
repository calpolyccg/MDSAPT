

import nglview as nv

from .sapt import TrajectorySAPT
from .reader import InputReader


class Viewer:
    def __init__(self, trj_sapt: TrajectorySAPT):
        self._unv = trj_sapt._unv
        self._sel = trj_sapt._sel


