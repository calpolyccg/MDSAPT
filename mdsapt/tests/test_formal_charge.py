from typing import Optional
import pytest
from ..utils.formal_charge import _get_fc, calculate_electron_info, ElectronInfo


@pytest.mark.parametrize('element, bonds, fc, fields', [
    # fmt: off

    ('H',  1, None, dict(valence=1, fc= 0, lone=0, radical=0)),

    # Typical C bonding
    ('C',  2, None, dict(valence=4, fc= 0, lone=2, radical=0)),  # Carbene
    ('C',  3, None, dict(valence=4, fc=-1, lone=2, radical=0)),  # C in Carbon monoxide (CO)
    ('C',  4, None, dict(valence=4, fc= 0, lone=0, radical=0)),  # Normal C
    ('C',  1,    0, dict(valence=4, fc= 0, lone=3, radical=3)),  # Carbyne
    ('C',  3,   -1, dict(valence=4, fc=-1, lone=2, radical=0)),  # Carboanion
    ('C',  3,    1, dict(valence=4, fc= 1, lone=0, radical=0)),  # Carbocation

    ('N',  1, None, dict(valence=5, fc= 0, lone=4, radical=2)),  # Nitrene
    ('N',  2, None, dict(valence=5, fc=-1, lone=4, radical=0)),  # N in Nitric oxide
    ('N',  3, None, dict(valence=5, fc= 0, lone=2, radical=0)),  # Normal N
    ('N',  4, None, dict(valence=5, fc= 1, lone=0, radical=0)),  # N in NH4+

    ('O',  1, None, dict(valence=6, fc=-1, lone=6, radical=0)),  # O in Hydroxide ion (OH-)
    ('O',  2, None, dict(valence=6, fc= 0, lone=4, radical=0)),  # Normal O
    ('O',  3, None, dict(valence=6, fc= 1, lone=2, radical=0)),  # O in Hydronium (H3O+)
    ('O',  1,    0, dict(valence=6, fc= 0, lone=5, radical=1)),  # O in Hydroxyl radical (OH)

    ('P',  5, None, dict(valence=5, fc= 0, lone=0, radical=0)),  # P in Phosphate
    ('P',  6, None, dict(valence=5, fc=-1, lone=0, radical=0)),  # P in PF6-

    # Halogens
    ('F',  1, None, dict(valence=7, fc= 0, lone=6, radical=0)),
    ('Cl', 1, None, dict(valence=7, fc= 0, lone=6, radical=0)),
    ('Br', 1, None, dict(valence=7, fc= 0, lone=6, radical=0)),
    ('I',  1, None, dict(valence=7, fc= 0, lone=6, radical=0)),

    # fmt: on
])
def test_calculate_electron_info(element: str, bonds: int, fc: Optional[int], fields: dict):
    actual = calculate_electron_info(element, bonds, fc)

    expected = ElectronInfo(element=element, bonds=bonds, **fields)
    assert actual == expected
