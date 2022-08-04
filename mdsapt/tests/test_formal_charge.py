import pytest
from ..utils.formal_charge import get_fc

@pytest.mark.parametrize('element, bonds, expected', [
    ('C', 2, 0),
    ('C', 3, -1),
    ('C', 4, 0),
    ('O', 2, 0),
    ('O', 3, 1),
    ('O', 1, -1),
    ('N', 2, -1),
    ('N', 3, 0),
    ('N', 4, 1),
])
def test_formal_charge_is_accurate(element: str, bonds: int, expected: int):
    assert get_fc(element, bonds) == expected
