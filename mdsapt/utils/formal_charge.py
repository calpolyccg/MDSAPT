r"""
:mod:`mdsapt.utils.formal_charge` -- Utilities for calculating charge-related properties of an atom
===================================================================================================

.. autofunction:: get_fc

.. autofunction:: get_lone_electrons
"""
from typing import Optional, NamedTuple


_ELEMENT_TO_VALENCE = {
    'H': 1,
    'C': 4,
    'N': 5,
    'O': 6,
    'P': 5,
    'S': 6,

    # Halogens
    'F': 7,
    'Cl': 7,
    'Br': 7,
    'I': 7,
}


class ElectronInfo(NamedTuple):
    """Various electron-related properties of an atom."""
    element: str
    valence: int
    bonds: int
    fc: int
    lone: int
    radical: int


def calculate_electron_info(element: str, bonds: int, fc: Optional[int] = None) -> ElectronInfo:
    """
    Calculates various electron-related properties of an atom.

    :param element: the element of the atom
    :param bonds: number of bonds attached
    :param fc:
        Formal charge of the atom, if known. If not provided, it will be automatically calculated.
        Our formal charge calculation is optimized for biochemistry simulations, so it might not produce the correct
        values for sulfur or phosphorous in more complicated molecules.
    """

    try:
        valence = _ELEMENT_TO_VALENCE[element]
    except KeyError:
        raise ValueError(f"Unsupported element {element}, please manually specify the charge.")

    if fc is None:
        fc = _get_fc(element, valence, bonds)

    # Note that we calculate lone back from FC for two reasons:
    # - charge is provided by the user
    # - we pin the charges of certain element/bond inputs
    lone = valence - fc - bonds

    # Calculation of radical electrons. This formula appears to work in most cases.
    if bonds == 1 and lone <= 4: 
        radical = [0, 1, 0, 3, 2][lone]  # lookup table
    else:
        radical = lone % 2

    return ElectronInfo(
        element=element,
        valence=valence,
        bonds=bonds,
        fc=fc,
        lone=lone,
        radical=radical,
    )


def calculate_spin_multiplicity(total_radicals: int) -> int:
    """
    Calculates the spin multiplicity of a molecule, given the total number of radical electrons.
    
    :param total_radicals: The number of radical electrons that exist across the molecule.
    """
    total_spin: int = total_radicals // 2
    spin_mult: int = total_spin + 1
    return spin_mult


def _get_fc(element: str, valence: int, bonds: int) -> int:
    """
    Calculates the formal charge of an atom given its number of bonds.

    This is optimized for biochemistry simulations, so it might not produce the correct
    values for sulfur or phosphorous in more complicated molecules.

    :param element: the element of the atom
    :param valence: number of valence electrons
    :param bonds: number of bonds attached to the atom
    """
    if element == 'H':
        return 0
    if element == 'P':
        if bonds == 5:
            return 0
        if bonds == 6:
            return -1

    unpaired = 8 - valence

    if bonds <= unpaired:
        lone = valence - bonds

        # Round lone electrons up to next even number
        if lone % 2 != 0:
            lone += 1
    else:
        # Every bond above unpaired donates 2 lone electrons
        lone = 8 - 2 * bonds

    return valence - lone - bonds
