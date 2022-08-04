ELEMENT_TO_VALENCE = {
    'C': 4,
    'N': 5,
    'O': 6,
    'F': 7,
    'P': 5,
    'S': 6,
    'Cl': 7,
    'I': 7,
}

def get_fc(element: str, bonds: int) -> int:
    """
    Calculates the formal charge of an atom given its number of bonds.

    This is optimized for biochemistry simulations, so it might not produce the correct
    values for sulfur in more complicated molecules.
    """
    if element == 'H':
        return 0
    if element == 'S':
        return get_fc('O', bonds)
    if element == 'P':
        if bonds == 5:
            return 0
        if bonds == 6:
            return 1

    try:
        valence = ELEMENT_TO_VALENCE[element]
    except KeyError:
        raise ValueError(f"Unsupported element {element}, please manually specify the charge.")

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
