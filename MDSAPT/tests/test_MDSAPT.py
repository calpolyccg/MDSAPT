"""
Unit and regression test for the MDSAPT package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import MDSAPT


def test_MDSAPT_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "MDSAPT" in sys.modules
