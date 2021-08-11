"""
Unit and regression test for the mdsapt package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import mdsapt


def test_MDSAPT_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "mdsapt" in sys.modules
