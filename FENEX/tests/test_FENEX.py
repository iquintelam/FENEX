"""
Unit and regression test for the FENEX package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import FENEX


def test_FENEX_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "FENEX" in sys.modules
