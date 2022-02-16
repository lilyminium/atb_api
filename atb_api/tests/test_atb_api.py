"""
Unit and regression test for the atb_api package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import atb_api


def test_atb_api_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "atb_api" in sys.modules
