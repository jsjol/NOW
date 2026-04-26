"""Shared test helpers for MATLAB fixture loading."""
import numpy as np
import pytest
from pathlib import Path

FIXTURE_DIR = Path(__file__).parent / 'fixtures'


def load_fixture(name):
    """Load a .mat fixture file. Returns None if file doesn't exist."""
    filepath = FIXTURE_DIR / name
    if not filepath.exists():
        return None
    from scipy.io import loadmat
    return loadmat(str(filepath), squeeze_me=True)


def matlab_fixture_available(name):
    """Check if a MATLAB fixture file exists."""
    return (FIXTURE_DIR / name).exists()


def skip_without_fixture(name):
    """Decorator to skip test if MATLAB fixture is not available."""
    return pytest.mark.skipif(
        not matlab_fixture_available(name),
        reason=f"MATLAB fixture {name} not generated yet. Run generate_test_fixtures.m in MATLAB.")
