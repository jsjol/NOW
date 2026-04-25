import pytest
import numpy as np
import sys
from pathlib import Path

# Add the NOW directory to path so `now` package can be imported
sys.path.insert(0, str(Path(__file__).parent.parent))


@pytest.fixture
def default_config():
    """Default STE configuration matching MATLAB defaults."""
    from now.config import NOW_config
    return NOW_config()


@pytest.fixture
def small_config():
    """Small N config for quick tests."""
    from now.config import NOW_config
    return NOW_config(N=20)


@pytest.fixture
def deterministic_x0():
    """Generate deterministic x0."""
    def _make(config):
        rng = np.random.RandomState(42)
        x0 = rng.randn(3 * config.N + 1)
        x0[:config.N] *= config.targetTensor[0, 0]
        x0[config.N:2 * config.N] *= config.targetTensor[1, 1]
        x0[2 * config.N:3 * config.N] *= config.targetTensor[2, 2]
        return x0
    return _make
