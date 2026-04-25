"""Tests for the objective function."""
import numpy as np
from now.optimize import objective


class TestObjective:
    def test_returns_tuple(self):
        x = np.array([1.0, 2.0, 3.0])
        f, g = objective(x)
        assert isinstance(f, float)
        assert g.shape == x.shape

    def test_value(self):
        x = np.array([1.0, 2.0, 5.0])
        f, g = objective(x)
        assert f == -5.0

    def test_gradient(self):
        x = np.array([1.0, 2.0, 3.0])
        f, g = objective(x)
        expected = np.array([0.0, 0.0, -1.0])
        np.testing.assert_array_equal(g, expected)

    def test_gradient_finite_diff(self):
        x = np.random.randn(31)
        f, g = objective(x)
        eps = 1e-7
        g_fd = np.zeros_like(x)
        for i in range(len(x)):
            xp = x.copy(); xp[i] += eps
            xm = x.copy(); xm[i] -= eps
            g_fd[i] = (objective(xp)[0] - objective(xm)[0]) / (2 * eps)
        np.testing.assert_array_almost_equal(g, g_fd)
