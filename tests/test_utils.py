"""Tests for utility functions."""
import numpy as np
from now.utils import (
    gamma, build_first_derivative_matrix, build_second_derivative_matrix,
    build_integration_weights, gwf_to_q, problem_to_name
)
from now.config import NOW_config


class TestGamma:
    def test_value(self):
        assert abs(gamma() - 42.576e6) < 1e-3


class TestIntegrationWeights:
    def test_shape(self):
        w = build_integration_weights(10)
        assert w.shape == (10, 1)

    def test_trapezoidal(self):
        w = build_integration_weights(5)
        assert w[0, 0] == 0.5 / 4
        assert w[-1, 0] == 0.5 / 4
        assert w[1, 0] == 1.0 / 4

    def test_sum(self):
        w = build_integration_weights(100)
        assert abs(w.sum() - 1.0) < 1e-10


class TestGwfToQ:
    def test_constant_gradient(self):
        gwf = np.ones((10, 3))
        dt = 0.5
        q = gwf_to_q(gwf, dt)
        assert q.shape == (10, 3)
        np.testing.assert_array_almost_equal(q[-1], [5.0, 5.0, 5.0])


class TestProblemToName:
    def test_ste(self):
        c = NOW_config()
        name = problem_to_name(c)
        assert 'STE' in name

    def test_lte(self):
        c = NOW_config(targetTensor=np.diag([1, 0, 0]))
        name = problem_to_name(c)
        assert 'LTE' in name

    def test_pte(self):
        c = NOW_config(targetTensor=np.diag([1, 1, 0]))
        name = problem_to_name(c)
        assert 'PTE' in name
