"""Tests for linear constraint assembly."""
import numpy as np
import pytest
from now.config import NOW_config
from now.constraints.linear import get_linear_constraint_matrices
from now.utils import build_first_derivative_matrix, build_second_derivative_matrix
from .helpers import load_fixture, skip_without_fixture


class TestDerivativeMatrices:
    """Test finite difference matrices against MATLAB."""

    def test_first_derivative_shape(self):
        A1 = build_first_derivative_matrix(20)
        assert A1.shape == (19, 20)

    def test_second_derivative_shape(self):
        A2 = build_second_derivative_matrix(20)
        assert A2.shape == (20, 20)

    def test_first_derivative_structure(self):
        A1 = build_first_derivative_matrix(5)
        expected = np.array([
            [-1, 1, 0, 0, 0],
            [0, -1, 1, 0, 0],
            [0, 0, -1, 1, 0],
            [0, 0, 0, -1, 1],
        ])
        np.testing.assert_array_equal(A1, expected)

    def test_second_derivative_structure(self):
        A2 = build_second_derivative_matrix(4)
        expected = np.array([
            [-2, 1, 0, 0],
            [1, -2, 1, 0],
            [0, 1, -2, 1],
            [0, 0, 1, -2],
        ])
        np.testing.assert_array_equal(A2, expected)

    @skip_without_fixture('constraints_ste_default.mat')
    def test_against_matlab(self):
        f = load_fixture('constraints_ste_default.mat')
        N = int(f['N'])
        A1 = build_first_derivative_matrix(N)
        A2 = build_second_derivative_matrix(N)
        np.testing.assert_array_almost_equal(A1, f['firstDerivativeMatrix'])
        np.testing.assert_array_almost_equal(A2, f['secondDerivativeMatrix'])


class TestLinearConstraints:
    """Test assembled linear constraint matrices."""

    def test_inequality_shape_default(self):
        c = NOW_config(N=20)
        A_ineq, b_ineq, _, _ = get_linear_constraint_matrices(c)
        # No max-norm → only slew rate: 2 * (3*N) rows
        assert A_ineq.shape == (2 * 3 * 20, 3 * 20 + 1)
        assert len(b_ineq) == A_ineq.shape[0]

    def test_inequality_shape_maxnorm(self):
        c = NOW_config(N=20, useMaxNorm=True)
        A_ineq, b_ineq, _, _ = get_linear_constraint_matrices(c)
        # max-norm: 2 * 3*(N-1) + 2 * 3*N rows
        expected_rows = 2 * 3 * 19 + 2 * 3 * 20
        assert A_ineq.shape[0] == expected_rows

    def test_equality_echo_condition(self):
        c = NOW_config(N=20)
        _, _, A_eq, b_eq = get_linear_constraint_matrices(c)
        # b_eq should be all zeros
        np.testing.assert_array_equal(b_eq, 0)
        # A_eq should enforce q(0)=0 and q(N-1)=0 for each of 3 axes
        # Plus zero-gradient constraints

    def test_equality_no_pause(self):
        c = NOW_config(N=20, durationZeroGradientRequested=0)
        _, _, A_eq, b_eq = get_linear_constraint_matrices(c)
        # Only echo condition: 2 rows per axis = 6 rows
        assert A_eq.shape[0] == 6

    @skip_without_fixture('constraints_ste_default.mat')
    def test_inequality_against_matlab(self):
        f = load_fixture('constraints_ste_default.mat')
        c = NOW_config()
        A_ineq, b_ineq, _, _ = get_linear_constraint_matrices(c)
        np.testing.assert_array_almost_equal(A_ineq, f['A_ineq'], decimal=10)
        np.testing.assert_array_almost_equal(b_ineq, f['b_ineq'], decimal=10)

    @skip_without_fixture('constraints_ste_default.mat')
    def test_equality_against_matlab(self):
        f = load_fixture('constraints_ste_default.mat')
        c = NOW_config()
        _, _, A_eq, b_eq = get_linear_constraint_matrices(c)
        np.testing.assert_array_almost_equal(A_eq, f['A_eq'], decimal=10)
        np.testing.assert_array_almost_equal(b_eq, f['b_eq'], decimal=10)

    @skip_without_fixture('constraints_ste_maxnorm.mat')
    def test_maxnorm_against_matlab(self):
        f = load_fixture('constraints_ste_maxnorm.mat')
        c = NOW_config(useMaxNorm=True)
        A_ineq, b_ineq, _, _ = get_linear_constraint_matrices(c)
        np.testing.assert_array_almost_equal(A_ineq, f['A_ineq'], decimal=10)
        np.testing.assert_array_almost_equal(b_ineq, f['b_ineq'], decimal=10)

    @skip_without_fixture('constraints_ste_motion.mat')
    def test_motion_comp_against_matlab(self):
        f = load_fixture('constraints_ste_motion.mat')
        mc = {'order': [1, 2], 'maxMagnitude': [0, 1e-4]}
        c = NOW_config(motionCompensation=mc)
        _, _, A_eq, b_eq = get_linear_constraint_matrices(c)
        np.testing.assert_array_almost_equal(A_eq, f['A_eq'], decimal=10)

    @skip_without_fixture('constraints_ste_bgcomp.mat')
    def test_bgcomp_against_matlab(self):
        f = load_fixture('constraints_ste_bgcomp.mat')
        mc = {'order': [1], 'maxMagnitude': [0]}
        c = NOW_config(doBackgroundCompensation=1, motionCompensation=mc)
        _, _, A_eq, b_eq = get_linear_constraint_matrices(c)
        np.testing.assert_array_almost_equal(A_eq, f['A_eq'], decimal=10)

    @skip_without_fixture('constraints_ste_symmetric.mat')
    def test_symmetry_against_matlab(self):
        f = load_fixture('constraints_ste_symmetric.mat')
        c = NOW_config(enforceSymmetry=True,
                       durationFirstPartRequested=25,
                       durationSecondPartRequested=25,
                       durationZeroGradientRequested=8)
        _, _, A_eq, b_eq = get_linear_constraint_matrices(c)
        np.testing.assert_array_almost_equal(A_eq, f['A_eq'], decimal=10)
