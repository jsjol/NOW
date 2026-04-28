"""Tests for the solver-independent problem representation."""
import numpy as np
from now.config import NOW_config
from now.problem import (
    build_problem, OptimizationProblem, LinearConstraints,
    NonlinearConstraints, ProblemParams,
)
from now.constraints.linear import get_linear_constraint_matrices
from now.constraints.nonlinear import evaluate_all_nonlinear


class TestBuildProblem:
    def test_returns_correct_type(self):
        c = NOW_config(N=20)
        problem = build_problem(c)
        assert isinstance(problem, OptimizationProblem)
        assert isinstance(problem.linear, LinearConstraints)
        assert isinstance(problem.nonlinear, NonlinearConstraints)
        assert isinstance(problem.params, ProblemParams)

    def test_n_vars(self):
        c = NOW_config(N=20)
        problem = build_problem(c)
        assert problem.n_vars == 3 * 20 + 1

    def test_linear_matches_matrices(self):
        c = NOW_config(N=20)
        problem = build_problem(c)
        A_ineq, b_ineq, A_eq, b_eq = get_linear_constraint_matrices(c)
        np.testing.assert_array_equal(problem.linear.A_ineq, A_ineq)
        np.testing.assert_array_equal(problem.linear.b_ineq, b_ineq)
        np.testing.assert_array_equal(problem.linear.A_eq, A_eq)
        np.testing.assert_array_equal(problem.linear.b_eq, b_eq)

    def test_nonlinear_matches(self):
        c = NOW_config(N=20)
        problem = build_problem(c)
        np.random.seed(42)
        x0 = np.random.randn(problem.n_vars)
        expected = evaluate_all_nonlinear(x0, c)
        actual = problem.nonlinear.fun(x0)
        np.testing.assert_array_equal(actual, expected)

    def test_nonlinear_jac_shape(self):
        c = NOW_config(N=20)
        problem = build_problem(c)
        np.random.seed(42)
        x0 = np.random.randn(problem.n_vars)
        jac = problem.nonlinear.jac(x0)
        assert jac.shape == (problem.nonlinear.n_constraints, problem.n_vars)

    def test_objective_callable(self):
        c = NOW_config(N=20)
        problem = build_problem(c)
        x = np.random.randn(problem.n_vars)
        f, grad = problem.objective(x)
        assert isinstance(f, float)
        assert grad.shape == (problem.n_vars,)
        assert f == -x[-1]

    def test_nonlinear_count_default(self):
        c = NOW_config(N=20)
        problem = build_problem(c)
        # 1 (tensor) + 19 (grad norm) + 3 (power) + 1 (maxwell) = 24
        assert problem.nonlinear.n_constraints == 24

    def test_nonlinear_count_maxnorm(self):
        c = NOW_config(N=20, useMaxNorm=True)
        problem = build_problem(c)
        # 1 (tensor) + 0 (no grad norm) + 3 (power) + 1 (maxwell) = 5
        assert problem.nonlinear.n_constraints == 5

    def test_nonlinear_count_motion_comp(self):
        mc = {'order': [1, 2], 'maxMagnitude': [0, 1e-4]}
        c = NOW_config(N=20, motionCompensation=mc)
        problem = build_problem(c)
        # 1 (tensor) + 19 (grad norm) + 3 (power) + 1 (maxwell) + 1 (nonlinear motion) = 25
        assert problem.nonlinear.n_constraints == 25
        assert problem.nonlinear.fun(np.zeros(problem.n_vars)).shape[0] == 25


class TestProblemParams:
    def test_from_config_preserves_values(self):
        c = NOW_config(N=30, gMax=100, sMax=200, eta=0.8)
        params = ProblemParams.from_config(c)
        assert params.N == c.N
        assert params.dt == c._dt
        assert params.gMaxConstraint == c._gMaxConstraint
        assert params.sMaxConstraint == c._sMaxConstraint
        assert params.integralConstraint == c._integralConstraint
        assert params.tolIsotropy == c._tolIsotropy
        assert params.tolMaxwell == c._tolMaxwell
        assert params.gMax == c.gMax
        assert params.totalTimeActual == c.totalTimeActual
        assert params.useMaxNorm == c.useMaxNorm
        assert params.enforceSymmetry == c.enforceSymmetry
        np.testing.assert_array_equal(params.targetTensor, c.targetTensor)
        np.testing.assert_array_equal(params.zeroGradientAtIndex, c.zeroGradientAtIndex)

    def test_from_config_copies_arrays(self):
        c = NOW_config(N=20)
        params = ProblemParams.from_config(c)
        assert params.targetTensor is not c.targetTensor
        assert params.zeroGradientAtIndex is not c.zeroGradientAtIndex

    def test_frozen(self):
        c = NOW_config(N=20)
        params = ProblemParams.from_config(c)
        import pytest
        with pytest.raises(AttributeError):
            params.N = 99
