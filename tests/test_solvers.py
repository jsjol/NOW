"""Tests for solver backends."""
import numpy as np
import pytest
from now.config import NOW_config
from now.problem import build_problem
from now.solvers import get_solver, register_solver
from now.solvers.protocol import SolverResult, SolverProtocol
from now.solvers.scipy_solver import ScipySolver


class TestSolverRegistry:
    def test_get_slsqp(self):
        solver = get_solver('SLSQP')
        assert isinstance(solver, ScipySolver)
        assert solver.method == 'SLSQP'

    def test_get_trust_constr(self):
        solver = get_solver('trust-constr')
        assert isinstance(solver, ScipySolver)
        assert solver.method == 'trust-constr'

    def test_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown solver"):
            get_solver('nonexistent')

    def test_register_custom(self):
        class DummySolver:
            def solve(self, problem, x0, options=None):
                return SolverResult(x=x0, fun=0.0, success=True,
                                    message='dummy', n_iterations=0)

        register_solver('dummy', lambda: DummySolver())
        solver = get_solver('dummy')
        assert isinstance(solver, DummySolver)


class TestScipySolver:
    def test_invalid_method(self):
        with pytest.raises(ValueError, match="ScipySolver supports"):
            ScipySolver('bogus')

    def test_implements_protocol(self):
        solver = ScipySolver('SLSQP')
        assert isinstance(solver, SolverProtocol)

    @pytest.mark.slow
    def test_slsqp_solves(self):
        c = NOW_config(N=15)
        problem = build_problem(c)
        solver = ScipySolver('SLSQP')
        np.random.seed(42)
        x0 = np.random.randn(problem.n_vars) * 0.1
        result = solver.solve(problem, x0,
                              options={'maxiter': 500, 'ftol': 1e-10})
        assert isinstance(result, SolverResult)
        assert result.x.shape == (problem.n_vars,)
        assert result.raw is not None

    @pytest.mark.slow
    def test_trust_constr_solves(self):
        c = NOW_config(N=15)
        problem = build_problem(c)
        solver = ScipySolver('trust-constr')
        np.random.seed(42)
        x0 = np.random.randn(problem.n_vars) * 0.1
        result = solver.solve(problem, x0,
                              options={'maxiter': 500})
        assert isinstance(result, SolverResult)
        assert result.x.shape == (problem.n_vars,)
