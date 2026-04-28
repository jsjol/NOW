"""Scipy-based solver backend (SLSQP and trust-constr)."""
from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint

from ..problem import OptimizationProblem
from .protocol import SolverResult


class ScipySolver:
    """Wraps scipy.optimize.minimize for SLSQP and trust-constr methods."""

    def __init__(self, method: str = 'SLSQP'):
        if method not in ('SLSQP', 'trust-constr'):
            raise ValueError(f"ScipySolver supports 'SLSQP' and 'trust-constr', got '{method}'")
        self.method = method

    def solve(self,
              problem: OptimizationProblem,
              x0: np.ndarray,
              options: dict | None = None) -> SolverResult:
        constraints = self._translate_constraints(problem)
        # trust-constr's BFGS update warns when consecutive gradients are
        # identical (delta_grad == 0). This is harmless — it just means the
        # quasi-Newton Hessian approximation is temporarily uninformative.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore', message='delta_grad == 0.0',
                category=UserWarning)
            res = minimize(
                problem.objective, x0,
                method=self.method,
                jac=True,
                constraints=constraints,
                options=options or {},
            )
        return SolverResult(
            x=res.x,
            fun=res.fun,
            success=res.success,
            message=res.message,
            n_iterations=getattr(res, 'nit', 0),
            raw=res,
        )

    def _translate_constraints(self, problem: OptimizationProblem) -> list:
        lin = problem.linear
        nl = problem.nonlinear

        if self.method == 'SLSQP':
            return [
                {'type': 'ineq',
                 'fun': lambda x, A=lin.A_ineq, b=lin.b_ineq: b - A @ x},
                {'type': 'eq',
                 'fun': lambda x, A=lin.A_eq, b=lin.b_eq: A @ x - b},
                {'type': 'ineq',
                 'fun': lambda x: -nl.fun(x),
                 'jac': lambda x: -nl.jac(x)},
            ]
        else:
            return [
                LinearConstraint(lin.A_ineq, lb=-np.inf, ub=lin.b_ineq),
                LinearConstraint(lin.A_eq, lb=lin.b_eq, ub=lin.b_eq),
                NonlinearConstraint(nl.fun, lb=-np.inf, ub=0., jac=nl.jac),
            ]
