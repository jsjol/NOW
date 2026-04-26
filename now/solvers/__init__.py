"""Solver registry."""
from __future__ import annotations

from .protocol import SolverProtocol, SolverResult
from .scipy_solver import ScipySolver

_SOLVERS: dict[str, type] = {
    'SLSQP': lambda: ScipySolver('SLSQP'),
    'trust-constr': lambda: ScipySolver('trust-constr'),
}


def get_solver(method: str) -> SolverProtocol:
    """Get a solver instance by method name."""
    if method not in _SOLVERS:
        available = ', '.join(sorted(_SOLVERS))
        raise ValueError(f"Unknown solver '{method}'. Available: {available}")
    return _SOLVERS[method]()


def register_solver(name: str, factory):
    """Register a custom solver factory."""
    _SOLVERS[name] = factory
