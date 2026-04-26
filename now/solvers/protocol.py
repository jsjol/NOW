"""Solver protocol and result type."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

from ..problem import OptimizationProblem


@dataclass
class SolverResult:
    """Solver-agnostic optimization result."""
    x: np.ndarray
    fun: float
    success: bool
    message: str
    n_iterations: int
    raw: object = None


@runtime_checkable
class SolverProtocol(Protocol):
    def solve(self,
              problem: OptimizationProblem,
              x0: np.ndarray,
              options: dict | None = None) -> SolverResult: ...
