"""Solver-independent optimization problem representation."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .constraints.linear import get_linear_constraint_matrices
from .constraints.nonlinear import evaluate_all_nonlinear, evaluate_all_nonlinear_jacobian


def objective(x):
    """Minimize -s (maximize b-value). Returns (f, grad)."""
    f = -x[-1]
    grad = np.zeros_like(x)
    grad[-1] = -1
    return f, grad


@dataclass(frozen=True)
class LinearConstraints:
    """Linear constraints in standard form: A_ineq @ x <= b_ineq, A_eq @ x == b_eq."""
    A_ineq: np.ndarray
    b_ineq: np.ndarray
    A_eq: np.ndarray
    b_eq: np.ndarray


@dataclass(frozen=True)
class NonlinearConstraints:
    """Nonlinear constraints: fun(x) <= 0, with analytical Jacobian."""
    fun: Callable[[np.ndarray], np.ndarray]
    jac: Callable[[np.ndarray], np.ndarray]
    n_constraints: int


@dataclass(frozen=True)
class ProblemParams:
    """Derived parameters needed by constraints and result building.

    Extracted from NOW_config to decouple the solver from the config class.
    """
    N: int
    dt: float
    gMaxConstraint: float
    sMaxConstraint: float
    integralConstraint: float
    tolIsotropy: float
    tolMaxwell: float
    targetTensor: np.ndarray
    useMaxNorm: bool
    signs: np.ndarray | None
    zeroGradientAtIndex: np.ndarray
    motionCompensation: dict
    totalTimeActual: float
    gMax: float
    enforceSymmetry: bool

    @staticmethod
    def from_config(config) -> ProblemParams:
        return ProblemParams(
            N=config.N,
            dt=config._dt,
            gMaxConstraint=config._gMaxConstraint,
            sMaxConstraint=config._sMaxConstraint,
            integralConstraint=config._integralConstraint,
            tolIsotropy=config._tolIsotropy,
            tolMaxwell=config._tolMaxwell,
            targetTensor=config.targetTensor.copy(),
            useMaxNorm=config.useMaxNorm,
            signs=config.signs.copy() if config.signs is not None else None,
            zeroGradientAtIndex=config.zeroGradientAtIndex.copy(),
            motionCompensation=config.motionCompensation,
            totalTimeActual=config.totalTimeActual,
            gMax=config.gMax,
            enforceSymmetry=config.enforceSymmetry,
        )


@dataclass(frozen=True)
class OptimizationProblem:
    """Complete, solver-independent NLP for gradient waveform optimization."""
    n_vars: int
    objective: Callable[[np.ndarray], tuple[float, np.ndarray]]
    linear: LinearConstraints
    nonlinear: NonlinearConstraints
    params: ProblemParams


def _count_nonlinear(config) -> int:
    """Count nonlinear constraints for a given config."""
    n = 1  # tensor encoding
    if not config.useMaxNorm:
        n += config.N - 1
    n += 3  # power
    if not np.isinf(config._tolMaxwell * config._dt ** 2):
        n += 1
    mc = config.motionCompensation
    if len(mc['linear']) > 0:
        n += int(np.sum(~mc['linear']))
    return n


def build_problem(config) -> OptimizationProblem:
    """Build a solver-independent OptimizationProblem from a NOW_config."""
    N = config.N
    n_vars = 3 * N + 1

    A_ineq, b_ineq, A_eq, b_eq = get_linear_constraint_matrices(config)
    linear = LinearConstraints(A_ineq=A_ineq, b_ineq=b_ineq,
                               A_eq=A_eq, b_eq=b_eq)

    def nl_fun(x):
        return evaluate_all_nonlinear(x, config)

    def nl_jac(x):
        return evaluate_all_nonlinear_jacobian(x, config)

    nonlinear = NonlinearConstraints(
        fun=nl_fun, jac=nl_jac,
        n_constraints=_count_nonlinear(config),
    )

    params = ProblemParams.from_config(config)

    return OptimizationProblem(
        n_vars=n_vars,
        objective=objective,
        linear=linear,
        nonlinear=nonlinear,
        params=params,
    )
