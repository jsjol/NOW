"""Export OptimizationProblem data to npz or mat format."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.io import savemat

from ..problem import OptimizationProblem


def export_problem(problem: OptimizationProblem,
                   path: str | Path,
                   fmt: str = 'npz') -> None:
    """Export the linear constraint matrices and problem metadata.

    Nonlinear constraints are callables and cannot be serialized directly.
    They are evaluated at a zero vector and included as sample values.

    Parameters
    ----------
    problem : OptimizationProblem
    path : str or Path
    fmt : 'npz' or 'mat'
    """
    path = Path(path)
    x0 = np.zeros(problem.n_vars)

    data = {
        'n_vars': np.array(problem.n_vars),
        'A_ineq': problem.linear.A_ineq,
        'b_ineq': problem.linear.b_ineq,
        'A_eq': problem.linear.A_eq,
        'b_eq': problem.linear.b_eq,
        'nonlinear_at_zero': problem.nonlinear.fun(x0),
        'nonlinear_jac_at_zero': problem.nonlinear.jac(x0),
        'n_nonlinear': np.array(problem.nonlinear.n_constraints),
    }

    if fmt == 'mat':
        savemat(str(path), data)
    elif fmt == 'npz':
        np.savez(str(path), **data)
    else:
        raise ValueError(f"Unsupported format '{fmt}'. Use 'npz' or 'mat'.")
