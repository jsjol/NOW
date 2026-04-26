"""Export Python NOW results for cross-validation against MATLAB.

Run from the NOW/ directory:
    python crossval/export_python.py
"""

import sys
import os

# Ensure the NOW package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from pathlib import Path

from now.config import NOW_config
from now.utils import build_first_derivative_matrix, build_second_derivative_matrix
from now.constraints.linear import get_linear_constraint_matrices
from now.constraints.nonlinear import evaluate_all_nonlinear, evaluate_all_nonlinear_jacobian
from now.optimize import objective


def main():
    # Create config with defaults (N=77, default timings)
    config = NOW_config()
    N = config.N

    # Deterministic x0 — linspace avoids RNG mismatch between MATLAB and Python
    x0 = np.linspace(-1, 1, 3 * N + 1)

    # Apply the same target tensor scaling as MATLAB's getInitialGuess.
    # MATLAB uses targetTensor(1), (5), (9) = diagonal elements (1,1), (2,2), (3,3).
    # MATLAB index ranges: (1:N), (N+1:2N+1), (2N+1:end-1) — element 2N+1 scaled twice.
    # With identity tensor (default) all factors are 1 so x0 is unchanged.
    x0[:N] *= config.targetTensor[0, 0]
    x0[N:2 * N + 1] *= config.targetTensor[1, 1]
    x0[2 * N:3 * N] *= config.targetTensor[2, 2]

    # Config-derived values
    dt = config._dt
    gMaxConstraint = config._gMaxConstraint
    sMaxConstraint = config._sMaxConstraint
    integralConstraint = config._integralConstraint
    tolMaxwell = config._tolMaxwell
    signs = config.signs if config.signs is not None else np.array([])
    zeroGradientAtIndex = config.zeroGradientAtIndex
    totalTimeActual = config.totalTimeActual

    # Derivative matrices
    firstDerivativeMatrix = build_first_derivative_matrix(N)
    secondDerivativeMatrix = build_second_derivative_matrix(N)

    # Linear constraint matrices
    A_ineq, b_ineq, A_eq, b_eq = get_linear_constraint_matrices(config)

    # Nonlinear constraints at x0
    c_nonlinear = evaluate_all_nonlinear(x0, config)
    J_nonlinear = evaluate_all_nonlinear_jacobian(x0, config)

    # Objective at x0
    fval, grad = objective(x0)

    # Save to npz
    out_dir = Path(__file__).parent / 'data'
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / 'python_results.npz'

    np.savez(
        out_path,
        x0=x0,
        dt=dt,
        gMaxConstraint=gMaxConstraint,
        sMaxConstraint=sMaxConstraint,
        integralConstraint=integralConstraint,
        tolMaxwell=tolMaxwell,
        signs=signs,
        zeroGradientAtIndex=zeroGradientAtIndex,
        totalTimeActual=totalTimeActual,
        firstDerivativeMatrix=firstDerivativeMatrix,
        secondDerivativeMatrix=secondDerivativeMatrix,
        A_ineq=A_ineq,
        b_ineq=b_ineq,
        A_eq=A_eq,
        b_eq=b_eq,
        c_nonlinear=c_nonlinear,
        J_nonlinear=J_nonlinear,
        fval=fval,
        grad=grad,
    )

    print(f'Exported Python results to {out_path}')
    print(f'  N = {N}')
    print(f'  dt = {dt}')
    print(f'  x0 shape = {x0.shape}')
    print(f'  c_nonlinear shape = {c_nonlinear.shape}')
    print(f'  J_nonlinear shape = {J_nonlinear.shape}')
    print(f'  A_ineq shape = {A_ineq.shape}')
    print(f'  A_eq shape = {A_eq.shape}')
    print(f'  zeroGradientAtIndex (0-based) = {zeroGradientAtIndex}')


if __name__ == '__main__':
    main()
