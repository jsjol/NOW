"""Compare MATLAB and Python NOW cross-validation results.

Run from the NOW/ directory:
    python crossval/compare_results.py
"""

import sys
import os
from pathlib import Path

import numpy as np
import scipy.io


def compare(name, py_val, mat_val, atol=1e-10, rtol=1e-10):
    """Compare two arrays, return (passed, message)."""
    py_val = np.atleast_1d(np.asarray(py_val, dtype=float)).ravel()
    mat_val = np.atleast_1d(np.asarray(mat_val, dtype=float)).ravel()

    if py_val.shape != mat_val.shape:
        return False, f'{name}: SHAPE MISMATCH  python={py_val.shape}  matlab={mat_val.shape}'

    if py_val.size == 0:
        return True, f'{name}: OK (both empty)'

    max_abs_diff = np.max(np.abs(py_val - mat_val))
    # Use the same formula as np.allclose: |a-b| <= atol + rtol * |b|
    close = np.allclose(py_val, mat_val, atol=atol, rtol=rtol)

    if close:
        return True, f'{name}: OK  (max |diff| = {max_abs_diff:.2e})'
    else:
        idx = np.argmax(np.abs(py_val - mat_val))
        return False, (f'{name}: FAIL  max |diff| = {max_abs_diff:.2e} at index {idx}  '
                       f'python={py_val[idx]:.10e}  matlab={mat_val[idx]:.10e}')


def compare_matrix(name, py_val, mat_val, atol=1e-10, rtol=1e-10):
    """Compare two 2D arrays."""
    py_val = np.asarray(py_val, dtype=float)
    mat_val = np.asarray(mat_val, dtype=float)

    if py_val.shape != mat_val.shape:
        return False, f'{name}: SHAPE MISMATCH  python={py_val.shape}  matlab={mat_val.shape}'

    max_abs_diff = np.max(np.abs(py_val - mat_val))
    close = np.allclose(py_val, mat_val, atol=atol, rtol=rtol)

    if close:
        return True, f'{name}: OK  (max |diff| = {max_abs_diff:.2e}, shape={py_val.shape})'
    else:
        idx = np.unravel_index(np.argmax(np.abs(py_val - mat_val)), py_val.shape)
        return False, (f'{name}: FAIL  max |diff| = {max_abs_diff:.2e} at {idx}  '
                       f'python={py_val[idx]:.10e}  matlab={mat_val[idx]:.10e}')


def main():
    data_dir = Path(__file__).parent / 'data'
    py_path = data_dir / 'python_results.npz'
    mat_path = data_dir / 'matlab_results.mat'

    if not py_path.exists():
        print(f'ERROR: Python results not found at {py_path}')
        print('Run:  python crossval/export_python.py')
        return 1

    if not mat_path.exists():
        print(f'ERROR: MATLAB results not found at {mat_path}')
        print('Run:  export_matlab.m in MATLAB')
        return 1

    py = np.load(py_path)
    mat = scipy.io.loadmat(str(mat_path), squeeze_me=True)

    results = []
    all_pass = True

    print('=' * 70)
    print('NOW Cross-Validation: MATLAB vs Python')
    print('=' * 70)

    # --- Scalars ---
    for name in ['dt', 'gMaxConstraint', 'sMaxConstraint',
                  'integralConstraint', 'tolMaxwell', 'totalTimeActual']:
        ok, msg = compare(name, py[name], mat[name])
        results.append((ok, msg))
        if not ok:
            all_pass = False

    # --- x0 ---
    # MATLAB x0 is a column vector; squeeze handles this
    ok, msg = compare('x0', py['x0'], mat['x0'].ravel())
    results.append((ok, msg))
    if not ok:
        all_pass = False

    # --- signs ---
    py_signs = py['signs']
    mat_signs = mat['signs']
    if py_signs.size == 0 and (mat_signs.size == 0 or
                                (hasattr(mat_signs, 'shape') and 0 in mat_signs.shape)):
        results.append((True, 'signs: OK (both empty)'))
    else:
        ok, msg = compare('signs', py_signs.ravel(), np.asarray(mat_signs).ravel())
        results.append((ok, msg))
        if not ok:
            all_pass = False

    # --- zeroGradientAtIndex ---
    # Python is 0-based, MATLAB is 1-based
    py_zgi = py['zeroGradientAtIndex']
    mat_zgi = np.asarray(mat['zeroGradientAtIndex']).ravel()
    if py_zgi.size == 0 and mat_zgi.size == 0:
        results.append((True, 'zeroGradientAtIndex: OK (both empty)'))
    else:
        # Convert Python 0-based to 1-based for comparison
        py_zgi_1based = py_zgi + 1
        ok, msg = compare('zeroGradientAtIndex (py+1 vs mat)', py_zgi_1based, mat_zgi)
        results.append((ok, msg))
        if not ok:
            all_pass = False

    # --- Derivative matrices ---
    ok, msg = compare_matrix('firstDerivativeMatrix',
                             py['firstDerivativeMatrix'],
                             mat['firstDerivativeMatrix'])
    results.append((ok, msg))
    if not ok:
        all_pass = False

    ok, msg = compare_matrix('secondDerivativeMatrix',
                             py['secondDerivativeMatrix'],
                             mat['secondDerivativeMatrix'])
    results.append((ok, msg))
    if not ok:
        all_pass = False

    # --- Linear constraints ---
    for name in ['A_ineq', 'b_ineq', 'A_eq', 'b_eq']:
        py_val = py[name]
        mat_val = np.asarray(mat[name])
        if py_val.ndim == 1 and mat_val.ndim == 1:
            ok, msg = compare(name, py_val, mat_val)
        else:
            # Ensure 2D
            if py_val.ndim == 1:
                py_val = py_val.reshape(-1, 1)
            if mat_val.ndim == 1:
                mat_val = mat_val.reshape(-1, 1)
            ok, msg = compare_matrix(name, py_val, mat_val)
        results.append((ok, msg))
        if not ok:
            all_pass = False

    # --- Nonlinear constraints ---
    # Use looser tolerances for nonlinear values (floating point differences)
    ok, msg = compare('c_nonlinear', py['c_nonlinear'],
                      mat['c_nonlinear'].ravel(), atol=1e-8, rtol=1e-8)
    results.append((ok, msg))
    if not ok:
        all_pass = False

    # Jacobian: MATLAB gradc transposed to match Python layout
    py_jac = py['J_nonlinear']
    mat_jac = np.asarray(mat['J_nonlinear'])
    ok, msg = compare_matrix('J_nonlinear', py_jac, mat_jac, atol=1e-8, rtol=1e-8)
    results.append((ok, msg))
    if not ok:
        all_pass = False

    # --- Objective ---
    ok, msg = compare('fval', py['fval'], mat['fval'])
    results.append((ok, msg))
    if not ok:
        all_pass = False

    ok, msg = compare('grad', py['grad'], mat['grad'].ravel())
    results.append((ok, msg))
    if not ok:
        all_pass = False

    # --- Print report ---
    print()
    for ok, msg in results:
        status = 'PASS' if ok else 'FAIL'
        print(f'  [{status}] {msg}')

    print()
    n_pass = sum(1 for ok, _ in results if ok)
    n_total = len(results)
    print(f'Result: {n_pass}/{n_total} checks passed.')

    if all_pass:
        print('ALL CHECKS PASSED.')
        return 0
    else:
        print('SOME CHECKS FAILED.')
        return 1


if __name__ == '__main__':
    sys.exit(main())
