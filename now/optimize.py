import numpy as np
import time
from scipy.optimize import minimize

from .config import NOW_config
from .constraints import get_constraints
from .result import build_result
from .utils import build_first_derivative_matrix, build_second_derivative_matrix


def objective(x):
    """Minimize -s (maximize b-value). Returns (f, grad)."""
    f = -x[-1]
    grad = np.zeros_like(x)
    grad[-1] = -1
    return f, grad


def get_initial_guess(config, iteration):
    """Generate initial guess, matching MATLAB getInitialGuess."""
    if config.initialGuess == 'user-provided' and config.x0 is not None:
        if iteration > 0:
            x0 = np.random.randn(3 * config.N + 1)
        else:
            return config.x0.copy()
    else:
        x0 = np.random.randn(3 * config.N + 1)

    N = config.N
    # Scale by target tensor diagonal (MATLAB uses linear indexing: (1,1), (2,2), (3,3))
    x0[:N] *= config.targetTensor[0, 0]
    x0[N:2 * N] *= config.targetTensor[1, 1]
    x0[2 * N:3 * N] *= config.targetTensor[2, 2]
    return x0


def now_optimize(config=None, method='SLSQP', max_attempts=10, verbose=True):
    """Run the NOW optimization.

    Parameters
    ----------
    config : NOW_config or None
        Problem configuration. If None, uses defaults.
    method : str
        Scipy optimizer method: 'SLSQP' or 'trust-constr'.
    max_attempts : int
        Maximum optimization attempts with random restarts.
    verbose : bool
        Print progress information.

    Returns
    -------
    result : NOWResult
    config : NOW_config
    """
    if config is None:
        config = NOW_config()

    constraints = get_constraints(config, method=method)

    success = False
    best_x = None
    best_fval = np.inf
    total_time = 0

    for iteration in range(max_attempts):
        try:
            x0 = get_initial_guess(config, iteration)

            if verbose:
                if iteration == 0:
                    print(f'Optimizing {config.name}')
                else:
                    print(f'Optimizing {config.name}, attempt {iteration + 1}')

            t0 = time.perf_counter()

            if method == 'SLSQP':
                res = minimize(objective, x0, method='SLSQP', jac=True,
                               constraints=constraints,
                               options={'maxiter': config.MaxIter, 'disp': verbose,
                                        'ftol': 1e-10})
            else:
                res = minimize(objective, x0, method='trust-constr', jac=True,
                               constraints=constraints,
                               options={'maxiter': config.MaxIter, 'disp': verbose})

            elapsed = time.perf_counter() - t0
            total_time += elapsed

            if verbose:
                print(f'Optimization took {elapsed:.3f}s.')

            if res.success or res.fun < best_fval:
                best_x = res.x
                best_fval = res.fun
                if res.success:
                    success = True

            if success:
                break

            if not config.redoIfFailed:
                if verbose and not success:
                    print('Optimization failed but will not be repeated!')
                break

        except Exception as e:
            if verbose:
                print(f'Attempt {iteration + 1} failed: {e}')
            continue

    if best_x is None:
        raise RuntimeError(f'All {max_attempts} optimization attempts failed.')

    A1 = build_first_derivative_matrix(config.N)
    A2 = build_second_derivative_matrix(config.N)
    x0_final = get_initial_guess(config, 0)

    result = build_result(best_x, x0_final, config, A1, A2,
                          optimization_time=total_time,
                          n_iter=iteration + 1,
                          optimizer_output=res if best_x is not None else None)

    return result, config
