import numpy as np
from scipy.optimize import NonlinearConstraint

from ..utils import build_first_derivative_matrix, build_integration_weights


def _unpack(x, N):
    """Split optimization vector into q-trajectory and scaling parameter."""
    q = x[:-1]
    s = x[-1]
    Q = q.reshape(N, 3, order='F')
    return Q, s


def evaluate_all_nonlinear(x, config):
    """Evaluate all nonlinear constraint values. Returns c vector (fmincon convention: c <= 0).

    This matches MATLAB's nonlconAnalytic return order:
    [c1, c2(N-1), c3, c4, c5, c6, c7(...)]
    """
    N = config.N
    Q, s = _unpack(x, N)

    A1 = build_first_derivative_matrix(N)
    w = build_integration_weights(N)
    g = A1 @ Q

    # Tensor encoding (c1)
    wQ = w * Q
    B = Q.T @ wQ
    c1 = np.trace((B - s * config.targetTensor).T @ (B - s * config.targetTensor)) \
         - (s * config._tolIsotropy) ** 2

    # Gradient norm (c2) — only if not using max-norm
    if not config.useMaxNorm:
        c2 = (np.sum(g ** 2, axis=1) - config._gMaxConstraint ** 2)
    else:
        c2 = np.array([])

    # Power constraints (c3, c4, c5)
    c3 = g[:, 0] @ g[:, 0] - config._integralConstraint
    c4 = g[:, 1] @ g[:, 1] - config._integralConstraint
    c5 = g[:, 2] @ g[:, 2] - config._integralConstraint

    # Maxwell compensation (c6)
    tolMaxwell_scaled = config._tolMaxwell * config._dt ** 2
    if np.isinf(tolMaxwell_scaled):
        c6 = np.array([])
    else:
        signedg = g * config.signs
        M = g.T @ signedg
        m = np.sqrt(np.trace(M.T @ M))
        c6 = np.array([m - tolMaxwell_scaled])

    # Nonlinear motion compensation (c7)
    mc = config.motionCompensation
    if len(mc['linear']) > 0:
        nonlinear_ind = np.where(~mc['linear'])[0]
    else:
        nonlinear_ind = np.array([], dtype=int)

    if len(nonlinear_ind) == 0:
        c7 = np.array([])
    else:
        t = (np.arange(1, N + 1) - 0.5) * config._dt
        gamma_rad = 2.6751e+08
        c7 = np.zeros(len(nonlinear_ind))
        for i, ni in enumerate(nonlinear_ind):
            order = mc['order'][ni]
            moment_weighting = -order * config._dt * t ** (order - 1)
            moment_vector = moment_weighting @ Q
            c7[i] = np.sum(moment_vector ** 2) - \
                    (mc['maxMagnitude'][ni] * 1000 ** order / (gamma_rad * 1e-6)) ** 2

    return np.concatenate([np.atleast_1d(c1), c2,
                           np.atleast_1d(c3), np.atleast_1d(c4), np.atleast_1d(c5),
                           c6, c7])


def evaluate_all_nonlinear_jacobian(x, config):
    """Evaluate analytical Jacobian of all nonlinear constraints.

    Returns (n_constraints, 3N+1) matrix. Matches MATLAB gradc' (transposed).
    """
    N = config.N
    Q, s = _unpack(x, N)

    A1 = build_first_derivative_matrix(N)
    w = build_integration_weights(N)
    g = A1 @ Q
    wQ = w * Q
    B = Q.T @ wQ

    # --- Tensor encoding Jacobian (dc1/dx) ---
    dc1_dB = 2 * B - 2 * s * config.targetTensor
    firstTerm = np.kron(np.eye(3), wQ.T)
    secondTerm = firstTerm.reshape(3, 3, 3 * N, order='F')
    secondTerm = secondTerm.transpose(1, 0, 2)
    secondTerm = secondTerm.reshape(9, 3 * N)
    dB_dq = firstTerm + secondTerm

    dc1_dq = dc1_dB.reshape(1, 9) @ dB_dq
    Bt = config.targetTensor
    dc1_ds = 2 * s * np.trace(Bt.T @ Bt) - 2 * np.trace(B.T @ Bt) - 2 * s * config._tolIsotropy ** 2
    dc1_dx = np.concatenate([dc1_dq.ravel(), [dc1_ds]])

    # --- Gradient norm Jacobian (dc2/dx) ---
    if not config.useMaxNorm:
        dc2_dq = 2 * np.hstack([
            A1 * g[:, 0:1],
            A1 * g[:, 1:2],
            A1 * g[:, 2:3],
        ])
        dc2_dx = np.c_[dc2_dq, np.zeros((N - 1, 1))]
    else:
        dc2_dx = np.zeros((0, 3 * N + 1))

    # --- Power Jacobians (dc3/dx, dc4/dx, dc5/dx) ---
    dc3_dx = np.zeros(3 * N + 1)
    dc3_dx[:N] = 2 * g[:, 0] @ A1

    dc4_dx = np.zeros(3 * N + 1)
    dc4_dx[N:2 * N] = 2 * g[:, 1] @ A1

    dc5_dx = np.zeros(3 * N + 1)
    dc5_dx[2 * N:3 * N] = 2 * g[:, 2] @ A1

    # --- Maxwell Jacobian (dc6/dx) ---
    tolMaxwell_scaled = config._tolMaxwell * config._dt ** 2
    if np.isinf(tolMaxwell_scaled):
        dc6_dx = np.zeros((0, 3 * N + 1))
    else:
        signedg = g * config.signs
        M = g.T @ signedg
        m = np.sqrt(np.trace(M.T @ M))
        dc6_dM = (1 / m) * M

        weightedQ_maxwell = A1.T @ signedg
        firstTerm_m = np.kron(np.eye(3), weightedQ_maxwell.T)
        secondTerm_m = firstTerm_m.reshape(3, 3, 3 * N, order='F')
        secondTerm_m = secondTerm_m.transpose(1, 0, 2)
        secondTerm_m = secondTerm_m.reshape(9, 3 * N)
        dM_dq = firstTerm_m + secondTerm_m

        dc6_dq = dc6_dM.reshape(1, 9) @ dM_dq
        dc6_dx = np.concatenate([dc6_dq.ravel(), [0]]).reshape(1, -1)

    # --- Motion compensation Jacobian (dc7/dx) ---
    mc = config.motionCompensation
    if len(mc['linear']) > 0:
        nonlinear_ind = np.where(~mc['linear'])[0]
    else:
        nonlinear_ind = np.array([], dtype=int)

    if len(nonlinear_ind) == 0:
        dc7_dx = np.zeros((0, 3 * N + 1))
    else:
        t = (np.arange(1, N + 1) - 0.5) * config._dt
        dc7_dx = np.zeros((len(nonlinear_ind), 3 * N + 1))
        for i, ni in enumerate(nonlinear_ind):
            order = mc['order'][ni]
            moment_weighting = -order * config._dt * t ** (order - 1)
            moment_vector = moment_weighting @ Q
            dc7_dx[i, :3 * N] = 2 * np.kron(moment_vector, moment_weighting)

    # Stack all rows
    jac = np.vstack([
        dc1_dx.reshape(1, -1),
        dc2_dx,
        dc3_dx.reshape(1, -1),
        dc4_dx.reshape(1, -1),
        dc5_dx.reshape(1, -1),
        dc6_dx,
        dc7_dx,
    ])
    return jac


def get_nonlinear_constraints(config, method=None):
    """Build scipy-compatible nonlinear constraints."""

    def fun(x):
        return evaluate_all_nonlinear(x, config)

    def jac(x):
        return evaluate_all_nonlinear_jacobian(x, config)

    n_constraints = _count_nonlinear_constraints(config)

    if method == 'SLSQP':
        return [{'type': 'ineq', 'fun': lambda x: -fun(x), 'jac': lambda x: -jac(x)}]
    else:
        return [NonlinearConstraint(fun, lb=-np.inf, ub=0.,
                                    jac=jac)]


def _count_nonlinear_constraints(config):
    n = 1  # tensor encoding
    if not config.useMaxNorm:
        n += config.N - 1  # gradient norm
    n += 3  # power (x, y, z)
    if not np.isinf(config._tolMaxwell * config._dt ** 2):
        n += 1  # Maxwell
    mc = config.motionCompensation
    if len(mc['linear']) > 0:
        n += np.sum(~mc['linear'])
    return n
