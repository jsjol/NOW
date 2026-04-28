import numpy as np

from ..utils import build_first_derivative_matrix, build_second_derivative_matrix


def get_linear_constraint_matrices(config):
    """Return raw (A_ineq, b_ineq, A_eq, b_eq) for testing against MATLAB fixtures."""
    N = config.N
    A1_mat = build_first_derivative_matrix(N)
    A2_mat = build_second_derivative_matrix(N)

    # Inequality
    if config.useMaxNorm:
        A_grad = np.kron(np.eye(3), A1_mat)
        A_grad = np.c_[A_grad, np.zeros((A_grad.shape[0], 1))]
        b_grad = config._gMaxConstraint * np.ones(A_grad.shape[0])
        A1_block = A_grad
        b1_block = b_grad
    else:
        A1_block = np.zeros((0, 3 * N + 1))
        b1_block = np.array([])

    A_slew = np.kron(np.eye(3), A2_mat)
    A_slew = np.c_[A_slew, np.zeros((A_slew.shape[0], 1))]
    b_slew = config._sMaxConstraint * np.ones(A_slew.shape[0])

    A_ineq = np.vstack([A1_block, -A1_block, A_slew, -A_slew]) if A1_block.shape[0] > 0 \
        else np.vstack([A_slew, -A_slew])
    b_ineq = np.concatenate([b1_block, b1_block, b_slew, b_slew]) if len(b1_block) > 0 \
        else np.concatenate([b_slew, b_slew])

    # Equality constraints
    mc = config.motionCompensation
    linear_ind = np.where(mc['linear'])[0] if len(mc['linear']) > 0 else np.array([], dtype=int)

    n_eq_rows = (2
                 + len(config.zeroGradientAtIndex)
                 + len(linear_ind)
                 + (1 if config.doBackgroundCompensation > 0 else 0))

    Aeq_single = np.zeros((n_eq_rows, N))
    Aeq_single[0, 0] = 1
    Aeq_single[1, N - 1] = 1

    if len(config.zeroGradientAtIndex) > 0:
        Aeq_single[2:2 + len(config.zeroGradientAtIndex), :] = \
            A1_mat[config.zeroGradientAtIndex.astype(int), :]

    row_offset = 2 + len(config.zeroGradientAtIndex)
    if len(linear_ind) > 0:
        t = (np.arange(1, N + 1) - 0.5) * config._dt
        for i, li in enumerate(linear_ind):
            order = mc['order'][li]
            Aeq_single[row_offset + i, :] = -order * config._dt * t ** (order - 1)

    if config.doBackgroundCompensation > 0:
        row_bg = row_offset + len(linear_ind)
        s = config.startTime
        signs_with_first = np.concatenate([[1], config.signs.ravel()])
        H = np.cumsum(signs_with_first) * config._dt + s
        Aeq_single[row_bg, 0] = H[0] / 2
        Aeq_single[row_bg, 1:-1] = H[1:-1]
        Aeq_single[row_bg, -1] = H[-1] / 2

    if config.enforceSymmetry:
        if len(config.zeroGradientAtIndex) == 0:
            indicesBefore = np.arange(0, N // 2)
            indicesAfter = np.arange(N // 2, N // 2 + (N - N // 2))
        else:
            first_zero_matlab = config.zeroGradientAtIndex[0] + 1
            last_zero_matlab = config.zeroGradientAtIndex[-1] + 1
            indicesBefore = np.arange(0, first_zero_matlab)
            indicesAfter = np.arange(last_zero_matlab, N)

        if len(indicesBefore) != len(indicesAfter):
            raise ValueError('Cannot enforce symmetry.')

        Nact = len(indicesBefore)
        sym_rows = np.zeros((Nact, N))
        sym_rows[:, indicesBefore] = np.eye(Nact)[:, ::-1]
        sym_rows[:, indicesAfter] = -np.eye(Nact)
        Aeq_single = np.vstack([Aeq_single, sym_rows])

    Aeq = np.kron(np.eye(3), Aeq_single)
    Aeq = np.c_[Aeq, np.zeros((Aeq.shape[0], 1))]
    beq = np.zeros(Aeq.shape[0])

    return A_ineq, b_ineq, Aeq, beq
