import numpy as np
from scipy import sparse


NOW_GAMMA = 42.576e6  # Hz/T for hydrogen


def gamma():
    return NOW_GAMMA


def build_first_derivative_matrix(N):
    """First derivative (forward difference) matrix, (N-1) x N.
    Matches MATLAB: -diag(ones(N,1)) + diag(ones(N-1,1),1), then rows 1:end-1.
    """
    D = -np.eye(N) + np.diag(np.ones(N - 1), 1)
    return D[:-1, :]


def build_second_derivative_matrix(N):
    """Second derivative (central difference) matrix, N x N.
    Matches MATLAB: diag(ones(N-1,1),-1) - 2*diag(ones(N,1)) + diag(ones(N-1,1),1).
    """
    return np.diag(np.ones(N - 1), -1) - 2 * np.eye(N) + np.diag(np.ones(N - 1), 1)


def build_integration_weights(N):
    """Trapezoidal integration weights, shape (N, 1). Normalized by (N-1)."""
    w = np.ones((N, 1))
    w[0] = 0.5
    w[-1] = 0.5
    w /= (N - 1)
    return w


def gwf_to_q(gwf, dt):
    """Convert gradient waveform to q-vector via cumulative integration."""
    return np.cumsum(gwf, axis=0) * dt


def maxwell_coeff(gwf, rf, dt):
    """Compute Maxwell coefficient and index from gradient waveform.

    Parameters
    ----------
    gwf : (T, 3) array — gradient waveform in T/m
    rf : (T,) array — spin dephasing direction (+1/-1)
    dt : float — time step in seconds

    Returns
    -------
    M : (3, 3) array — Maxwell matrix
    maxwell_index : float — Frobenius norm of M
    """
    signed_gwf = gwf * rf[:, None]
    M = gwf.T @ signed_gwf * dt
    maxwell_index = np.sqrt(np.trace(M.T @ M))
    return M, maxwell_index


def write_wf(result, filename):
    """Write gradient waveform to text file(s), compatible with md-dMRI."""
    gwf = result.gwf
    rf = result.rf
    dt = result.dt

    sign_change = np.where(np.diff(rf.ravel()) != 0)[0]

    if len(sign_change) > 0:
        split_idx = sign_change[len(sign_change) // 2] + 1
        gwf_a = gwf[:split_idx]
        gwf_b = gwf[split_idx:]
        np.savetxt(filename + '_A.txt', gwf_a, fmt='%.6e')
        np.savetxt(filename + '_B.txt', gwf_b, fmt='%.6e')
    else:
        np.savetxt(filename + '.txt', gwf, fmt='%.6e')


def problem_to_name(config):
    """Generate a descriptive filename from problem parameters."""
    tensor = config.targetTensor
    eigvals = np.sort(np.linalg.eigvalsh(tensor))[::-1]

    if np.allclose(eigvals, [1, 1, 1]):
        shape = 'STE'
    elif np.allclose(eigvals, [1, 1, 0]):
        shape = 'PTE'
    elif np.allclose(eigvals, [1, 0, 0]):
        shape = 'LTE'
    else:
        shape = 'custom'

    name = (f"{shape}_N{config.N}_gMax{config.gMax}_sMax{config.sMax}"
            f"_pre{config.durationFirstPartRequested}"
            f"_post{config.durationSecondPartRequested}"
            f"_pause{config.durationZeroGradientRequested}")
    return name
