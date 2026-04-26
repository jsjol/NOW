import numpy as np
from dataclasses import dataclass, field
from .utils import gamma as now_gamma_hz


@dataclass
class NOWResult:
    """Result of a NOW optimization.

    The fields ``gwf``, ``rf``, and ``dt`` are compatible with the md-dMRI
    framework (https://github.com/markus-nilsson/md-dmri).

    Attributes
    ----------
    q : (N, 3) ndarray — q-space trajectory [1/m]
    g : (N+2, 3) ndarray — gradient waveform [mT/m], zero-padded
    slew : (N+2, 3) ndarray — slew rate [T/m/s], zero-padded
    b : float — b-value [s/mm²]
    B : (3, 3) ndarray — full b-tensor [s/m²]
    kappa : float — encoding efficiency
    etaOpt : float — realised energy/efficacy ratio
    gwf : (N+2, 3) ndarray — gradient waveform [T/m], md-dMRI compatible
    rf : (N+2, 1) ndarray — spin dephasing direction
    dt : float — time step [s]
    """
    q: np.ndarray
    g: np.ndarray
    slew: np.ndarray
    b: float
    B: np.ndarray
    kappa: float
    etaOpt: float
    q0: np.ndarray
    rawq: np.ndarray
    zind: np.ndarray
    rf: np.ndarray
    gwf: np.ndarray
    dt: float
    optimizationTime: float = 0.0
    iter: int = 0
    optimizerOutput: object = None


def build_result(x, x0, config, A1, A2, optimization_time=0.0, n_iter=0, optimizer_output=None):
    """Build NOWResult from optimization solution, matching MATLAB optimize.m post-processing."""
    gamma = now_gamma_hz() * 2 * np.pi  # rad/(s*T)
    N = config.N

    q_raw = x[:3 * N].reshape(N, 3, order='F')
    g = np.vstack([np.zeros((1, 3)), A1 @ q_raw, np.zeros((1, 3))]) / config._dt
    slew = np.vstack([np.zeros((1, 3)), A2 @ q_raw, np.zeros((1, 3))]) / config._dt ** 2

    q = gamma * 1e-6 * q_raw  # SI units
    q0_raw = x0[:3 * N].reshape(N, 3, order='F')
    q0 = gamma * 1e-6 * q0_raw

    B = config._dt * 1e-3 * (q.T @ q)
    b_val = np.trace(B) * 1e-6  # s/mm^2

    C = b_val / (1e-6 * gamma ** 2 * (config._gMaxConstraint * 1e-3 / config._dt) ** 2
                 * (config.totalTimeActual * 1e-3) ** 3)
    kappa = 4 * C
    etaOpt = config._dt * np.max(np.diag(g.T @ g)) / \
             (config._gMaxConstraint ** 2 * config.totalTimeActual)

    # md-dMRI compatible output
    rf = np.ones((g.shape[0], 1))
    if len(config.zeroGradientAtIndex) > 0 and config.signs is not None:
        signs_flat = config.signs.ravel()
        rf = np.concatenate([[signs_flat[0]], signs_flat, [signs_flat[-1]]]).reshape(-1, 1)
        rf = rf / np.max(np.abs(rf))

    gwf = (g / 1000) * rf  # T/m

    return NOWResult(
        q=q, g=g, slew=slew, b=b_val, B=B, kappa=kappa, etaOpt=etaOpt,
        q0=q0, rawq=x[:3 * N], zind=config.zeroGradientAtIndex,
        rf=rf, gwf=gwf, dt=config._dt / 1000,
        optimizationTime=optimization_time, iter=n_iter,
        optimizerOutput=optimizer_output)
