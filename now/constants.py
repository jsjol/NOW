"""Physical constants for diffusion MRI gradient waveform design."""
import numpy as np

GAMMA_HZ: float = 42.576e6
"""Gyromagnetic ratio for hydrogen [Hz/T]."""

GAMMA_RAD: float = GAMMA_HZ * 2 * np.pi
"""Gyromagnetic ratio for hydrogen [rad/(s·T)]."""
