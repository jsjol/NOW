"""NOW — Numerical Optimization of gradient Waveforms for diffusion MRI.

Quick start::

    import numpy as np
    from now import NOW_config, now_optimize

    config = NOW_config(N=50, targetTensor=np.eye(3))
    result, config = now_optimize(config)
    print(f"b-value: {result.b:.2f} s/mm²")
"""
from .config import NOW_config
from .optimize import now_optimize, objective
from .result import NOWResult
from .constraints import get_constraints

__all__ = ['NOW_config', 'now_optimize', 'objective', 'NOWResult', 'get_constraints']
