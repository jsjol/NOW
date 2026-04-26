"""NOW — Numerical Optimization of gradient Waveforms for diffusion MRI.

Quick start::

    import numpy as np
    from now import NOW_config, now_optimize

    config = NOW_config(N=50, targetTensor=np.eye(3))
    result, config = now_optimize(config)
    print(f"b-value: {result.b:.2f} s/mm²")
"""
from .config import NOW_config
from .optimize import now_optimize
from .problem import objective, OptimizationProblem, build_problem
from .result import NOWResult
from .solvers import get_solver, register_solver, SolverResult
from .io import save_config, load_config, export_problem

__all__ = [
    'NOW_config', 'now_optimize', 'objective', 'NOWResult',
    'OptimizationProblem', 'build_problem',
    'get_solver', 'register_solver', 'SolverResult',
    'save_config', 'load_config', 'export_problem',
]
