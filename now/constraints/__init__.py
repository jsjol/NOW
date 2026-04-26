from .linear import get_linear_constraints
from .nonlinear import get_nonlinear_constraints


def get_constraints(config, method=None):
    return get_linear_constraints(config, method=method) + \
           get_nonlinear_constraints(config, method=method)
