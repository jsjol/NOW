"""Constraint assembly for NOW optimization."""
from .linear import get_linear_constraint_matrices
from .nonlinear import evaluate_all_nonlinear, evaluate_all_nonlinear_jacobian
