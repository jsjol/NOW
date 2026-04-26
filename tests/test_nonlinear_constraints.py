"""Tests for nonlinear constraints and their Jacobians."""
import numpy as np
import pytest
from now.config import NOW_config
from now.constraints.nonlinear import (
    evaluate_all_nonlinear,
    evaluate_all_nonlinear_jacobian,
)
from .helpers import load_fixture, skip_without_fixture


class TestNonlinearConstraintValues:
    """Test nonlinear constraint function values."""

    def test_output_length_default(self):
        c = NOW_config(N=20)
        x = np.random.randn(3 * 20 + 1)
        vals = evaluate_all_nonlinear(x, c)
        # 1 (tensor) + 19 (grad norm) + 3 (power) + 1 (maxwell) = 24
        assert len(vals) == 24

    def test_output_length_maxnorm(self):
        c = NOW_config(N=20, useMaxNorm=True)
        x = np.random.randn(3 * 20 + 1)
        vals = evaluate_all_nonlinear(x, c)
        # 1 (tensor) + 0 (no grad norm) + 3 (power) + 1 (maxwell) = 5
        assert len(vals) == 5

    def test_output_length_no_maxwell(self):
        c = NOW_config(N=20, doMaxwellComp=False)
        x = np.random.randn(3 * 20 + 1)
        vals = evaluate_all_nonlinear(x, c)
        # 1 (tensor) + 19 (grad norm) + 3 (power) + 0 (no maxwell) = 23
        assert len(vals) == 23

    def test_output_length_motion_comp(self):
        mc = {'order': [1, 2], 'maxMagnitude': [0, 1e-4]}
        c = NOW_config(N=20, motionCompensation=mc)
        x = np.random.randn(3 * 20 + 1)
        vals = evaluate_all_nonlinear(x, c)
        # 1 + 19 + 3 + 1 + 1 (nonlinear motion, order 2) = 25
        assert len(vals) == 25

    @skip_without_fixture('constraints_ste_default.mat')
    def test_values_against_matlab(self):
        f = load_fixture('constraints_ste_default.mat')
        c = NOW_config()
        vals = evaluate_all_nonlinear(f['x0'], c)
        np.testing.assert_array_almost_equal(vals, f['nonlinear_c'], decimal=8)

    @skip_without_fixture('constraints_ste_nomax.mat')
    def test_values_no_maxwell_against_matlab(self):
        f = load_fixture('constraints_ste_nomax.mat')
        c = NOW_config(doMaxwellComp=False)
        vals = evaluate_all_nonlinear(f['x0'], c)
        np.testing.assert_array_almost_equal(vals, f['nonlinear_c'], decimal=8)

    @skip_without_fixture('constraints_lte.mat')
    def test_values_lte_against_matlab(self):
        f = load_fixture('constraints_lte.mat')
        c = NOW_config(targetTensor=np.diag([1, 0, 0]))
        vals = evaluate_all_nonlinear(f['x0'], c)
        np.testing.assert_array_almost_equal(vals, f['nonlinear_c'], decimal=8)


class TestNonlinearJacobian:
    """Test analytical Jacobians of nonlinear constraints."""

    def test_jacobian_shape_default(self):
        c = NOW_config(N=20)
        x = np.random.randn(3 * 20 + 1)
        jac = evaluate_all_nonlinear_jacobian(x, c)
        vals = evaluate_all_nonlinear(x, c)
        assert jac.shape == (len(vals), 3 * 20 + 1)

    def test_jacobian_vs_finite_differences(self):
        """Compare analytical Jacobian to finite differences."""
        c = NOW_config(N=10, doMaxwellComp=False)
        np.random.seed(123)
        x = np.random.randn(3 * 10 + 1) * 0.1

        jac_analytical = evaluate_all_nonlinear_jacobian(x, c)

        eps = 1e-6
        n = len(x)
        vals0 = evaluate_all_nonlinear(x, c)
        jac_fd = np.zeros((len(vals0), n))
        for i in range(n):
            x_plus = x.copy()
            x_plus[i] += eps
            x_minus = x.copy()
            x_minus[i] -= eps
            jac_fd[:, i] = (evaluate_all_nonlinear(x_plus, c) -
                            evaluate_all_nonlinear(x_minus, c)) / (2 * eps)

        np.testing.assert_allclose(jac_analytical, jac_fd, atol=1e-3, rtol=1e-3)

    def test_jacobian_vs_finite_differences_with_maxwell(self):
        """Compare analytical Jacobian to finite differences with Maxwell."""
        c = NOW_config(N=10)
        np.random.seed(456)
        x = np.random.randn(3 * 10 + 1) * 0.5

        jac_analytical = evaluate_all_nonlinear_jacobian(x, c)

        eps = 1e-6
        n = len(x)
        vals0 = evaluate_all_nonlinear(x, c)
        jac_fd = np.zeros((len(vals0), n))
        for i in range(n):
            x_plus = x.copy()
            x_plus[i] += eps
            x_minus = x.copy()
            x_minus[i] -= eps
            jac_fd[:, i] = (evaluate_all_nonlinear(x_plus, c) -
                            evaluate_all_nonlinear(x_minus, c)) / (2 * eps)

        np.testing.assert_allclose(jac_analytical, jac_fd, atol=1e-3, rtol=1e-3)

    def test_jacobian_vs_jax(self):
        """Compare analytical Jacobian to JAX autodiff."""
        try:
            import jax
            import jax.numpy as jnp
        except ImportError:
            pytest.skip("JAX not available")

        c = NOW_config(N=10, doMaxwellComp=False)
        np.random.seed(789)
        x = np.random.randn(3 * 10 + 1) * 0.1

        # JAX version of constraint evaluation
        @jax.jit
        def jax_constraints(x_jax):
            N = 10
            Q = x_jax[:-1].reshape(N, 3, order='F')
            s = x_jax[-1]

            A1 = jnp.array(np.array([[-1 if i == j else (1 if i + 1 == j else 0)
                                       for j in range(N)] for i in range(N - 1)], dtype=float))
            w = jnp.ones((N, 1))
            w = w.at[0].set(0.5)
            w = w.at[-1].set(0.5)
            w = w / (N - 1)

            g = A1 @ Q
            wQ = w * Q
            B = Q.T @ wQ

            Bt = jnp.eye(3)
            tol = 0.5e-2

            c1 = jnp.trace((B - s * Bt).T @ (B - s * Bt)) - (s * tol) ** 2
            c2 = jnp.sum(g ** 2, axis=1) - c._gMaxConstraint ** 2
            c3 = g[:, 0] @ g[:, 0] - c._integralConstraint
            c4 = g[:, 1] @ g[:, 1] - c._integralConstraint
            c5 = g[:, 2] @ g[:, 2] - c._integralConstraint

            return jnp.concatenate([jnp.atleast_1d(c1), c2,
                                    jnp.atleast_1d(c3), jnp.atleast_1d(c4), jnp.atleast_1d(c5)])

        jax_jac_fn = jax.jacobian(jax_constraints)
        jac_jax = np.array(jax_jac_fn(jnp.array(x)))
        jac_analytical = evaluate_all_nonlinear_jacobian(x, c)

        np.testing.assert_array_almost_equal(jac_analytical, jac_jax, decimal=6)

    @skip_without_fixture('constraints_ste_default.mat')
    def test_jacobian_against_matlab(self):
        f = load_fixture('constraints_ste_default.mat')
        c = NOW_config()
        jac = evaluate_all_nonlinear_jacobian(f['x0'], c)
        # MATLAB gradc is transposed (columns = constraints in fmincon)
        # Our fixture stores gradc already in (3N+1, n_constraints) form
        matlab_jac = f['nonlinear_gradc'].T  # transpose to (n_constraints, 3N+1)
        np.testing.assert_array_almost_equal(jac, matlab_jac, decimal=8)

    def test_jacobian_maxnorm(self):
        """With useMaxNorm=True, gradient norm rows should be absent."""
        c = NOW_config(N=10, useMaxNorm=True)
        np.random.seed(111)
        x = np.random.randn(3 * 10 + 1) * 0.1
        jac = evaluate_all_nonlinear_jacobian(x, c)
        vals = evaluate_all_nonlinear(x, c)
        assert jac.shape == (len(vals), 3 * 10 + 1)
        # 1 (tensor) + 0 (no grad norm) + 3 (power) + 1 (maxwell) = 5
        assert jac.shape[0] == 5

    def test_jacobian_maxnorm_vs_finite_differences(self):
        """Verify maxnorm Jacobian against finite differences."""
        c = NOW_config(N=10, useMaxNorm=True)
        np.random.seed(222)
        x = np.random.randn(3 * 10 + 1) * 0.1
        jac_analytical = evaluate_all_nonlinear_jacobian(x, c)

        eps = 1e-6
        vals0 = evaluate_all_nonlinear(x, c)
        jac_fd = np.zeros((len(vals0), len(x)))
        for i in range(len(x)):
            x_plus = x.copy(); x_plus[i] += eps
            x_minus = x.copy(); x_minus[i] -= eps
            jac_fd[:, i] = (evaluate_all_nonlinear(x_plus, c) -
                            evaluate_all_nonlinear(x_minus, c)) / (2 * eps)

        np.testing.assert_allclose(jac_analytical, jac_fd, atol=1e-3, rtol=1e-3)

    def test_nonlinear_motion_compensation(self):
        """Nonlinear motion compensation constraints and Jacobian."""
        mc = {'order': [2], 'maxMagnitude': [1e-4]}
        c = NOW_config(N=10, motionCompensation=mc)
        np.random.seed(333)
        x = np.random.randn(3 * 10 + 1) * 0.1

        vals = evaluate_all_nonlinear(x, c)
        jac = evaluate_all_nonlinear_jacobian(x, c)
        # 1 (tensor) + 9 (grad norm) + 3 (power) + 1 (maxwell) + 1 (motion) = 15
        assert len(vals) == 15
        assert jac.shape == (15, 3 * 10 + 1)

        # Verify motion Jacobian against finite differences
        eps = 1e-6
        jac_fd = np.zeros_like(jac)
        for i in range(len(x)):
            x_plus = x.copy(); x_plus[i] += eps
            x_minus = x.copy(); x_minus[i] -= eps
            jac_fd[:, i] = (evaluate_all_nonlinear(x_plus, c) -
                            evaluate_all_nonlinear(x_minus, c)) / (2 * eps)

        np.testing.assert_allclose(jac, jac_fd, atol=1e-3, rtol=1e-3)
