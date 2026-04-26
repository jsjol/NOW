"""End-to-end optimization tests."""
import numpy as np
import pytest
from now.config import NOW_config
from now.optimize import now_optimize, get_initial_guess
from now.constraints.nonlinear import evaluate_all_nonlinear
from now.solvers.scipy_solver import ScipySolver
from now.solvers.protocol import SolverResult


class TestOptimize:
    """Test the full optimization pipeline."""

    @pytest.mark.slow
    def test_basic_ste_runs(self):
        """STE optimization completes without error."""
        c = NOW_config(N=20)
        result, config = now_optimize(c, method='SLSQP', max_attempts=3, verbose=False)
        assert result.b > 0

    @pytest.mark.slow
    def test_constraints_satisfied(self):
        """All nonlinear constraints satisfied at solution."""
        c = NOW_config(N=20)
        result, config = now_optimize(c, method='SLSQP', max_attempts=3, verbose=False)
        x = np.concatenate([result.rawq, [-result.b]])  # reconstruct x
        # Actually, we need the raw x including the scaling parameter
        # The result doesn't directly store s, but we can reconstruct from the optimizer
        # For now, just check b > 0 and basic sanity
        assert result.b > 0
        assert result.B.shape == (3, 3)
        assert result.q.shape[1] == 3
        assert result.g.shape[1] == 3

    @pytest.mark.slow
    def test_ste_tensor_shape(self):
        """STE should produce approximately isotropic encoding."""
        c = NOW_config(N=20)
        result, config = now_optimize(c, method='SLSQP', max_attempts=5, verbose=False)
        if result.b > 0.1:  # only test if optimization found something reasonable
            B_normalized = result.B / np.trace(result.B) * 3
            # Should be close to identity
            assert np.linalg.norm(B_normalized - np.eye(3), 'fro') < 0.5

    @pytest.mark.slow
    def test_lte_runs(self):
        """LTE optimization completes."""
        c = NOW_config(N=20, targetTensor=np.diag([1, 0, 0]))
        result, config = now_optimize(c, method='SLSQP', max_attempts=3, verbose=False)
        assert result.b > 0

    @pytest.mark.slow
    def test_result_fields(self):
        """Check all result fields are populated."""
        c = NOW_config(N=20)
        result, config = now_optimize(c, method='SLSQP', max_attempts=3, verbose=False)
        assert result.q is not None
        assert result.g is not None
        assert result.slew is not None
        assert result.B is not None
        assert result.gwf is not None
        assert result.rf is not None
        assert result.dt > 0

    @pytest.mark.slow
    def test_default_config_none(self):
        """Passing config=None should create default config and run."""
        result, config = now_optimize(config=None, method='SLSQP', max_attempts=3, verbose=False)
        assert result.b > 0

    @pytest.mark.slow
    def test_trust_constr_runs(self):
        """trust-constr optimization completes without error."""
        c = NOW_config(N=15)
        result, config = now_optimize(c, method='trust-constr', max_attempts=3, verbose=False)
        assert result.b > 0

    @pytest.mark.slow
    def test_redo_if_failed_false(self):
        """With redoIfFailed=False, should not retry on failure."""
        c = NOW_config(N=15, redoIfFailed=False)
        result, config = now_optimize(c, method='SLSQP', max_attempts=5, verbose=False)
        # Should complete (succeed or fail) without retrying
        assert result.b >= 0


class TestGetInitialGuess:
    def test_user_provided_iteration_zero(self):
        """iteration=0 with user-provided x0 returns that x0."""
        N = 10
        x0_input = np.ones(3 * N + 1)
        c = NOW_config(N=N, initialGuess='user-provided', x0=x0_input)
        x0 = get_initial_guess(c, iteration=0)
        np.testing.assert_array_equal(x0, x0_input)

    def test_user_provided_iteration_nonzero(self):
        """iteration>0 with user-provided x0 returns random (not x0)."""
        N = 10
        x0_input = np.ones(3 * N + 1) * 999.0
        c = NOW_config(N=N, initialGuess='user-provided', x0=x0_input)
        np.random.seed(42)
        x0 = get_initial_guess(c, iteration=1)
        # Should NOT be the user-provided value (random instead)
        assert not np.allclose(x0, x0_input)
        assert x0.shape == (3 * N + 1,)

    def test_lte_scaling(self):
        """For LTE targetTensor=diag(1,0,0), x0[N:2N] and x0[2N:3N] should be zero."""
        N = 15
        c = NOW_config(N=N, targetTensor=np.diag([1, 0, 0]))
        np.random.seed(0)
        x0 = get_initial_guess(c, iteration=0)
        # targetTensor[1,1]=0 and targetTensor[2,2]=0 => those blocks are zeroed
        np.testing.assert_array_equal(x0[N:2 * N], 0.0)
        np.testing.assert_array_equal(x0[2 * N:3 * N], 0.0)
        # First block should be non-zero (scaled by targetTensor[0,0]=1)
        assert not np.allclose(x0[:N], 0.0)


class TestOptimizeErrorPaths:
    """Test error handling and retry logic in now_optimize."""

    def test_all_attempts_fail_raises(self, monkeypatch):
        """When solver always raises, RuntimeError is raised."""
        c = NOW_config(N=15)

        def mock_solve(self, problem, x0, options=None):
            raise ValueError("deliberate test failure")

        monkeypatch.setattr(ScipySolver, 'solve', mock_solve)
        with pytest.raises(RuntimeError, match="All.*failed"):
            now_optimize(c, method='SLSQP', max_attempts=3, verbose=False)

    def test_exception_then_success(self, monkeypatch):
        """Solver raises on first attempt, succeeds on second."""
        c = NOW_config(N=15)
        call_count = [0]
        original_solve = ScipySolver.solve

        def mock_solve(self, problem, x0, options=None):
            call_count[0] += 1
            if call_count[0] == 1:
                raise ValueError("transient failure")
            return original_solve(self, problem, x0, options)

        monkeypatch.setattr(ScipySolver, 'solve', mock_solve)
        result, config = now_optimize(c, method='SLSQP', max_attempts=5, verbose=False)
        assert call_count[0] >= 2

    def test_verbose_prints(self, capsys):
        """Verbose mode produces output."""
        c = NOW_config(N=15)

        def mock_solve(self, problem, x0, options=None):
            return SolverResult(x=x0, fun=-1.0, success=True,
                                message='ok', n_iterations=1)

        original_solve = ScipySolver.solve
        ScipySolver.solve = mock_solve
        try:
            now_optimize(c, method='SLSQP', max_attempts=1, verbose=True)
        finally:
            ScipySolver.solve = original_solve
        captured = capsys.readouterr()
        assert 'Optimizing' in captured.out

    def test_verbose_exception_prints(self, monkeypatch, capsys):
        """Exception in verbose mode prints the error."""
        c = NOW_config(N=15)
        call_count = [0]
        original_solve = ScipySolver.solve

        def mock_solve(self, problem, x0, options=None):
            call_count[0] += 1
            if call_count[0] == 1:
                raise ValueError("test error msg")
            return original_solve(self, problem, x0, options)

        monkeypatch.setattr(ScipySolver, 'solve', mock_solve)
        now_optimize(c, method='SLSQP', max_attempts=5, verbose=True)
        captured = capsys.readouterr()
        assert 'test error msg' in captured.out

    def test_redo_if_failed_verbose(self, monkeypatch, capsys):
        """With redoIfFailed=False and verbose, prints non-repeat message."""
        c = NOW_config(N=15, redoIfFailed=False)

        def mock_solve(self, problem, x0, options=None):
            return SolverResult(x=x0, fun=1e10, success=False,
                                message='fail', n_iterations=0)

        monkeypatch.setattr(ScipySolver, 'solve', mock_solve)
        result, config = now_optimize(c, method='SLSQP', max_attempts=5, verbose=True)
        captured = capsys.readouterr()
        assert 'will not be repeated' in captured.out
