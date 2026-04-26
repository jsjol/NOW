"""Tests for utility functions."""
import numpy as np
from now.utils import (
    gamma, build_first_derivative_matrix, build_second_derivative_matrix,
    build_integration_weights, gwf_to_q, problem_to_name, maxwell_coeff, write_wf
)
from now.config import NOW_config


class TestGamma:
    def test_value(self):
        assert abs(gamma() - 42.576e6) < 1e-3

    def test_gamma_rad(self):
        from now.constants import GAMMA_HZ, GAMMA_RAD
        assert abs(GAMMA_RAD - GAMMA_HZ * 2 * np.pi) < 1
        assert abs(GAMMA_RAD - 2.6751e8) < 1e4


class TestIntegrationWeights:
    def test_shape(self):
        w = build_integration_weights(10)
        assert w.shape == (10, 1)

    def test_trapezoidal(self):
        w = build_integration_weights(5)
        assert w[0, 0] == 0.5 / 4
        assert w[-1, 0] == 0.5 / 4
        assert w[1, 0] == 1.0 / 4

    def test_sum(self):
        w = build_integration_weights(100)
        assert abs(w.sum() - 1.0) < 1e-10


class TestGwfToQ:
    def test_constant_gradient(self):
        gwf = np.ones((10, 3))
        dt = 0.5
        q = gwf_to_q(gwf, dt)
        assert q.shape == (10, 3)
        np.testing.assert_array_almost_equal(q[-1], [5.0, 5.0, 5.0])


class TestProblemToName:
    def test_ste(self):
        c = NOW_config()
        name = problem_to_name(c)
        assert 'STE' in name

    def test_lte(self):
        c = NOW_config(targetTensor=np.diag([1, 0, 0]))
        name = problem_to_name(c)
        assert 'LTE' in name

    def test_pte(self):
        c = NOW_config(targetTensor=np.diag([1, 1, 0]))
        name = problem_to_name(c)
        assert 'PTE' in name

    def test_custom_tensor(self):
        c = NOW_config(targetTensor=np.diag([2, 1, 0]))
        name = problem_to_name(c)
        assert 'custom' in name


class TestMaxwellCoeff:
    def test_identity_gwf(self):
        """Constant gwf with all-positive rf gives known M."""
        T = 5
        gwf = np.ones((T, 3))
        rf = np.ones(T)
        dt = 0.1
        M, mi = maxwell_coeff(gwf, rf, dt)
        # M = gwf.T @ (gwf * rf[:,None]) * dt = ones(3,T) @ ones(T,3) * 0.1 = T * 0.1 * ones(3,3)
        expected_M = T * dt * np.ones((3, 3))
        np.testing.assert_array_almost_equal(M, expected_M)
        expected_mi = np.sqrt(np.trace(expected_M.T @ expected_M))
        assert abs(mi - expected_mi) < 1e-10

    def test_shape(self):
        gwf = np.random.randn(8, 3)
        rf = np.ones(8)
        M, mi = maxwell_coeff(gwf, rf, 0.5)
        assert M.shape == (3, 3)
        assert isinstance(mi, (float, np.floating))

    def test_alternating_rf(self):
        """Alternating rf signs should change the result."""
        gwf = np.ones((4, 3))
        rf_pos = np.ones(4)
        rf_alt = np.array([1, -1, 1, -1], dtype=float)
        dt = 1.0
        M_pos, _ = maxwell_coeff(gwf, rf_pos, dt)
        M_alt, _ = maxwell_coeff(gwf, rf_alt, dt)
        # With alternating signs, signed_gwf columns alternate sign,
        # so M_alt = gwf.T @ signed_gwf * dt where signed_gwf rows alternate sign
        # sum of rf = 0 so off-diagonal / diagonal should be 0
        expected_alt = np.eye(3) * 0.0  # (1+(-1)+1+(-1)) * dt = 0
        np.testing.assert_array_almost_equal(M_alt, expected_alt)


class TestWriteWf:
    def test_single_file_no_sign_change(self, tmp_path):
        """When rf has no sign changes, writes single .txt file."""
        class MockResult:
            gwf = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
            rf = np.array([[1.0], [1.0]])
            dt = 0.001

        result = MockResult()
        fname = str(tmp_path / "test_wf")
        write_wf(result, fname)

        outfile = tmp_path / "test_wf.txt"
        assert outfile.exists()
        data = np.loadtxt(str(outfile))
        np.testing.assert_array_almost_equal(data, result.gwf)

    def test_split_files_with_sign_change(self, tmp_path):
        """When rf has sign changes, writes _A.txt and _B.txt."""
        class MockResult:
            gwf = np.array([[1.0, 0.0, 0.0],
                            [2.0, 0.0, 0.0],
                            [3.0, 0.0, 0.0],
                            [4.0, 0.0, 0.0]])
            rf = np.array([[1.0], [1.0], [-1.0], [-1.0]])
            dt = 0.001

        result = MockResult()
        fname = str(tmp_path / "test_wf")
        write_wf(result, fname)

        file_a = tmp_path / "test_wf_A.txt"
        file_b = tmp_path / "test_wf_B.txt"
        assert file_a.exists()
        assert file_b.exists()

        data_a = np.loadtxt(str(file_a))
        data_b = np.loadtxt(str(file_b))
        # Split at median sign change index + 1
        assert data_a.shape[0] + data_b.shape[0] == result.gwf.shape[0]
