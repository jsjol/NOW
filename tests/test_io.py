"""Tests for config and problem I/O."""
import json
import numpy as np
import pytest
from scipy.io import loadmat
from now.config import NOW_config
from now.problem import build_problem
from now.io.config_io import save_config, load_config
from now.io.problem_io import export_problem


class TestConfigDict:
    """Test NOW_config.to_dict / from_dict round-trip."""

    def test_default_roundtrip(self):
        c = NOW_config()
        d = c.to_dict()
        c2 = NOW_config.from_dict(d)
        assert c2.N == c.N
        assert c2.gMax == c.gMax
        assert c2.sMax == c.sMax
        assert c2.eta == c.eta
        np.testing.assert_array_equal(c2.targetTensor, c.targetTensor)

    def test_derived_values_match(self):
        c = NOW_config(N=50, gMax=60, eta=0.8)
        c2 = NOW_config.from_dict(c.to_dict())
        assert c2._dt == pytest.approx(c._dt)
        assert c2._gMaxConstraint == pytest.approx(c._gMaxConstraint)
        assert c2._sMaxConstraint == pytest.approx(c._sMaxConstraint)
        assert c2._integralConstraint == pytest.approx(c._integralConstraint)
        assert c2._tolIsotropy == pytest.approx(c._tolIsotropy)
        assert c2._tolMaxwell == pytest.approx(c._tolMaxwell)
        assert c2.totalTimeActual == pytest.approx(c.totalTimeActual)
        np.testing.assert_array_equal(c2.zeroGradientAtIndex, c.zeroGradientAtIndex)

    def test_custom_target_tensor(self):
        T = np.diag([1.0, 1.0, 0.0])
        c = NOW_config(targetTensor=T)
        c2 = NOW_config.from_dict(c.to_dict())
        np.testing.assert_array_equal(c2.targetTensor, T)

    def test_motion_compensation_roundtrip(self):
        mc = {'order': [1, 2], 'maxMagnitude': [0, 1e-4]}
        c = NOW_config(motionCompensation=mc)
        c2 = NOW_config.from_dict(c.to_dict())
        np.testing.assert_array_almost_equal(c2.motionCompensation['order'], [1, 2])
        np.testing.assert_array_almost_equal(c2.motionCompensation['maxMagnitude'], [0, 1e-4])

    def test_no_motion_compensation(self):
        c = NOW_config()
        d = c.to_dict()
        assert d['motionCompensation'] is None
        c2 = NOW_config.from_dict(d)
        assert len(c2.motionCompensation['order']) == 0

    def test_symmetry_roundtrip(self):
        c = NOW_config(enforceSymmetry=True,
                       durationFirstPartRequested=25,
                       durationSecondPartRequested=25,
                       durationZeroGradientRequested=8)
        c2 = NOW_config.from_dict(c.to_dict())
        assert c2.enforceSymmetry is True
        np.testing.assert_array_equal(c2.zeroGradientAtIndex, c.zeroGradientAtIndex)

    def test_x0_roundtrip(self):
        x0 = np.random.randn(3 * 20 + 1)
        c = NOW_config(N=20, x0=x0, initialGuess='user-provided')
        c2 = NOW_config.from_dict(c.to_dict())
        np.testing.assert_array_almost_equal(c2.x0, x0)
        assert c2.initialGuess == 'user-provided'

    def test_to_dict_is_json_serializable(self):
        c = NOW_config()
        d = c.to_dict()
        s = json.dumps(d)
        assert isinstance(s, str)


class TestConfigFileIO:
    """Test save_config / load_config."""

    def test_json_roundtrip(self, tmp_path):
        c = NOW_config(N=30, gMax=70, eta=0.9)
        path = tmp_path / 'config.json'
        save_config(c, path)
        c2 = load_config(path)
        assert c2.N == 30
        assert c2.gMax == 70
        assert c2.eta == pytest.approx(0.9)

    def test_json_with_motion_comp(self, tmp_path):
        mc = {'order': [1, 2], 'maxMagnitude': [0, 1e-4]}
        c = NOW_config(motionCompensation=mc)
        path = tmp_path / 'config.json'
        save_config(c, path)
        c2 = load_config(path)
        np.testing.assert_array_almost_equal(c2.motionCompensation['order'], [1, 2])

    def test_yaml_roundtrip(self, tmp_path):
        pytest.importorskip('yaml')
        c = NOW_config(N=40, gMax=65)
        path = tmp_path / 'config.yaml'
        save_config(c, path)
        c2 = load_config(path)
        assert c2.N == 40
        assert c2.gMax == 65

    def test_yaml_import_error(self, tmp_path, monkeypatch):
        import builtins
        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == 'yaml':
                raise ImportError("no yaml")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, '__import__', mock_import)
        c = NOW_config()
        with pytest.raises(ImportError, match="pyyaml"):
            save_config(c, tmp_path / 'config.yaml')


class TestProblemExport:
    """Test export_problem to npz and mat formats."""

    def test_npz_roundtrip(self, tmp_path):
        c = NOW_config(N=15)
        problem = build_problem(c)
        path = tmp_path / 'problem.npz'
        export_problem(problem, path, fmt='npz')
        data = np.load(str(path))
        np.testing.assert_array_equal(data['A_ineq'], problem.linear.A_ineq)
        np.testing.assert_array_equal(data['b_ineq'], problem.linear.b_ineq)
        np.testing.assert_array_equal(data['A_eq'], problem.linear.A_eq)
        np.testing.assert_array_equal(data['b_eq'], problem.linear.b_eq)
        assert int(data['n_vars']) == problem.n_vars

    def test_mat_loadable(self, tmp_path):
        c = NOW_config(N=15)
        problem = build_problem(c)
        path = tmp_path / 'problem.mat'
        export_problem(problem, path, fmt='mat')
        data = loadmat(str(path))
        np.testing.assert_array_almost_equal(
            data['A_ineq'], problem.linear.A_ineq)
        np.testing.assert_array_almost_equal(
            data['A_eq'], problem.linear.A_eq)

    def test_nonlinear_sample_included(self, tmp_path):
        c = NOW_config(N=15)
        problem = build_problem(c)
        path = tmp_path / 'problem.npz'
        export_problem(problem, path)
        data = np.load(str(path))
        x0 = np.zeros(problem.n_vars)
        np.testing.assert_array_almost_equal(
            data['nonlinear_at_zero'], problem.nonlinear.fun(x0))

    def test_unsupported_format(self, tmp_path):
        c = NOW_config(N=15)
        problem = build_problem(c)
        with pytest.raises(ValueError, match="Unsupported format"):
            export_problem(problem, tmp_path / 'problem.hdf5', fmt='hdf5')
