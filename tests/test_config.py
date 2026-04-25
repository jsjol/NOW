"""Tests for NOW_config and getActualTimings."""
import numpy as np
import pytest
from now.config import NOW_config, getActualTimings
from .helpers import load_fixture, skip_without_fixture


class TestGetActualTimings:
    """Test getActualTimings against known values and MATLAB fixtures."""

    def test_default_timings_positive(self):
        dur1, dur0, dur2, total, zgi = getActualTimings(28, 8, 22, 77, False)
        assert dur1 > 0
        assert dur0 > 0
        assert dur2 > 0
        assert abs(dur1 + dur0 + dur2 - total) < 1e-10

    def test_total_time_consistency(self):
        dur1, dur0, dur2, total, zgi = getActualTimings(28, 8, 22, 77, False)
        assert abs(total - (28 + 8 + 22)) < 2  # within 2ms of requested

    def test_zero_gradient_indices_exist(self):
        dur1, dur0, dur2, total, zgi = getActualTimings(28, 8, 22, 77, False)
        assert len(zgi) > 0

    def test_zero_gradient_no_pause(self):
        dur1, dur0, dur2, total, zgi = getActualTimings(30, 0, 28, 77, False)
        assert len(zgi) == 0

    def test_symmetry_equal_durations(self):
        dur1, dur0, dur2, total, zgi = getActualTimings(25, 8, 25, 77, True)
        assert abs(dur1 - dur2) < 1e-10

    def test_symmetry_unequal_requested(self):
        """When forceSymmetry=True, min(dur1, dur2) is used."""
        dur1, dur0, dur2, total, zgi = getActualTimings(28, 8, 22, 77, True)
        assert abs(dur1 - dur2) < 1e-10

    @skip_without_fixture('config_ste_default.mat')
    def test_default_against_matlab(self):
        f = load_fixture('config_ste_default.mat')
        dur1, dur0, dur2, total, zgi = getActualTimings(28, 8, 22, 77, False)
        assert abs(dur1 - f['durationFirstPartActual']) < 1e-10
        assert abs(dur0 - f['durationZeroGradientActual']) < 1e-10
        assert abs(dur2 - f['durationSecondPartActual']) < 1e-10
        assert abs(total - f['totalTimeActual']) < 1e-10
        # Compare indices (MATLAB 1-based → Python 0-based)
        matlab_zgi = f['zeroGradientAtIndex']
        if matlab_zgi.ndim == 0:
            matlab_zgi = np.array([matlab_zgi.item()])
        np.testing.assert_array_equal(zgi, matlab_zgi - 1)

    @skip_without_fixture('config_ste_symmetric.mat')
    def test_symmetric_against_matlab(self):
        f = load_fixture('config_ste_symmetric.mat')
        dur1, dur0, dur2, total, zgi = getActualTimings(25, 8, 25, 77, True)
        assert abs(dur1 - f['durationFirstPartActual']) < 1e-10
        assert abs(dur0 - f['durationZeroGradientActual']) < 1e-10
        assert abs(dur2 - f['durationSecondPartActual']) < 1e-10


class TestNOWConfig:
    """Test NOW_config initialization and derived quantities."""

    def test_default_creation(self):
        c = NOW_config()
        assert c.N == 77
        assert c.gMax == 80
        assert c.sMax == 100
        np.testing.assert_array_equal(c.targetTensor, np.eye(3))

    def test_derived_quantities(self):
        c = NOW_config()
        assert c._dt > 0
        assert c._gMaxConstraint > 0
        assert c._sMaxConstraint > 0
        assert c._integralConstraint > 0
        assert c._tolIsotropy == 0.5e-2

    def test_maxwell_off(self):
        c = NOW_config(doMaxwellComp=False)
        assert np.isinf(c._tolMaxwell)

    def test_maxwell_on(self):
        c = NOW_config(doMaxwellComp=True)
        assert not np.isinf(c._tolMaxwell)
        assert c._tolMaxwell > 0

    def test_signs_vector(self):
        c = NOW_config()
        assert c.signs is not None
        assert c.signs.shape == (c.N - 1, 1)
        # Should have +1, 0, -1 pattern
        assert np.any(c.signs > 0)
        assert np.any(c.signs < 0)

    def test_no_zero_gradient_no_signs(self):
        c = NOW_config(durationZeroGradientRequested=0)
        assert c.signs is None

    def test_motion_compensation_validation(self):
        mc = {'order': np.array([1, 2]), 'maxMagnitude': np.array([0, 1e-4])}
        c = NOW_config(motionCompensation=mc)
        assert c.motionCompensation['linear'][0] == True  # maxMag=0 → linear
        assert c.motionCompensation['linear'][1] == False  # maxMag>0 → nonlinear

    def test_background_comp_adds_velocity(self):
        """doBackgroundCompensation=1 should force velocity compensation if not present."""
        c = NOW_config(doBackgroundCompensation=1)
        mc = c.motionCompensation
        assert 1.0 in mc['order']
        idx = np.where(mc['order'] == 1)[0][0]
        assert mc['linear'][idx] == True

    def test_lte_config(self):
        c = NOW_config(targetTensor=np.diag([1, 0, 0]))
        np.testing.assert_array_equal(c.targetTensor, np.diag([1, 0, 0]))

    @skip_without_fixture('config_ste_default.mat')
    def test_derived_against_matlab(self):
        f = load_fixture('config_ste_default.mat')
        c = NOW_config()
        assert abs(c._dt - f['dt']) < 1e-10
        assert abs(c._gMaxConstraint - f['gMaxConstraint']) < 1e-10
        assert abs(c._sMaxConstraint - f['sMaxConstraint']) < 1e-10
        assert abs(c._integralConstraint - f['integralConstraint']) < 1e-10
        assert abs(c._tolMaxwell - f['tolMaxwell']) < 1e-10

    @skip_without_fixture('config_ste_default.mat')
    def test_signs_against_matlab(self):
        f = load_fixture('config_ste_default.mat')
        c = NOW_config()
        matlab_signs = f['signs'].reshape(-1, 1)
        np.testing.assert_array_almost_equal(c.signs, matlab_signs)

    @skip_without_fixture('config_ste_motion.mat')
    def test_motion_comp_against_matlab(self):
        f = load_fixture('config_ste_motion.mat')
        mc = {'order': [1, 2], 'maxMagnitude': [0, 1e-4]}
        c = NOW_config(motionCompensation=mc)
        assert abs(c._dt - f['dt']) < 1e-10
