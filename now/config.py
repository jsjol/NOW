import numpy as np


def getActualTimings(durationFirstPartRequested,
                     durationZeroGradientRequested,
                     durationSecondPartRequested,
                     discretizationSteps,
                     forceSymmetry):
    """Compute actual timings after discretization.

    Returns (durationFirstPartActual, durationZeroGradientActual,
             durationSecondPartActual, totalTimeActual, zeroGradientAtIndex)

    zeroGradientAtIndex uses 0-based Python indexing (MATLAB indices - 1).
    """
    if forceSymmetry:
        durBothRequested = min(durationFirstPartRequested, durationSecondPartRequested)
        totalTime = 2 * durBothRequested + durationZeroGradientRequested
        dt = totalTime / discretizationSteps

        num_zero = int(np.ceil(durationZeroGradientRequested / dt))
        if (discretizationSteps - num_zero) % 2 != 0:
            num_zero -= 1

        # MATLAB 1-based: startZeroGradientsIndex = (N - num_zero)/2
        # Python 0-based: subtract 1
        startZeroGradientsIndex_matlab = (discretizationSteps - num_zero) // 2
        startSecondPartIndex_matlab = startZeroGradientsIndex_matlab + num_zero

        durationFirstPartActual = startZeroGradientsIndex_matlab * dt
        durationSecondPartActual = (discretizationSteps - startSecondPartIndex_matlab) * dt
        durationZeroGradientActual = (startSecondPartIndex_matlab - startZeroGradientsIndex_matlab) * dt
        totalTimeActual = durationFirstPartActual + durationSecondPartActual + durationZeroGradientActual

        if abs(durationFirstPartActual - durationSecondPartActual) > 1e-10 * totalTimeActual:
            raise ValueError('The first and second halves need to be equally long to enforce symmetry.')

        if durationZeroGradientActual > 0:
            # MATLAB: (startZeroGradientsIndex:startSecondPartIndex) — 1-based inclusive
            # Python: convert to 0-based
            zeroGradientAtIndex = np.arange(
                startZeroGradientsIndex_matlab - 1,
                startSecondPartIndex_matlab,  # -1+1 = 0 offset for inclusive end
                dtype=int)
        else:
            zeroGradientAtIndex = np.array([], dtype=int)
    else:
        totalTime = durationFirstPartRequested + durationSecondPartRequested + durationZeroGradientRequested
        dt = totalTime / discretizationSteps

        # MATLAB uses 1-based indexing
        startZeroGradientsIndex_matlab = int(np.floor(durationFirstPartRequested / dt))
        startSecondPartIndex_matlab = int(np.floor(
            (durationFirstPartRequested + durationZeroGradientRequested) / dt))

        durationFirstPartActual = startZeroGradientsIndex_matlab * dt
        durationSecondPartActual = (discretizationSteps - startSecondPartIndex_matlab) * dt
        durationZeroGradientActual = (startSecondPartIndex_matlab - startZeroGradientsIndex_matlab) * dt
        totalTimeActual = durationFirstPartActual + durationSecondPartActual + durationZeroGradientActual

        if durationZeroGradientActual > 0:
            # MATLAB: (startZeroGradientsIndex:startSecondPartIndex) — 1-based inclusive
            # Python: convert to 0-based
            zeroGradientAtIndex = np.arange(
                startZeroGradientsIndex_matlab - 1,
                startSecondPartIndex_matlab,  # -1+1 = 0 offset for inclusive end
                dtype=int)
        else:
            zeroGradientAtIndex = np.array([], dtype=int)

    return (durationFirstPartActual,
            durationZeroGradientActual,
            durationSecondPartActual,
            totalTimeActual,
            zeroGradientAtIndex)


class NOW_config:
    """Configuration for a NOW optimization problem.

    Parameters
    ----------
    targetTensor : (3, 3) array, optional
        Target b-tensor shape. Default ``np.eye(3)`` (spherical encoding).
        Use ``np.diag([1,0,0])`` for linear or ``np.diag([1,1,0])`` for planar.
    N : int
        Number of time-discretization points (default 77).
    gMax : float
        Maximum gradient amplitude in mT/m (default 80).
    sMax : float
        Maximum slew rate in T/m/s (default 100).
    durationFirstPartRequested : float
        Requested duration before the zero-gradient pause, in ms (default 28).
    durationSecondPartRequested : float
        Requested duration after the zero-gradient pause, in ms (default 22).
    durationZeroGradientRequested : float
        Requested duration of the zero-gradient pause, in ms (default 8).
    eta : float
        Energy/efficacy balance in (0, 1]. Lower values penalise coil heating
        more aggressively (default 1).
    useMaxNorm : bool
        If True, constrain gradient amplitude per axis (max-norm).
        If False (default), constrain L2-norm across axes.
    enforceSymmetry : bool
        Force waveform symmetry about the zero-gradient interval (default False).
    doMaxwellComp : bool
        Enable Maxwell (concomitant gradient) compensation (default True).
    MaxwellIndex : float
        Maxwell term threshold in (mT/m)² ms (default 100).
    motionCompensation : dict or None
        If not None, a dict with keys ``'order'`` (array of int) and
        ``'maxMagnitude'`` (array of float). An order with maxMagnitude=0
        is enforced exactly (linear constraint); maxMagnitude>0 allows that
        deviation (nonlinear constraint).
    doBackgroundCompensation : int
        0 = off, 1 = general timing condition (requires velocity nulling),
        2 = specific timing condition (requires ``startTime``).
    startTime : float
        Time from excitation to first gradient sample in ms. Used when
        ``doBackgroundCompensation=2``.

    Attributes (derived, read-only)
    -------------------------------
    durationFirstPartActual, durationSecondPartActual, durationZeroGradientActual : float
        Actual durations after time-discretization, in ms.
    totalTimeActual : float
        Total actual encoding time in ms.
    zeroGradientAtIndex : ndarray
        0-based indices where the gradient must be zero.
    signs : ndarray or None
        Spin dephasing direction vector, shape (N-1, 1).
    """
    def __init__(
        self,
        targetTensor=None,
        N=77,
        useMaxNorm=False,
        gMax=80,
        sMax=100,
        durationFirstPartRequested=28,
        durationSecondPartRequested=22,
        durationZeroGradientRequested=8,
        eta=1,
        enforceSymmetry=False,
        redoIfFailed=True,
        name='NOW',
        initialGuess='random',
        x0=None,
        doMaxwellComp=True,
        MaxwellIndex=100,
        MaxFunEval=1e5,
        MaxIter=5e3,
        motionCompensation=None,
        doBackgroundCompensation=0,
        startTime=0):

        if targetTensor is None:
            targetTensor = np.eye(3)
        self.targetTensor = np.asarray(targetTensor, dtype=float)
        self.N = N
        self.useMaxNorm = useMaxNorm
        self.gMax = gMax
        self.sMax = sMax
        self.durationFirstPartRequested = durationFirstPartRequested
        self.durationSecondPartRequested = durationSecondPartRequested
        self.durationZeroGradientRequested = durationZeroGradientRequested
        self.eta = eta
        self.enforceSymmetry = enforceSymmetry
        self.redoIfFailed = redoIfFailed
        self.name = name
        self.initialGuess = initialGuess
        self.doMaxwellComp = doMaxwellComp
        self.MaxwellIndex = MaxwellIndex
        self.MaxFunEval = int(MaxFunEval)
        self.MaxIter = int(MaxIter)
        self.doBackgroundCompensation = doBackgroundCompensation
        self.startTime = startTime

        # Motion compensation: default to empty struct matching MATLAB
        if motionCompensation is None:
            self.motionCompensation = {
                'order': np.array([], dtype=float),
                'maxMagnitude': np.array([], dtype=float),
                'linear': np.array([], dtype=bool),
            }
        else:
            self.motionCompensation = {
                'order': np.asarray(motionCompensation['order'], dtype=float),
                'maxMagnitude': np.asarray(motionCompensation['maxMagnitude'], dtype=float),
                'linear': np.array([], dtype=bool),  # will be inferred below
            }

        # Validate motion compensation
        mc = self.motionCompensation
        if len(mc['maxMagnitude']) != len(mc['order']):
            raise ValueError('motionCompensation.maxMagnitude must have the same size as motionCompensation.order.')

        if len(mc['maxMagnitude']) == 0:
            mc['linear'] = np.array([], dtype=bool)
        else:
            mc['linear'] = mc['maxMagnitude'] <= 0

        # Cross-term compensation
        if self.doBackgroundCompensation == 1:
            ind_velo = np.where(mc['order'] == 1)[0]
            if len(ind_velo) == 0:
                idx = len(mc['order'])
                mc['order'] = np.append(mc['order'], 1)
                mc['linear'] = np.append(mc['linear'], True)
                mc['maxMagnitude'] = np.append(mc['maxMagnitude'], 0)
            else:
                if not mc['linear'][ind_velo[0]]:
                    raise ValueError('Cross-term-compensation for a general timing requires velocity compensation!')
        elif self.doBackgroundCompensation == 2:
            if self.startTime == 0:
                import warnings
                warnings.warn('Start time for waveform is t = 0, which is an unlikely setting. Please check!')
            if self.startTime < 0:
                raise ValueError('Start time cannot be smaller than zero!')
        elif self.doBackgroundCompensation != 0:
            raise ValueError('Selection for Cross-term-compensation not recognized! Use value 0, 1 or 2.')

        # Initial guess
        if x0 is not None:
            self.x0 = np.asarray(x0, dtype=float)
        else:
            self.x0 = None

        # Compute actual timings
        (self.durationFirstPartActual,
         self.durationZeroGradientActual,
         self.durationSecondPartActual,
         self.totalTimeActual,
         self.zeroGradientAtIndex) = getActualTimings(
            self.durationFirstPartRequested,
            self.durationZeroGradientRequested,
            self.durationSecondPartRequested,
            self.N,
            self.enforceSymmetry)

        # Derived constraints
        self._dt = self.totalTimeActual / self.N
        self._gMaxConstraint = self.gMax * self._dt
        self._sMaxConstraint = self.sMax * self._dt ** 2
        self._integralConstraint = self.eta * self._gMaxConstraint ** 2 * self.totalTimeActual / self._dt
        self._tolIsotropy = 0.5e-2

        if self.doMaxwellComp:
            self._tolMaxwell = self.MaxwellIndex / self._dt
        else:
            self._tolMaxwell = np.inf

        # Create spin dephasing direction vector
        self.signs = None
        if len(self.zeroGradientAtIndex) > 0:
            signs = np.ones((self.N - 1, 1))
            # MATLAB: mi = median(zeroGradientAtIndex) — 1-based
            # Python: our indices are 0-based, but median should still give
            # the midpoint. MATLAB signs(round(mi):end) = -1 uses 1-based index.
            # We need to match the same physical position.
            # MATLAB zeroGradientAtIndex is 1-based; ours is 0-based.
            # MATLAB mi = median(matlab_indices). Python mi = median(python_indices) + 1 = median(matlab_indices).
            # So we need: signs[round(mi_matlab) - 1:] = -1 in Python.
            mi_0based = np.median(self.zeroGradientAtIndex)
            mi_matlab = mi_0based + 1  # convert to MATLAB 1-based
            idx = int(np.round(mi_matlab)) - 1  # back to 0-based for array indexing
            signs[idx:] = -1
            if mi_matlab == np.round(mi_matlab):
                signs[idx] = 0
            self.signs = signs
