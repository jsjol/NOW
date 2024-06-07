# %%
from warnings import warn
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass


# %%
@dataclass
class motionCompensation:
    order: int
    maxMagnitude: float
    linear: bool

def getActualTimings(durationFirstPartRequested,
                     durationZeroGradientRequested,
                     durationSecondPartRequested,
                     discretizationSteps,
                     forceSymmetry):
  if forceSymmetry:
    raise NotImplementedError
  else:
    totalTime = durationFirstPartRequested + durationSecondPartRequested + durationZeroGradientRequested
    dt = totalTime/discretizationSteps
    
    startZeroGradientsIndex  = np.floor(durationFirstPartRequested/dt) - 1
    startSecondPartIndex = np.floor((durationFirstPartRequested+durationZeroGradientRequested)/dt) - 1 
    
    durationFirstPartActual = (startZeroGradientsIndex + 1) * dt
    durationSecondPartActual = (discretizationSteps - startSecondPartIndex - 1) * dt
    durationZeroGradientActual  = (startSecondPartIndex - startZeroGradientsIndex) * dt
    totalTimeActual = durationFirstPartActual + durationSecondPartActual + durationZeroGradientActual
    
    if durationZeroGradientActual > 0:
        zeroGradientAtIndex = np.arange(startZeroGradientsIndex, startSecondPartIndex + 1)
    else:
        zeroGradientAtIndex = None
        
  return (durationFirstPartActual, 
          durationZeroGradientActual,
          durationSecondPartActual,
          totalTimeActual,
          zeroGradientAtIndex)

class NOW:
  def __init__(
    self,
    targetTensor = np.eye(3), # Isotropic encoding tensor,
    N = 77,
    initialGuess = 'random',
    useMaxNorm = False,
    gMax = 80,
    sMax = 100,
    durationFirstPartRequested = 28,
    durationSecondPartRequested = 22,
    durationZeroGradientRequested = 8,
    eta = 1,
    enforceSymmetry = False,
    redoIfFailed = True,
    name = 'NOW',
    x0 = None,
    doMaxwellComp = True,
    MaxwellIndex = 100,
    MaxFunEval = 1e5,
    MaxIter    = 5e3,
    motionCompensation = None,
    doBackgroundCompensation = 0, # 0 = off; 1 = general timing cond.; 2 = specific timing cond.
    startTime = 0): # Time from the excitation (t=0) to the first gradient waveform sample in ms.
    
    # Parse kwargs
    self.targetTensor = targetTensor
    self.N = N
    self.initialGuess = initialGuess
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
    self.x0 = x0
    self.doMaxwellComp = doMaxwellComp
    self.MaxwellIndex = MaxwellIndex
    self.MaxFunEval = MaxFunEval
    self.MaxIter = MaxIter
    self.motionCompensation = motionCompensation
    self.doBackgroundCompensation = doBackgroundCompensation
    self.startTime = startTime

    # Get actual times after discretization
    self.durationFirstPartActual, self.durationZeroGradientActual, self.durationSecondPartActual, self.totalTimeActual, self.zeroGradientAtIndex = \
      getActualTimings(self.durationFirstPartRequested, self.durationZeroGradientRequested, self.durationSecondPartRequested, self.N, self.enforceSymmetry)
    
    self._dt = self.totalTimeActual/self.N # Time step in milliseconds. Division by N instead of N-1 due to half step shift in gradients.
    self._gMaxConstraint = self.gMax * self._dt
    self._sMaxConstraint = self.sMax * self._dt**2
    self._integralConstraint = self.eta * self._gMaxConstraint**2 * self.totalTimeActual / self._dt
    
    if self.doMaxwellComp:
      self._tolMaxwell = self.MaxwellIndex/self._dt
    else:
      self._tolMaxwell = np.inf
    

# %%
problem = NOW()
print('------------ Requested timing parameters: ------------ \n')
print('DurPre = {:.3f} DurPost = {:.3f}  DurPi = {:.3f}  [ms]\n\n'.format(problem.durationFirstPartRequested, problem.durationSecondPartRequested, problem.durationZeroGradientRequested))
print('------------   Actual timing parameters:  ------------ \n')
print('DurPre = {:.3f} DurPost = {:.3f}  DurPi = {:.3f}  [ms]\n\n'.format(problem.durationFirstPartActual, problem.durationSecondPartActual, problem.durationZeroGradientActual))

# %%
