import numpy as np

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
        zeroGradientAtIndex = np.arange(startZeroGradientsIndex, startSecondPartIndex + 1, dtype=int)
    else:
        zeroGradientAtIndex = None 
        
  return (durationFirstPartActual, 
          durationZeroGradientActual,
          durationSecondPartActual,
          totalTimeActual,
          zeroGradientAtIndex)

def obj_fun(x):
    f = -x[-1] #minimize negative of b
    g = np.zeros_like(x)
    g[-1] = -1.0
    return f, g


class NOW_config:
  def __init__(
    self,
    targetTensor = np.eye(3), # Isotropic encoding tensor,
    N = 50, #77
    useMaxNorm = False,
    gMax = 80, #OBS OBS This is in mT/m
    sMax = 100, #OBS OBS This is in T/m/s
    durationFirstPartRequested = 28,
    durationSecondPartRequested = 22,
    durationZeroGradientRequested = 8,
    eta = 1, #1
    enforceSymmetry = False,
    name = 'NOW',
    x0 = None,
    doMaxwellComp = True,
    MaxwellIndex = 100,
    MaxFunEval = 1e5,
    MaxIter    = 5e3,
    motionCompensation = {'order': [], 'maxMagnitude': [], 'linear': []},
    doBackgroundCompensation = 0, # 0 = off; 1 = general timing cond.; 2 = specific timing cond.
    startTime = 0): # Time from the excitation (t=0) to the first gradient waveform sample in ms.
    
    # Parse kwargs
    self.targetTensor = targetTensor #; self.targetTensor[-1,-1] = 0 ; self.targetTensor[1,1] = 0
    self.N = N
    self.useMaxNorm = useMaxNorm
    self.gMax = gMax
    self.sMax = sMax
    self.durationFirstPartRequested = durationFirstPartRequested
    self.durationSecondPartRequested = durationSecondPartRequested
    self.durationZeroGradientRequested = durationZeroGradientRequested
    self.eta = eta
    self.enforceSymmetry = enforceSymmetry
    self.name = name
    self.doMaxwellComp = doMaxwellComp
    self.MaxwellIndex = MaxwellIndex
    self.MaxFunEval = MaxFunEval
    self.MaxIter = MaxIter
    self.motionCompensation = motionCompensation
    self.doBackgroundCompensation = doBackgroundCompensation
    self.startTime = startTime

    
    if x0 is None:
      x0 = np.random.randn(3 * self.N + 1)
    self.x0 = x0

    # Get actual times after discretization
    self.durationFirstPartActual, self.durationZeroGradientActual, self.durationSecondPartActual, self.totalTimeActual, self.zeroGradientAtIndex = \
      getActualTimings(self.durationFirstPartRequested, self.durationZeroGradientRequested, self.durationSecondPartRequested, self.N, self.enforceSymmetry)
    
    self._dt = self.totalTimeActual/self.N # Time step in milliseconds. Division by N instead of N-1 due to half step shift in gradients.
    self._gMaxConstraint = self.gMax * self._dt
    self._sMaxConstraint = self.sMax * self._dt**2
    self._integralConstraint = self.eta * self._gMaxConstraint**2 * self.totalTimeActual / self._dt
    self._tolIsotropy = .5e-2
    
    if self.doMaxwellComp:
      self._tolMaxwell = self.MaxwellIndex/self._dt
    else:
      self._tolMaxwell = np.inf
    
    # Create spin dephasing direction vector
    if self.zeroGradientAtIndex is not None:
      signs = np.ones((self.N - 1, 1)) # Ghost points excluded during opt
      
      # Assume that sign change happens in the middle of the pause
      mi = np.median(self.zeroGradientAtIndex).astype(int)
      signs[np.round(mi):] = -1
      if mi == round(mi):
        signs[mi] = 0
      self.signs = signs
    
    if self.motionCompensation['order'] != []:
      self.motionCompensation['linear'] = (np.array(self.motionCompensation['maxMagnitude']) <= 0)
    else:
      self.motionCompensation['linear'] = []
         

    if self.doBackgroundCompensation != 0:
      raise NotImplementedError
    

if __name__ == "__main__":
  prob = NOW_config(N=11,
                    motionCompensation={'order': [1, 2], 'maxMagnitude': [0, 1e-4], 'linear': []})
  print(prob.motionCompensation['linear'])


  
