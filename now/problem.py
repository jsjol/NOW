# %% 
import numpy as np
from scipy.optimize import LinearConstraint, NonlinearConstraint
from scipy import sparse
from config import NOW_config

def  defineLinearInequalityConstraints(config):

  constraints = []

  # Constraints on the gradients are linear if we use max-norm
  if config.useMaxNorm: 
      A1 = np.diff(np.eye(config.N)).T
      A1 = np.kron(np.eye(3), A1)
      A1 = np.c_[A1, np.zeros((np.shape(A1)[0], 1))] # Add column of zeros for s
      b1 = config._gMaxConstraint * np.ones(A1.shape[0])
      constraints.append(LinearConstraint(sparse.csc_array(A1), lb=-b1, ub=b1))
  
  # Constraint on change in gradients (slew rate)
  A2 = sparse.diags_array([1, -2, 1], offsets=[-1, 0, 1], shape=(3 * config.N, 3 * config.N)).toarray()
  A2 = np.c_[A2, np.zeros((np.shape(A2)[0], 1))] # Add column of zeros for s
  b2 = config._sMaxConstraint * np.ones(A2.shape[0])
  constraints.append(LinearConstraint(sparse.csc_array(A2), lb=-b2, ub=b2))

  return constraints


# %%
conf = NOW_config(useMaxNorm=True)
print('------------ Requested timing parameters: ------------ \n')
print('DurPre = {:.3f} DurPost = {:.3f}  DurPi = {:.3f}  [ms]\n\n'.format(conf.durationFirstPartRequested, conf.durationSecondPartRequested, conf.durationZeroGradientRequested))
print('------------   Actual timing parameters:  ------------ \n')
print('DurPre = {:.3f} DurPost = {:.3f}  DurPi = {:.3f}  [ms]\n\n'.format(conf.durationFirstPartActual, conf.durationSecondPartActual, conf.durationZeroGradientActual))

constraints = defineLinearInequalityConstraints(conf)
