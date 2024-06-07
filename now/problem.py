# %% 
import numpy as np
from scipy.optimize import LinearConstraint, NonlinearConstraint, minimize
from scipy import sparse
from config import NOW_config
import matplotlib.pyplot as plt

def  get_linear_constraints(config):

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

  # Allocate matrix Aeq that will act on a single component
  Aeq = np.zeros((2 + config.zeroGradientAtIndex.size, config.N)) 
  # + nnz(problem.motionCompensation.linear) + 1 * problem.doBackgroundCompensation, problem.N);

  # Require start and end in q-space origin (echo condition)
  Aeq[0, 0] = 1
  Aeq[1, config.N - 1] = 1

  # Require zero gradient at the specified indices
  firstDerivativeMatrix = np.diff(np.eye(config.N)).T
  Aeq[2 + np.arange(config.zeroGradientAtIndex.size, dtype=int), :] = \
    firstDerivativeMatrix[config.zeroGradientAtIndex.astype(int),:]

  if config.motionCompensation:
    raise NotImplementedError
  
  if config.doBackgroundCompensation:
     raise NotImplementedError

  if config.enforceSymmetry:
     raise NotImplementedError
     
  Aeq = np.kron(np.eye(3), Aeq)
  Aeq = np.c_[Aeq, np.zeros((np.shape(Aeq)[0], 1))] # Add column of zeros for s
  beq = np.zeros(np.shape(Aeq)[0])
  constraints.append(LinearConstraint(sparse.csc_array(Aeq), lb=beq, ub=beq))

  return constraints

class tensor_encoding_fun:
  def __init__(self, config):
     self.targetTensor = config.targetTensor
     self.tolIsotropy = config._tolIsotropy
     self.N = config.N

     integrationWeights = np.ones(self.N)
     integrationWeights[0] = 0.5
     integrationWeights[-1] = 0.5
     integrationWeights = integrationWeights / (self.N - 1)
     self.integrationsweights = integrationWeights.reshape(-1, 1)
     
  def __call__(self, x):
    q = x[0:-1]
    s = x[-1]
    Q = q.reshape(self.N, 3)
    B = Q.T @ (self.integrationWeights * Q)

    out = (np.linalg.norm(B - s * self.targetTensor, 'fro')**2
           - (s * self._tolIsotropy)**2)
    return out

def instantaneous_gradient_norm_squared(x):
  Q = x[0:-1].reshape(-1, 3)
  g = np.diff(Q, axis=0)
  out = np.sum(g**2, axis=1)
  return out

def power_constraint(x):
  Q = x[0:-1].reshape(-1, 3)
  g = np.diff(Q, axis=0)
  out = np.array([g[:, i].T @ g[:, i] for i in [0, 1, 2]])
  return out

def get_nonlinear_constraints(config):
  
  constraints = [NonlinearConstraint(tensor_encoding_fun(config),
                                     lb=-np.inf, ub=0.)]
  
  if not config.useMaxNorm:
    constraints.append(
      NonlinearConstraint(instantaneous_gradient_norm_squared,
                          lb=-np.inf, ub=config._gMaxConstraint**2))

  constraints.append(
    NonlinearConstraint(power_constraint,
                        lb=-np.inf, ub=config._integralConstraint))
  return constraints

def objective(x):
  f = -x[-1]
  grad = np.zeros_like(x)
  grad[-1] = -1
  return f, grad

# %%
conf = NOW_config(useMaxNorm=False)
print('------------ Requested timing parameters: ------------ \n')
print('DurPre = {:.3f} DurPost = {:.3f}  DurPi = {:.3f}  [ms]\n\n'.format(conf.durationFirstPartRequested, conf.durationSecondPartRequested, conf.durationZeroGradientRequested))
print('------------   Actual timing parameters:  ------------ \n')
print('DurPre = {:.3f} DurPost = {:.3f}  DurPi = {:.3f}  [ms]\n\n'.format(conf.durationFirstPartActual, conf.durationSecondPartActual, conf.durationZeroGradientActual))

constraints = get_linear_constraints(conf) + \
              get_nonlinear_constraints(conf)

# %%
res = minimize(objective, conf.x0, method=None, jac=None,
               bounds=None, constraints=constraints,
               tol=None, options=None)


# %%
