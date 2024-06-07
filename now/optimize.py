# %%
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from problem import objective, get_constraints
from config import NOW_config

# %%
conf = NOW_config(N=50, useMaxNorm=True, eta=1)
print('------------ Requested timing parameters: ------------ \n')
print('DurPre = {:.3f} DurPost = {:.3f}  DurPi = {:.3f}  [ms]\n\n'.format(conf.durationFirstPartRequested, conf.durationSecondPartRequested, conf.durationZeroGradientRequested))
print('------------   Actual timing parameters:  ------------ \n')
print('DurPre = {:.3f} DurPost = {:.3f}  DurPi = {:.3f}  [ms]\n\n'.format(conf.durationFirstPartActual, conf.durationSecondPartActual, conf.durationZeroGradientActual))


# %%
# The only applicable methods in scipy.opimize are SLSQP and trust-constr.

method = None

if method == 'SLSQP':
  # Halts with positive directional derivative for linesearch
  res = minimize(objective,
                 conf.x0,
                 method=method,
                 jac=True,
                 constraints=get_constraints(conf, method),
                 options= {'maxiter': 5000, 'disp': True})
else:
  # Horribly slow, but runs. Uses the method in:
  # Byrd, Richard H., Mary E. Hribar, and Jorge Nocedal. 1999. An interior point algorithm for large-scale nonlinear programming. SIAM Journal on Optimization 9.4: 877-900.
  res = minimize(objective,
                conf.x0*1e3,
                method='trust-constr',
                jac=True,
                constraints=get_constraints(conf),
                options= {'maxiter': 100000, 'disp': True,
                          'sparse_jacobian': False, 'verbose': 2,
                          'initial_constr_penalty': 1,
                          'initial_tr_radius': 1,
                          'initial_barrier_parameter': 0.1})


# %%
q = res.x[:-1].reshape(-1, 3)
g = np.diff(q, axis=0) / conf._dt
s = res.x[-1]
t = np.linspace(0.5 * conf._dt, conf.totalTimeActual - 0.5 * conf._dt, endpoint=True, num=conf.N - 1)

plt.plot(t, g[:, 0], label='x')
plt.plot(t, g[:, 1], label='y')
plt.plot(t, g[:, 2], label='z')
plt.legend()
plt.xlabel('Time [ms]')
plt.ylabel('Magnetic field gradient')
