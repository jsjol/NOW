# %% 
import numpy as np
from config import NOW_config
# %%
problem = NOW_config()
print('------------ Requested timing parameters: ------------ \n')
print('DurPre = {:.3f} DurPost = {:.3f}  DurPi = {:.3f}  [ms]\n\n'.format(problem.durationFirstPartRequested, problem.durationSecondPartRequested, problem.durationZeroGradientRequested))
print('------------   Actual timing parameters:  ------------ \n')
print('DurPre = {:.3f} DurPost = {:.3f}  DurPi = {:.3f}  [ms]\n\n'.format(problem.durationFirstPartActual, problem.durationSecondPartActual, problem.durationZeroGradientActual))
