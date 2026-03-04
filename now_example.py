#######
# Exameple script for optimization of gradient waveforms 
# Based the Matlab implementation https://github.com/jsjol/NOW
# Written by Jakob Eriksson and Oliver Näslund for Project course at Uppsala University
#######

from optimizationProblem import NOW_config
from optimize import optimize
import numpy as np

def main():
    N = 100 #Raster resolution
    target = np.eye(3) #Target tensor, default isotropic

    # Define problem, similar to Matlab implementation
    prob = NOW_config(N=N,
                      targetTensor=target,
                      durationFirstPartRequested = 28,
                      durationSecondPartRequested = 22,
                      durationZeroGradientRequested = 8,
                      gMax = 80,
                      sMax = 100,
                      useMaxNorm=False,
                      doMaxwellComp=True,
                      eta = 0.9,
                      motionCompensation={'order': [1,2], 'maxMagnitude': [0, 1e-4]})

    # Insert problem formulation into optimization routine
    opt = optimize(problem=prob) 

    # Run the optimization routine
    # early_stopping = True, sets a different termination condition which focuses on tolerance on b-value while perserving feasibility
    opt.run(early_stopping=True)

    # Plot results
    # animate = True, if subiterations wants to be visualized
    opt.plot_results(animate=False)

    # Get arrays of interest
    q, g, b, B = opt.get_results()


if __name__ == "__main__":
    main()
    