"""Example: optimize a spherical tensor encoding (STE) gradient waveform."""
import numpy as np
from now import NOW_config, now_optimize
from now.visualization import plot_result

# Create problem configuration (matching MATLAB scripted_NOW_Example.m)
config = NOW_config(
    gMax=80,               # mT/m
    sMax=100,              # T/m/s
    N=50,
    durationFirstPartRequested=32,
    durationSecondPartRequested=27,
    durationZeroGradientRequested=8,
    targetTensor=np.eye(3),  # Spherical tensor encoding
    eta=0.9,
    MaxwellIndex=100,
    doMaxwellComp=True,
)

print(f"Requested timing: pre={config.durationFirstPartRequested}ms, "
      f"pause={config.durationZeroGradientRequested}ms, "
      f"post={config.durationSecondPartRequested}ms")
print(f"Actual timing:    pre={config.durationFirstPartActual:.2f}ms, "
      f"pause={config.durationZeroGradientActual:.2f}ms, "
      f"post={config.durationSecondPartActual:.2f}ms")

# Run optimization
result, config = now_optimize(config, method='SLSQP', max_attempts=5)

print(f"\nb-value: {result.b:.4f} s/mm²")
print(f"Encoding efficiency κ: {result.kappa:.4f}")
print(f"B-tensor:\n{result.B}")

# Plot results
plot_result(result, config)
