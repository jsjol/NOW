# Numerical optimization of gradient waveforms (NOW) for tensor-valued dMRI
A package for optimization of gradient waveforms that yield b-tensors of arbitrary shape that can be tailored to pulse sequence timing and hardware restrictions. Available in **MATLAB** and **Python**.

The optimizer supports the following:
* Arbitrary b-tensor shape
* Variable raster resolution
* Asymmetric sequence timing
* Control for 0th-moment balance
* Control for energy consumption/heating
* L2 and max-norm amplitude constraints
* Nulling of concomitant gradient effects
* Nulling of motion encoding

## Python

### Installation

Requires Python ≥ 3.10. Install dependencies:

```bash
pip install numpy scipy jax matplotlib
pip install pytest  # for running tests
```

### Quick start

```python
import numpy as np
from now import NOW_config, now_optimize
from now.visualization import plot_result

# 1. Create a configuration
config = NOW_config(
    gMax=80,                    # max gradient amplitude [mT/m]
    sMax=100,                   # max slew rate [T/m/s]
    N=50,                       # number of discretization points
    durationFirstPartRequested=32,   # duration before pause [ms]
    durationSecondPartRequested=27,  # duration after pause [ms]
    durationZeroGradientRequested=8, # pause duration [ms]
    targetTensor=np.eye(3),     # spherical tensor encoding (STE)
    eta=0.9,                    # energy/efficacy balance (0, 1]
)

# 2. Run optimization
result, config = now_optimize(config)

# 3. Inspect results
print(f"b-value: {result.b:.2f} s/mm²")
print(f"B-tensor:\n{result.B}")

# 4. Plot
plot_result(result, config)
```

A complete example is provided in `now_example.py`.

### Configuration parameters

| Parameter | Default | Description |
|---|---|---|
| `targetTensor` | `np.eye(3)` | 3×3 target encoding tensor. `eye(3)` = STE, `diag([1,0,0])` = LTE, `diag([1,1,0])` = PTE |
| `N` | 77 | Number of discretization points. More = smoother waveforms, slower optimization |
| `gMax` | 80 | Maximum gradient amplitude [mT/m] |
| `sMax` | 100 | Maximum slew rate [T/m/s] |
| `durationFirstPartRequested` | 28 | Duration before the zero-gradient pause [ms] |
| `durationSecondPartRequested` | 22 | Duration after the zero-gradient pause [ms] |
| `durationZeroGradientRequested` | 8 | Duration of the zero-gradient pause [ms] |
| `eta` | 1 | Heat dissipation / energy balance in (0, 1] |
| `useMaxNorm` | `False` | Use max-norm instead of L2-norm for gradient amplitude |
| `doMaxwellComp` | `True` | Enable Maxwell (concomitant gradient) compensation |
| `MaxwellIndex` | 100 | Threshold for Maxwell terms [(mT/m)² ms] |
| `enforceSymmetry` | `False` | Force waveform symmetry about the pause |
| `motionCompensation` | `None` | Dict with keys `'order'` and `'maxMagnitude'` (see below) |
| `doBackgroundCompensation` | 0 | 0 = off, 1 = general timing, 2 = specific timing |
| `startTime` | 0 | Time from excitation to first sample [ms] (for `doBackgroundCompensation=2`) |

#### Motion compensation

```python
config = NOW_config(
    motionCompensation={
        'order': [1, 2],           # compensate 1st and 2nd order motion
        'maxMagnitude': [0, 1e-4], # 0 = exact nulling (linear constraint),
    },                             # >0 = allowed deviation (nonlinear constraint)
)
```

#### Result fields

| Field | Units | Description |
|---|---|---|
| `result.b` | s/mm² | b-value |
| `result.B` | s/m² | Full 3×3 b-tensor |
| `result.g` | mT/m | Gradient waveform, shape (N+2, 3) |
| `result.q` | 1/m | q-space trajectory, shape (N, 3) |
| `result.slew` | T/m/s | Slew rate, shape (N+2, 3) |
| `result.kappa` | — | Encoding efficiency |
| `result.gwf` | T/m | Gradient waveform, md-dMRI compatible |
| `result.rf` | — | Spin dephasing direction |
| `result.dt` | s | Time step |

#### Optimization methods

```python
result, config = now_optimize(config, method='SLSQP')       # default, fastest
result, config = now_optimize(config, method='trust-constr') # interior point
```

### Running tests

```bash
cd NOW
python -m pytest tests/ -v
```

The test suite validates the Python implementation against MATLAB reference data at every intermediate computation step: config values, constraint matrices, nonlinear constraint values, and analytical Jacobians. Jacobians are additionally verified against JAX automatic differentiation.

### Correspondence with MATLAB

| MATLAB | Python |
|---|---|
| `optimizationProblem()` | `NOW_config()` |
| `NOW_RUN(problem)` | `now_optimize(config)` |
| `result.gwf` | `result.gwf` |
| `result.b` | `result.b` |
| `fmincon` (SQP) | `scipy.optimize.minimize` (SLSQP) |

Because the MATLAB and Python solvers are different implementations, optimized waveforms will generally differ (different local optima), but both satisfy the same constraints and produce similar b-values.

## MATLAB

### Getting started
First, download (clone or fork) this repository and open it in MATLAB. You can generate waveforms in a graphical interface by calling `NOW_GUI.m`. To access all optimization controls, you can run the optimization via a script. An example script is provided in `scripted_NOW_example.m`.  

Setting up the optimizer always follows these steps:
1. Create object that specifies optimization problem by calling `pObj = optimizationProblem()`
2. Modify `pObj` to your specification
3. Update the derived parameters of `pObj` by feeding it through the function `pObj = optimizationProblem(pObj)`
4. Call the optimizer run function using `pObj` as input, according to `[result, pObj] = NOW_RUN(pObj)`

The `result` structure contains several fields; among them is the gradient waveform. Notably, the fields `gwf`, `rf` and `dt`, are compatible with the [multidimensional diffusion (md-dMRI) framework](https://github.com/markus-nilsson/md-dmri) format. With the md-dMRI framework installed, an overview of the resulting gradient waveform can be plotted by calling `gwf_plot_all(result.gwf, result.rf, result.dt)`.

## References to NOW and its components and extensions
The optimization framework contains contains several sub-functions, all of which are part of the master branch of this repository. Please consider citing the following papers if you use NOW in your research or applications.

* The underlying optimization framework by Sjölund et al. (2015), as described here:  
[Sjölund J, Szczepankiewicz F, Nilsson M, Topgaard D, Westin C-F, and Knutsson H. _Constrained optimization of gradient waveforms for generalized diffusion encoding._ Journal of Magnetic Resonance 261 (2015), 157-168.](https://doi.org/10.1016/j.jmr.2015.10.012)

* Concomitant gradient compensation (Maxwell compensation) by Szczepankiewicz et al. (2019), [patent](https://www.freepatentsonline.com/y2020/0284865.html), as described here:  
[Szczepankiewicz F, Westin C‐F, and Nilsson M. _Maxwell‐compensated design of asymmetric gradient waveforms for tensor‐valued diffusion encoding._ Magn Reson Med 82 (2019) 1424–1437](https://doi.org/10.1002/mrm.27828)

* Motion compensation by Szczepankiewicz et al. (2021), as described here:  
[Szczepankiewicz F, Sjölund J, Dall’Armellina E, Plein S, Schneider E J, Teh I, and Westin C-F, _Motion-compensated gradient waveforms for tensor-valued diffusion encoding by constrained numerical optimization._ Magn Reson Med (2021)](https://onlinelibrary.wiley.com/doi/10.1002/mrm.28551)

* Cross-term compensation by Szczepankiewicz and Sjölund (2021), as described here:  
[Szczepankiewicz F and Sjölund J, _Cross-term-compensated gradient waveform design for tensor-valued diffusion MRI._ Journal of Magnetic Resonance (2021)](https://doi.org/10.1016/j.jmr.2021.106991)

## Review paper on gradient waveform design
Gradient waveforms can be intended for many kinds of purposes and deployed on vastly different hardware. The general design of gradient waveforms for dMRI, with special focus on tensor-valued encoding, has been described in:
* General gradient waveform design by Szczepankiewicz et al. (2021):  
[Szczepankiewicz F, Westin C-F, and Nilsson, M, _Gradient waveform design for tensor-valued encoding in diffusion MRI._ Journal of Neuroscience Methods (2021)](https://doi.org/10.1016/j.jneumeth.2020.109007)

## Free waveform encoding pulse sequence
To run user-defined gradient waveforms a special MRI pulse sequence is usually required. Pulse sequences are available for multiple vendors and scanner software versions. Please refer to the [free waveform (FWF) encoding resource site](https://github.com/filip-szczepankiewicz/fwf_seq_resources) for more information.

## Help improve NOW
We welcome contributions from anyone! NOW is an open-source project that aspires to be community-driven. 
If you have a feature request or discover a bug, please submit an issue or - even better - fix it and make a pull request.
If you plan to make major changes, it might be a good idea to first submit an issue to discuss what you would like to change.

## MATLAB dependencies
The following toolboxes are called during optimization:  
* Curve Fitting Toolbox
* System Identification Toolbox
* Optimization Toolbox
* Simulink Control Design
* Statistics and Machine Learning Toolbox
* Computer Vision Toolbox

List was derived from `pList.Name` by calling: `[fList,pList] = matlab.codetools.requiredFilesAndProducts('scripted_NOW_Example.m');`

## License
This work is published under the BSD 3-Clause License.



