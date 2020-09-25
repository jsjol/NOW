# Numerical optimization of gradient waveforms (NOW) for tensor-valued dMRI
A MATLAB package for optimization of gradient waveforms that yield b-tensors of arbitrary shape that can be tailored to pulse sequence timing and hardware restrictions.  

The optimizer supports the following:
* Arbitrary b-tensor shape
* Variable raster resolution
* Asymmetric sequence timing
* Control for 0th-moment balance
* Control for energy consumption/heating
* L2 and max-norm amplitude constraints
* Nulling of concomitant gradient effects
* Nulling of motion encoding

## Getting started
First, download (clone or fork) this repository and open it in MATLAB. You can generate waveforms in a graphical interface by calling `NOW_GUI.m`. To access all optimization controls, you can run the optimization via a script. An example script is provided in `scripted_NOW_example.m`.  

Setting up the optimizer always follows these steps:
1. Create object that specifies optimization problem by calling `pObj = optimizationProblem()`
2. Modify `pObj` to your specification
3. Update the derived parameters of `pObj` by feeding it through the function `pObj = optimizationProblem(pObj)`
4. Call the optimizer run function using `pObj` as input, according to `[result, pObj] = NOW_RUN(pObj)`

The `result` structure contains several fields; among them is the gradient waveform. Notably, the fields `gwf`, `rf` and `dt`, are compatible with the [multidimensional diffusion (md-dMRI) framework](https://github.com/markus-nilsson/md-dmri) format. With the md-dMRI framework installed, an overview of the resulting gradient waveform can be plotted by calling `gwf_plot_all(result.gwf, result.rf, result.dt)`.

## References to NOW and its components
The optimization framework contains contains several sub-functions, all of which are part of the master branch of this repository. Please consider citing the following papers if you use NOW in your research or applications.

* The underlying optimization framework by Sjölund et al. (2015), as described here:  
[Sjölund J, Szczepankiewicz F, Nilsson M, Topgaard D, Westin C-F, and Knutsson H. _Constrained optimization of gradient waveforms for generalized diffusion encoding._ Journal of Magnetic Resonance 261 (2015), 157-168.](https://doi.org/10.1016/j.jmr.2015.10.012)

* Concomitant gradient compensation (Maxwell compensation) by Szczepankiewicz et al. (2019), [patent pending](https://www.freepatentsonline.com/y2020/0284865.html), as described here:  
[Szczepankiewicz F, Westin C‐F, and Nilsson M. _Maxwell‐compensated design of asymmetric gradient waveforms for tensor‐valued diffusion encoding._ Magn Reson Med 82 (2019) 1424–1437](https://doi.org/10.1002/mrm.27828)

* Motion compensation by Szczepankiewicz et al. (2020), as described here:  
[Szczepankiewicz F, Sjölund J, Dall’Armellina E, Plein S, Schneider E J, Teh I, and Westin C-F, _Motion-compensated gradient waveforms for tensor-valued diffusion encoding by constrained numerical optimization._ Magn Reson Med (2020) (in press)]()

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



