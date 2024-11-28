
# ARCGen - General, Feature-Based Corridor Generation

[![View ARCGen - Arc-length-based averaging and statistics on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/116975-arcgen-arc-length-based-averaging-and-statistics)

![ARCGen Logo](./Assets/ARCGen.svg)

> [!NOTE]
> [ARCGen is also available for Python! :snake:](https://github.com/IMMC-UWaterloo/ARCGen-Python)

ARCGen is a general, robust methodology providing feature-based assessment of average response and variability in the form of a characteristic average response and statistical response corridors. In particular, ARCGen is well suited to tackling the challenging types of signals common in biomechanics, such as:
- Monotonic signals that do not share a common termination point
- Highly oscillatory signals, such as those that capture head or limb kinematics
- Hysteretic signals or signals that are non-monotonic in both axes

ARCGen is released under the open-sourced GNU GPL v3 license. No warranty or guarantee of support is provided. The authors hold no responsibility for the validity, accuracy, or applicability of any results obtained from this code.

# Installation

Download the ARCGen Toolbox from the [MathWorks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/116975-arcgen-arc-length-based-averaging-and-statistics) or download and install directly from within MATLAB by selecting "Add-Ons" from the "Home" tab of the main toolbar and searching for "ARCGen" in the Add-On Explorer. Installing from ARCGenthe toolbox ensures that the `arcgen()` function is added to MATLAB's execution path. 

To test that MATLAB is installed correctly, type `open arcgen` into the MATLAB command window. If the `arcgen.m` file opens in the editor, ARCGen is installed correctly. 

ARCGen is provided with several example datasets that explore different types of signals and different ways to use ARCGen. These example datasets are packaged with the downloaded toolbox and in this repository's [`TestCases` folder](https://github.com/IMMC-UWaterloo/ARCGen/tree/main/TestCases).


# Usage
ARCGen is used like a typical MATLAB function with one required input variable (`inputSignals`) and a host of optional name-value pair arguments that provide important control over arc-length re-parameterization and signal registration. An example of ARCGen used in code with typical name-value pair arguments is provided below. 

````matlab
[charAvg, innerCorr, outerCorr] = arcgen(inputSignals, ...
                                         'nResamplePoints', 250, ...
                                         'CorridorRes', 250, ...
                                         'nWarpCtrlPts', 2, ...
                                         'WarpingPenalty', 1e-2, ...
                                         'UseParallel', 'on')
````

Complete documentation of all of ARCGen's input arguments and optional outputs is available both in the code and below. 

## Mandatory Inputs
`inputSignals`: Contains the signals from which a characteristic average and corridors can be computed. Signals must be two-dimensional (i.e. acceleration-time, force-displacement) but do not need to contain the same number of points or be sampled at the same sampling frequency. There is no limit to the number of input signals ARCGen can accommodate, but runtime will increase as the number of signals increases. `inputSignals` must be defined in one of the following three formats. The script ['PreProcessInputSignals.m'](https://github.com/IMMC-UWaterloo/ARCGen/blob/main/PreProcessInputSignals.m) is provided to help users place their data into an appropriate format for use with ARCGen. 
- a [nSignal, 2] dimensional structured array consisting of the following two fields. Field names are case-sensitive
  - data: a two-column array defining the input signal
  - specId: a character array containing an identifier for the associated signal
- a [nSignal, 1] dimensional structured array consisting of a single field
  - data: a two-column array defining the input signal
- a [nSignal, 1] cell array containing 'nSignal' two column arrays corresponding to each input signal. 

## Mandatory Outputs
`charAvg`: A two-column array containing the computed characteristic average

`innerCorr`: A two-column array defining the inner or lower portion of the corridor

`outerCorr`: A two-column array defining the outer or upper portion of the corridor

## Optional Inputs
ARCGen contains many optional inputs that allow the user to control many of the internal parameters ARCGen uses. As with all name-value pair arguments in MATLAB, the order in which these arguments are defined is unimportant. 

`nResamplePoints`: An integer defining the number of points re-parameterized signals will contain. Increasing this value will provide a smoother characteristic average and corridors. Default: 100. 

`CorridorRes`: An integer defining the number of grid points the marching squares algorithm uses to extract the corridor. This sampling grid is automatically calculated based on the most extreme possible corridor values. Increasing this value will improve how accurately the corridor is extracted. Default: 100. 

`NormalizeSignals`: A character array used to control whether magnitude normalization is performed before arc-length re-parameterization. It is highly recommended that this option be enabled. Input options: 'on' (default), 'off'

`EllipseKFact`: A float used to scale the size of the ellipses defined by the confidence region at each point of the characteristic average. A value of 1.0 creates a corridor that is one standard deviation wide at each point of the characteristic average. Default: 1.0 (or corridors of plus and minus one standard deviation).

`MinCorridorWidth`: A float value used to enforce a minimum corridor width based on maximum standard deviation. This option can be useful if corridors do not extend to the beginning or end of the characteristic average due to low variability between input signals. However, enforcing a minimum width is not representative of the actual variability of underlying data. Default: 0.0 (does not enforce a minimum corridor width). 

`nWarpCtrlPts`: An integer defining the number of interior control points used for signal registration. A value of 0 disables signal registration. Default: 2 (some warping). 

`WarpingPenalty`: A float defining the penalty factor used during signal registration. Default: 1e-2. 

`UseParallel`: A character array used to control if the Parallel Computing Toolbox is used to accelerate signal registration and envelope extraction. This option requires MATLAB's Parallel Computing Toolbox. Significantly reduces runtime for signals of 100K+ points or using 500+ resampling points and corridor resolution points. Input options: 'on', 'off' (default).

`Diagnostics`: A character array used to activate diagnostic plots. Useful for tracing anomalous behaviours. Options: 'off' (default), 'on', 'detailed'.

## Optional Outputs
In addition to the three mandatory outputs documented above, ARGCen also has two optional outputs. These optional outputs are position-sensitive in the output array. 

`processedSignalData`: a structured array containing information about the processed signals, such as basic statistics, warping control points, and re-parameterized signals. 

`debugData`: A structured array containing various pieces of information that may be useful for debugging, including raw average and standard deviation information, correlation scores used during signal registration, and other information. 


# Dependencies

ARCGen was developed on MATLAB R2020b. While later versions should work, compatibility with other versions, especially earlier versions, is not guaranteed. ARCGen requires three MATLAB toolboxes. 
- Optimization Toolbox: Required for signal registration
- Mapping Toolbox: Required for corridor extraction
- Parallel Computing Toolbox: (optional) Required to accelerate signal registration and corridor generation

ARCGen uses compiled MEX code for accelerated execution. Compilation was performed for MATLAB. Other operating systems are not supported or gauranteed to work. 

# Referencing

If you use ARCGen in published research, please use the following citation in your research. 

Hartlen D.C. and Cronin D.S. (2022), "Arc-Length Re-Parametrization and Signal Registration to Determine a Characteristic Average and Statistical Response Corridors of Biomechanical Data." *Frontiers in Bioengineering and Biotechnology* 10:843148. doi: 10.3389/fbioe.2022.843148

Bibtex format:
```
@article{Hartlen_Cronin_2022,
 AUTHOR={Hartlen, Devon C. and Cronin, Duane S.},   
 TITLE={Arc-Length Re-Parametrization and Signal Registration to Determine a Characteristic Average and Statistical Response Corridors of Biomechanical Data},      
 JOURNAL={Frontiers in Bioengineering and Biotechnology},      
 VOLUME={10},      
 YEAR={2022},      
 URL={https://www.frontiersin.org/article/10.3389/fbioe.2022.843148},       
 DOI={10.3389/fbioe.2022.843148},       
 ISSN={2296-4185},   
}
```
# Contributing

If you discover any bugs or issues or would like to suggest improvements, please create an [issue on this GitHub repostiory](https://github.com/DCHartlen/ARCGen/issues). You are invited free to submit pull requests to integrate fixes and features directly into ARCGen, although pull requests will be reviewed before integration. 

Anyone is free to fork ARCGen for their work as long as they follow the GNU GPL v3 license requirements.

# Overview of Operation
For a detailed description of how ARCGen operates, please refer to [Hartlen and Cronin (2022)](https://www.frontiersin.org/article/10.3389/fbioe.2022.843148). Therefore, only a high-level overview is presented here. In general, the operation of ARCGen can be divided into three stages: arc-length re-parameterization, signal registration, and statistical analysis. 

## Arc-length Re-parameterization
Arc-length re-parameterization allows ARCGen to handle input signals that are non-monotonic in either measurement axis (a behaviour called hysteresis). Arc-length provides a convenient means to define input points with respect to a strictly monotonic value inherently tied to the shape of the signal. 

Before computing the arc-length of each curve, all signals are scaled to eliminate issues of differing magnitude between the axes. Scaling ensures that the shape of the resulting average and corridors are reflective of the inputted signals. All signals are scaled based on the mean extreme values of both measurement axes to eliminate magnitude differences between axes when calculating arc-length. These scaled values are only used for arc-length calculations and are not used later in signal registration or statistical analysis.

Once signals have been scaled, the arc-length of each signal is computed for each point in the signal. The arc-length is then normalized to the total arc-length of each signal, such that all signals have a normalized arc-length of 1.0. Finally, signals are resampled such that all signals have points at the same normalized arc-lengths. 

## Signal Registration
One of the underlying assumptions of arc-length re-parametrization is that critical features in the signal appear at approximately the same normalized arc-length. However, the resulting characteristic average can be distorted if said features are not perfectly aligned. Additionally, features such as significant variability or noise can dramatically affect the arc-length calculation, changing where critical features occur with respect to normalized arc-length. 

Signal registration is applied to help align critical features of signals. This process introduces a warping function for each input signal that subtly alters how each signal is resampled with respect to arc-length to align critical features. These warping functions (strictly monotonic, piecewise Hermite cubic splines) are determined for all signals simultaneously by maximizing cross-correlation between all signals. To ensure that warping functions do not produce highly skewed, unrealistic responses, a penalty function is used to limit the amount of warping introduced. 

Signal registration is an optional process and is not needed if input signals are very highly correlated or strictly monotonic. While some experimentation is needed to select the best number of control points, a rule of thumb would be to set the number of interior control points to the number of inflection points expected in the characteristic average. 

## Statistical Analysis
Following arc-length re-parameterization and registration, all input signals will have the same number of points and have features aligned with respect to a consistent normalized arc-length. Statistical analysis is undertaken in a point-wise fashion at each normalized arc-length. This analysis assumes points are uncorrelated and distributed according to a two-dimensional normal distribution. Based on this assumption, an elliptical confidence region can be constructed at each normalized arc-length. The length of these ellipses' major and minor axes is proportional to the standard deviation at each arc-length. 

The characteristic average of the input signals is defined as the mean value at each normalized arc-length. The response corridors are the envelope of all ellipses. As there is no closed-form way of extracting this envelope, a marching-squares algorithm is used to extract this envelope numerically. Because the envelope is extracted numerically, it is important that the number of resampling points (`nResamplePoints`) is large enough to ensure that ellipses are sufficiently overlapped to provide a smooth, realistic envelope. Similarly, the resolution of the marching squares grid (`CorridorRes`) should be fine enough to capture the shape of the ellipses correctly. This last feature is similar to ensuring that the mesh of a finite element or computational fluid dynamics simulation is fine enough to resolve features. 

# Change Log

## V 1.5.0
Version 1.5.0 of ARCGen sees incremental improvement over previous versions. In addition to switching from calendar versioning to semantic versioning, ARCGen 1.5.0 now performs signal registration with 2 control points by default if `nWarpCtrlPts` is not otherwise defined. Additionally, the default number of resampling points and corridor grid resolution are both increased to 250 points from 150. All three of these parameters can be changed using the appropriate name-value pair arguments. These changes do not change any of the underlying behaviour of ARCGen, but are intended to ensure new users get better results more quickly. 

Additionally, the `README.md` has been updated to increase clarity and get new users up and running more quickly. 

 - `nResamplePoints` default value changed to 250 from 100. 
 - `CorridorRes` default value changed to 250 from 100.
 - `nWarpCtrlPts` default value changed to 2 from 0 (signal registration disabled)

## R2023a
R2023a introduces a completely revised algorithm for envelope splitting. The revised algorithm is unconditionally stable and operates on the concept of ray-polygon interception. In short, if the envelope of all ellipses found with the marching squares algorithm does not intercept the computed characteristic average, rays are projected from the start and/or end of the characteristic average to determine where the envelope should be split into inner and outer corridors. Previously, a simpler algorithm was used that was shown to be effective in most cases, but was only guaranteed to function in all circumstances. 

Additionally, MATLAB code suggestions or 'intellisense' has been added to provide the user with on-the-fly help on inputs and name-value arguments. 

- envelope splitting now uses a ray-polygon interception algorithm
- Added code suggestions. 

## R2022b
R2022 is a hotfix to correct how corridors are divided between inner and outer components. 

## R2022a
R2022a contains a number of small bug fixes but largely consists of documentation updates. R2022a also represents the first release to be packaged as a MATLAB toolbox for single file installation. 
- 'arcgen.m' updated to handle file paths more robustly. 
- All scripts in `Test Cases` updated to act as ARCGen tutorials

## R2021d
R2021d represents a significant performance improvement over previous versions. Corridor extraction has been completely refactored to reduce runtime in that part of the code by upwards of 10x. Corridor extraction is also far more robust than previous versions, which should eliminate crashes in some use cases. Runtimes associated with signal registration and corridor extraction can be further reduced through the integration of MATLAB's Parallel Computing Toolbox. Additionally, R2021d represents the first release with a formal license, GNU GPL v3
- Added the 'UseParallel' to accelerate signal registration and envelope extraction. Requires the MATLAB Parallel Computing Toolbox. 
- Accelerated signal registration for signals with a very large number of points using pre-compiled MEX code. 
- Completely redeveloped envelop extraction algorithm to significantly reduce computational expense. 
- Envelope splitting to produce inner and outer corridors is now far more robust and should eliminate any 'iIntStart' or 'iIntEnd' that have cropped up in the past. 
- Integration of the GNU GPL v3 license. 
  