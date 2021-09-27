# ARCGen - Arc-length Response Corridor Generation

Biofidelity response corridors are commonly used to assess the performance of surrogates such as computational models or anthropomorphic test devices while capturing the variability of experimental data. ARCGen presents a set of methodologies for computing response corridors and the characteristic average of experimental data capable of accommodating all types of input signals, including experimental data that is time-based, cross-variable, non-monotonic, and/or hysteretic. ARCGen is distributed as a single MATLAB function.

This document provides information on how to use ARCGen as well as a high-level overview of the methodologies that ARCGen uses. For a more detailed description of how ARCGen operates, please refer to Hartlen et al. (202x).

ARCGen is released under the open-sourced GPL v3 license. _TODO: add specifics_. 

# Usage
ARCGen is distributed as a MATLAB function and makes extensive use of name-value pair arguments to define options. In its most basic, ARCGen can be run with the following mandatory inputs and outputs. While all mandatory and optional arguments are described below, more information is available in in-code documentation. 

````matlab
[charAvg, innerCorr, outerCorr] = arcgen(inputSignals)
````

## Manditory Inputs
`inputSignals`: Contains the signals from which a characteristic average and corridors can be computed. Signals must be two-dimensional (i.e. acceleration-time, force-displacement) but do not need to contain the same number of points or be sampled at the same sampling frequency. There is no limit to the number of input signals ARCGen can accommodate, but runtime will increase as the number of signals increases. `inputSignals` must be defined in one of the following three formats. The script 'PreProcessInputSignals.m' is provided to help users place their data into an appropriate format for use with ARCGen. 
- a [nSignal, 2] dimensional structured array consisting of the following two fields. Field names are case-sensitive
  - data: a two column array defining the input signal
  - specId: a character array containing an identifier for the associated signal
- a [nSignal, 1] dimensional structured array consisting of a single field
  - data: an two column array defining the input signal
- a [nSignal, 1] cell array containing 'nSignal' two column arrays corresponding to each input signal. 

## Manditory Outputs
`charAvg`: A two-column array containing the computed characteristic average

`innerCorr`: A two-column array defining the inner or lower portion of the corridor

`outerCorr`: A two-column array defining the outer or upper portion of the corridor

## Optional Inputs
ARCGen contains many optional inputs that allow the user to control many of the internal parameters ARCGen uses. As with all name-value pair arguments in MATLAB, the order in which these arguments are defined is unimportant. 

`nResamplePoints`: An integer defining the number of points re-parameterized signals will contain. Increasing this value will provide a smoother characteristic average and corridors. Default: 100. 

`CorridorRes`: An integer defining the number of grid points the marching squares algorithm uses to extract the corridor. This sampling grid is automatically calculated based on the most extreme possible corridor values. Increasing this value will improve how accurately the corridor is extracted. Default: 100. 

`NormalizeSignals`: A character array used to control whether magnitude normalization is performed before arc-length re-parameterization. It is highly recommended that this option is enabled. Input options: 'on' (default), 'off'

`EllipseKFact`: A float used to scale the size of the ellipses defined by the confidence region at each point of the characteristic average. A value of 1.0 creates a corridor that is one standard deviation wide at each point of the characteristic average. However, this only corresponds to approximately 38% of variability for two-dimensional data. The square root of the chi-squared quantile function can be used to define response corridors with respect to p-value. Default: 1.0 (or corridors of plus and minus one standard deviation).

`MinCorridorWidth`: A float value used to enforce a minimum corridor width based on maximum standard deviation. This option can be useful if corridors do not extend to the beginning or end of the characteristic average due to low variability between input signals. However, enforcing a minimum width is not representative of the actual variability of underlying data. Default: 0.0 (does not enforce a minimum corridor width). 

`nWarpCtrlPts`: An integer defining the number of interior control points used for signal registration. Default: 0 (disables signal registration)

`WarpingPenalty`: Afloat defining the penalty factor used during signal registration. Default: 1e-2. 

`Diagnostics`: A character array used to activate diagnostic plots. Useful for tracing anomalous behaviours. Options: 'off' (default), 'on', 'detailed'.

## Optional Outputs
In addition to the three mandatory outputs documented above, ARGCen also has two optional outputs. These optional outputs are position-sensitive in the output array. 

`processedSignalData`: a structured array containing information about the processed signals, such as basic statistics, warping control points, and re-parameterized signals. 

`debugData`: A structured array containing various pieces of information that may be useful for debugging, including raw average and standard deviation information, correlation scores used during signal registration, and other information. 

# Overview of Operation
For a detailed description of how ARCGen operates, please refer to Hartlen et al. (202x). Only a high-level overview is presented here. 

The operation of ARCGen can be broken into three stages: arc-length re-parameterization (from which ARCGen draws its name), signal registration, and statistical analysis. 

## Arc-length Re-parameterization
Arc-length re-parameterization is the critical feature that allows ARCGen to handle input signals that are non-monotonic in both axes (a behaviour called hysteresis). Arc-length provides a convenient means to define input points with respect to a strictly monotonic value inherently tied to the shape of the signal. 

Before computing the arc-length of each curve, all signals are scaled to eliminate issues of differing magnitude between the axes. This scaling is based on the mean extreme values of all input signals, such that the relative magnitude of each scaled signal relative to one another is maintained. These scaled values are only used for arc-length calculations and are not used later in signal registration or statistical analysis. 

Once signals have been normalized, the arc-length of each signal is computed for each point in the signal. The arc-length is then normalized to itself, such that all signals have a normalized arc-length of 1.0. Finally, signals are resampled such that all signals have points at the same normalized arc-lengths. 

## Signal Registration
One of the underlying assumptions of arc-length re-parametrization is that critical features in the signal appear at approximately the same normalized arc-length. However, if said features are not perfectly aligned, an average can skew or smear the resulting characteristic average. Additionally, features such as significant variability or noise can dramatically affect arc-length calculation, changing where critical features occur with respect to normalized arc-length. 

Signal registration (or curve registration) can be applied to help align critical features of signals. This process introduces a warping function for each input signal that subtly alters how each signal is re-sampled with respect to arc-length to align critical features. These warping functions (strictly monotonic, piecewise Hermite cubic splines) are determined for all signals simultaneously by maximizing cross-correlation between all signals. To ensure that warping functions do not produce highly skewed, unrealistic responses, a penalty function is used to limit the amount of warping introduced. 

Signal registration is an optional process and is not needed if input signals are very highly correlated or strictly monotonic. While some experimentation is needed to select the best number of control points, a rule of thumb would be to set the number of interior control points to the number of inflection points expected in the characteristic average. 

## Statistical Analysis
Following arc-length reparameterization, all input signals will have the same number of points at the same normalized arc-length. If signal registration has been performed, the registered points will be used for statistical analysis. Statistical analysis is undertaken in a point-wise fashion at each normalized arc-length. This analysis assumes points are uncorrelated and distributed according to a two-dimensional normal distribution. Based on this assumption, an elliptical confidence region can be constructed at each normalized arc-length. The length of the major and minor axes of these ellipses are proportional to the standard deviation at each arc-length. However, unlike a one-dimensional confidence region, where a plus and minus one standard deviation region will account for 68% of variability, the same two-dimensional elliptical region only accounts for 38%. To control the size of the region based on variance, the quantile function of the chi-squared distribution at the desired variance or p-value can be used to scale the size of the ellipse (optional input option `EllipseKFact`).

The characteristic average of the input signals is defined as the mean value at each normalized arc-length. The response corridors are the envelope of all ellipses. As there is no closed-form way of extracting this envelope, a marching-squares algorithm is used to extract this envelope numerically. Because the envelope is extracted numerically, it is important that the number of resampling points (`nResamplePoints`) are large enough to ensure that ellipses are sufficiently overlapped to provide a smooth, realistic envelope. Similarly, the resolution of the marching squares grid (`CorridorRes`) should be fine enough to capture the shape of the ellipses correctly. This last feature is similar to ensuring that the mesh of a finite element or computational fluid dynamics simulation is fine enough to resolve features. 

# References
Hartlen et al. (202x) is currently being prepared. 

# License
_TODO: Add license closer to release date._
