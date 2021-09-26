# ARCGen - Arc-length Response Corridor Generation

Biofidelity response corridors are a commonly used technique to assesss the performance of surrogates such as computational models or anthropomorphic test devices while capturing the variablity of experimental data. ARCGen presets a new set of methodologies for computing response corridors and the characteristic average of experimental data capable of accomadating all types of input signals, including experimental data that is time-based, cross-variable, non-monotonic, and hysterytic. ARCGen is distributed, open-sourced, as a single MATLAB function.

This documentation provides information on how to use ARCGen as well as a high-level overview of the methodologies that ARCGen uses. For a more detailed description how ARCGen operates, please refer to Hartlen et al (20XX).

## Usage
ARCGen is distributed as a MATLAB function and makes extensive use of name-value pair arguments to define options.  In its most basic, ARGen can be run with the following manditory inputs and outputs. One can refer to in-code documentation for additional information. 

````matlab
[charAvg, innerCorr, outerCorr] = arcgen(inputSignals)
````

### Manditory Inputs
`inputSignals`: Contains the signals from which a characteristic average and corridors can be computed. Signals must be two-dimensional (i.e. acceleration-time, force-disp), but do not need to contain the same number of points or be sampled at the same sampling frequency. There is no limit to the number of input signals ARCGen can accomadate, but runtime will increase as the number of signals increases. `inputSignals` must be defined in one of the following three formats. 
- a [nSignal, 2] dimensional structured array consisting of the following two fields. Field names are case-sensitive
  - data: a two column array defining the input signal
  - specId: a character array containing an identifier for the associated signal
- a [nSignal, 1] dimensional structured array consisting of a single field
  - data: an two column array defining the input signal
- a [nSignal, 1] cell array containing 'nSignal' two column arrays corresponding to each input signal. 

### Manditory Outputs
`charAvg`: A two-column array containing the computed characteristic average

`innerCorr`: A two-column array defining the inner or lower portion of the corridor

`outerCorr`: A two-column array defnining the outer or upper portion of the corridor

### Optional Inputs
ARCGen contains a large number of optional inputs that allow the user to control many of the parameters ARCGen uses. As with all name-vlaue pair arguments in MATLAB, the order which which these arguments is included is not important. 

`nResamplePoints`: An integer defining the number of points re-parameterized signals will contain. Increasing this value will provide a smoother characteristic average and corridors. Default: 100. 

`CorridorRes`: An integer defining the number of grid points used by the marching squares algorithm to extract the corridor. This sampling grid is automatically calculated based on the most extreme possible corridor values. Increasing this value will improve how well the corridor is extracted. Default: 100. 

`NormalizeSignals`: A character array used to control whether signal normalization is performed prior to arc-length re-parameterization. It is hihgly recommended that this option is enabled. Input options: 'on' (default), 'off'

`EllipseKFact`: A float used to scale the size of the ellipses defined by the confidence region at each point of the characteristic average. A value of 1.0 creates a corridor that is one standard deviation wide at each point of the characterstic average. However, this only corresponds to approximately 39% of variablity for two-diemsnional data. The square root of the chi-squared quantile function can be used to define reponse corridors with respect to p-value. Default: 1.0 (or corridors of plus or minus one standard deviation).

`MinCorridorWidth`: A float value used to enforce a minimum corridor width based on maximum standard deviation. This option can be useful if corridors do not extend to the beginning or end of the charactstisic average due to low variablity between input signals. However, enforcing a minimum width is not representative of the actual underlying data. Default: 0.0 (does not enforce a minimum corridor width). 

`nWarpCtrlPts`: An integer defning the number of interior control points used for signal registration. Default: 0 (disables signal registration)

`WarpingPenalty`: A float defining the penalty factor used during signal registration. Default: 1e-2. 

`Diagnostics`: A character array used to activate diagnostic plots. Useful for tracing anomolous behavours. Options: 'off' (default), 'on', 'detailed'.

### Optional Outputs
In additon to the three manditory outputs documented above, ARGCen also has two optional outputs. These optional outputs are position-sensitive in the output array. 

`processedSignalData`: a structured array containing information about the processed signals, such as basic statistics, warping control points, and re-parameterized signals. 

`debugData`: A structured array containing various peices of information that may be useful for debugging including raw avarage and standard deviation information, correlation scores used during signal registration, and other information. 




