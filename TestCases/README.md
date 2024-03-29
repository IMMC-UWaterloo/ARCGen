# ARCGen Test Cases

Four datasets are provided with ARCGen to serve as examplars as well as ARCGen tutorials. A brief overview of each dataset and its corresponding literature source is given below. With exception of the Lessley Parabolas, a technical discussion of these datasets is presented in the accompanying paper (Hartlen & Cronin, 2022). 

## Lessley Parabolas
A set of four parabolic curves taken from the work of Lessley et al. (2004) that were used to demonstrate the corridor methodology developed in that paper. The first three curves are relatively similar and are useful to demonstrate the generation of corridors from data that is non-monotonic in at least one axis, as well as termining at different extents in x-axis. The fourth curve is much larger, and helps act as a test case to demonstrate how ARCGen and other corridor methods would react to an outlier. 

The script `TestCase_LessleyParabolas.m` can be considered an introductory tutorial to ARCGen and consists of three different calls of the function `arcgen()` with different input options. 
1. The first call exhibits the basic input options of `acrgen` as well as a quick explaination of `CorridorRes` and `nresamplePoints`.
2. The second call exhibits the effect of mangitude scaling. By setting the option `NormalizeSignals` to `off`, the x and y axis are not scaled the largest bounding box. Because the parabolas are two orders of mangitude smaller in the x-axis than the y-axis, re-parameterization does not capture the true shape of the underlying curve. One will note the resulting corridors exhibit a 'pinch' at the apex of the average parabola as a result. Running ARCGen without `NormalizeSignals` set to `on` is strongly discouraged with `on` being the default option if not otherwise specified. 
3. The third call demonstrates how ARCGen is influenced when an outlier is introduced. In this case, the resulting characteristic average and corridors are skewed, but still capture the underlying characteristic shape of the set of signals. Additionally, this third case demonstrates how additional information from for debugging can be plotted by setting the `Diagnostics` option to `on`. `Diagnostics` has an addition option called  `Detailed` providing more inforation pertaining to the marching squares algorithm. 

Dataset reference: 

Lessley, D., Crandall, J., Shaw, G., Kent, R., &#38; Funk, J. (2004). A normalization technique for developing corridors from individual subject responses. <i>SAE Technical Papers</i>. https://doi.org/10.4271/2004-01-0288

## Mattucci Ligament Data
This test case generates corridors for the Mattucci and Cronin (2015) cervical spine ligament data. This dataset consists of a number of force-displacement signals collected from the Anterior Longitudinal ligment, Posterior Longitudinal ligmanet, Interspinous Ligament, Capsular Ligment, and Ligamentum Flavum. All curves are truncated at the onset of post-traumatic damage. 

The script `TestCase_MattucciLigaments.m` does not deomstrate any new from the Lessley parameters. However, it does reproduce the plots used in the accompanying paper, Hartlen & Cronin (2022). ARCGen results are plotted against the averaging technique detailed in Mattucci & Cronin (2015)

Dataset reference: 

Mattucci, S. F. E., &#38; Cronin, D. S. (2015). A method to characterize average cervical spine ligament response based on raw data sets for implementation into injury biomechanics models. <i>Journal of the Mechanical Behavior of Biomedical Materials</i>, <i>41</i>, 251–260. https://doi.org/10.1016/j.jmbbm.2014.09.023

## NBDL 15g Frontal Acceleration Dataset
This dataset consists of head kinematics collected from human volunteers under 15g frontal deceleration (Ewing & Thomas, 1972). Head acceleration and displacement timeseries data is collected in the x and z axis as well as rotation about the y-axis. 

The scripts `TestCase_NBDL_15gFrontal_Accels.m` and `TestCase_NBDL_15gFrontal_Disp.m` demonstrate how signal registration options are used. The option `nWarpCtrlPts` is used to control the number of warping control points. The default value of this option is `0`, which disables signal registration. The option `warpingPenalty` specifies the warping penalty used to control warping. By default, this option is `1e-2`. Users are encouraged to change both `warpingPenalty` and `nWarpCtrlPts` in the provided scripts to understand how the resulting corridors are changed. 

Additionally, both scripts also make use of the parallel computing to accelerate the calculation of the characterisic average and response corridors. Use of the MATLAB Parallel Computing Toolbox is control with the input option `UseParallel` set to `on`. Parallel computing is not enabled if the Parallel Computing Toolbox is not installed. The impact of `UseParallel` becomes more apparent when `nResamplePoints` and `CorridorRes` are set to values greater than 500. ARCGen will terminate with an error if sets this option to `on` if the toolbox is not installed. 

Dataset Reference: 

Ewing, C. L., &#38; Thomas, D. J. (1972). <i>Human Head and Neck Response to Impact Acceleration.</i> Naval Aerospace Medical Research Lab Pensacola Fl.

National Highway Traffic Safety Administration. (2017). <i>Biomechanics Test Database</i>. https://www.nhtsa.gov/research-data/databases-and-software

## Kroell Thoracic Impact Database
This database consists of a set of force-displacement signals collected from PMHS undergoing impact to the thorax from a 50 lb pendulum travelling at 16 mph. This dataset was highly influential in the development of the first generation of anthropomorphic test devices. However, the load-unload or hysteretic nature of the signals makes traditional averaging and corridor methods difficult. Corridors have generally be defined using the "eyeball average" recorded in Lobdell et al (1972). 

The script `TestCase_KroellThoraxCorridors.m` demonstrates the effectiveness of ARCGen on signals that are non-monotonic in both axis, in addition to the signal registration and parallel computing introduced in the NBDL dataset. For hysteretic signals, the number of warping points can be chosen by examining the x and y components of the signals with respect to normalized arc-length. In this case, three warping control points are chosen for the single key feature in x-component and two critical features in the y-component. 

Dataset References:

Kroell, C. K., Schneider, D. C., &#38; Nahum, A. M. (1971). Impact Tolerance and Response of the Human Thorax. <i>SAE Technical Papers</i>. https://doi.org/10.4271/710851

Lobdell, T. E., Kroell, C. K., Schneider, D. C., Hering, W. E., &#38; Hahum, A. M. (1972). Impact Response of the Human Thorax. In W. F. King &#38; H. J. Mertz (Eds.), <i>Human Impact Response: Measurement and Simulation</i> (pp. 201–245). Springer Science+Business Media.