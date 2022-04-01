# ARCGen Test Cases

Four datasets are provided with ARCGen to serve as examplars as well as ARCGen tutorials. A brief overview of each dataset and its corresponding literature source is given below. With exception of the Lessley Parabolas, a technical discussion of these datasets is presented in the accompanying paper (Hartlen & Cronin, 2022). 

## Lessley Parabolas
A set of four parabolic curves taken from the work of Lessley et al. (2004) that were used to demonstrate the corridor methodology developed in that paper. The first three curves are relatively similar and are useful to demonstrate the generation of corridors from data that is non-monotonic in at least one axis, as well as termining at different extents in x-axis. The fourth curve is much larger, and helps act as a test case to demonstrate how ARCGen and other corridor methods would react to an outlier. 

The script `TestCase_LessleyParabolas.m` can be considered an introductory tutorial to ARCGen and consists of three different calls of the function `arcgen()` with different input options. 
1. The first call exhibits the basic input options of `acrgen` as well as a quick explaination of `CorridorRes` and `nresamplePoints`.
2. The second call exhibits the effect of mangitude scaling. By setting the option `NormalizeSignals` to `off`, the x and y axis are not scaled the largest bounding box. Because the parabolas are two orders of mangitude smaller in the x-axis than the y-axis, re-parameterization does not capture the true shape of the underlying curve. One will note the resulting corridors exhibit a 'pinch' at the apex of the average parabola as a result. Running ARCGen without `NormalizeSignals` set to `on` is strongly discouraged with `on` being the default option if not otherwise specified. 
3. The third call demonstrates how ARCGen is influenced when an outlier is introduced. In this case, the resulting characteristic average and corridors are skewed, but still capture the underlying characteristic shape of the set of signals. Additionally, this third case demonstrates how additional information from for debugging can be plotted by setting the `Diagnostics` option to `on`. `Diagnostics` has an addition option called  `Detailed` providing more inforation pertaining to the marching squares algorithm. 

Dataset reference: Lessley, D., Crandall, J., Shaw, G., Kent, R., &#38; Funk, J. (2004). A normalization technique for developing corridors from individual subject responses. <i>SAE Technical Papers</i>. https://doi.org/10.4271/2004-01-0288

## Mattucci Ligament Data
This test case generates corridors for the Mattucci and Cronin (2015) cervical spine ligament data. This dataset consists of a number of force-displacement signals collected from the Anterior Longitudinal ligment, Posterior Longitudinal ligmanet, Interspinous Ligament, Capsular Ligment, and Ligamentum Flavum. All curves are truncated at the onset of post-traumatic damage. 

The script `TestCase_MattucciLigaments.m` does not deomstrate any new from the Lessley parameters. However, it does reproduce the plots used in the accompanying paper, Hartlen & Cronin (2022). ARCGen results are plotted against the averaging technique detailed in Mattucci & Cronin (2015)

Dataset reference: Mattucci, S. F. E., &#38; Cronin, D. S. (2015). A method to characterize average cervical spine ligament response based on raw data sets for implementation into injury biomechanics models. <i>Journal of the Mechanical Behavior of Biomedical Materials</i>, <i>41</i>, 251–260. https://doi.org/10.1016/j.jmbbm.2014.09.023

## 