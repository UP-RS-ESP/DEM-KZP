# DEM processing, tectonic geomorphology, and knickpoint detection with DEM-KZP
Automatic Stream Pre-Processing, Selection of Knickzones, and optional calibration of knickzone selection parameters.
*The theoretical background and application of this code has been described in:*
Neely, A., Bookhagen B., Burbank, D.W. (2017): An automated knickzone selection algorithm (KZ-Picker) to analyze transient landscapes: Calibration and validation, JGR Earth Surface, doi:10.1002/2017JF004250, available at:
http://onlinelibrary.wiley.com/doi/10.1002/2017JF004250/full

Code developed by Al Neely (abn5031@psu.edu) and Bodo Bookhagen (bodo.bookhagen@uni-potsdam.de)
Latest update: 7-July-2017


There are two subdirectories containing different Matlab codes with different purposes:


The directory (KZ-Picker) contains a straight-forward Matlab code to calibrate and select knickpoints as described in Neely et al., 2017. Only Matlab, Topotoolbox and this code are required.

The directory KZ-Topo-Picker contains a more extensive Matlab code that generates several figures, maps, plots, longitudinal river profiles, shapefiles and additional data. However, this code requires installing addons (Ghostscript, GDAL, Topotoolbox, and scripts from Mathworks) to be fully functioning. This is intended for users with experience in DEM processing and knowledge in Matlab.

Both codes run on all operating systems and have been tested with various versions of Matlab.
