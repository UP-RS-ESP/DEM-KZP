# DEM processing, tectonic geomorphology, and knickpoint detection with DEM-KZP
Automatic Stream Pre-Processing, Selection of Knickzones, and optional calibration of knickzone selection parameters.

*The theoretical background and application of this code has been described in:*
Neely, A., Bookhagen B., Burbank, D.W. (2017): An automated knickzone selection algorithm (KZ-Picker) to analyze transient landscapes: Calibration and validation, JGR Earth Surface, doi:10.1002/2017JF004250, available at:
http://onlinelibrary.wiley.com/doi/10.1002/2017JF004250/full

Code developed by Al Neely (abn5031@psu.edu) and Bodo Bookhagen (bodo.bookhagen@uni-potsdam.de)

Latest update: 8-July-2017 (BB)

**There are two subdirectories containing different Matlab codes with different purposes:**

## KZ-Picker
The directory [KZ-Picker](KZ-Picker) contains a straight-forward Matlab code to calibrate and select knickpoints as described in Neely et al., 2017 (http://onlinelibrary.wiley.com/doi/10.1002/2017JF004250/full). Only Matlab, Topotoolbox, and this code are required. Only Matlab, Topotoolbox and this code are required. This code is particularly useful for calibrating individual basins and identifying best-fit parameters. See [KZ-Picker/Readme.md](KZ-Picker/Readme.md) for additional information. A box diagram (PDF) shows how the code operates [KZ-Picker_workflow.pdf](KZ-Picker/KZP-Picker_workflow.pdf).

## KZ-Topo-Picker
The directory [KZ-Topo-Picker](KZ-Topo-Picker) contains a more extensive Matlab code that generates several figures, maps, plots, longitudinal river profiles, shapefiles and additional data. It performs an automatic channel-steepness analysis for all subbasins of a given size (parameter) based an log area-log slope plots and chi plots and calculate channel curvature and steepness and can be adapted to process large DEMs.
However, this code requires installing addons (Topotoolbox, Ghostscript, GDAL, and scripts from Mathworks) to be fully functioning. This is intended for users with some experience in DEM processing and knowledge in Matlab. It generates high-resolution, editable PDF figures and longitudinal river profiles and knickpoint plots and mpas. This directory also contains some function that are still developed and improved, for example wind-gap and stream-capture detection algorithms. See [KZ-Topo-Picker/Readme.md](KZ-Topo-Picker/Readme.md) for additional information. 

Both codes run on Linux, Mac OSX, and Windows operating systems and have been tested with various versions of Matlab.

## Example Dataset
We provide example DEMs with parameter and calibration files:
1. _smugglers_1m_dem.7z_, an example dataset at containing a 1-m DEM of the northern part of Santa Cruz Island, southern California. This example GeoTIF with calibration files are available at https://www.dropbox.com/s/kpujj3mfezzd34b/smugglers_1m.7z?dl=0 or https://boxup.uni-potsdam.de/index.php/f/38152339

2. _pozo_1m_dem.7z_, an example dataset at containing a 1-m DEM of the southwestern part of Santa Cruz Island, southern California. 
