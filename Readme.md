# DEM processing, tectonic geomorphology, and knickpoint detection with DEM-KZP
Automatic Stream Pre-Processing, Selection of Knickzones, and optional calibration of knickzone selection parameters.

*The theoretical background and application of this code has been described in:*
Neely, A., Bookhagen B., Burbank, D.W. (2017): An automated knickzone selection algorithm (KZ-Picker) to analyze transient landscapes: Calibration and validation, JGR Earth Surface, doi:10.1002/2017JF004250, available at:
http://onlinelibrary.wiley.com/doi/10.1002/2017JF004250/full

Code developed by Al Neely (abn5031@psu.edu) and Bodo Bookhagen (bodo.bookhagen@uni-potsdam.de)

Latest update: 9-July-2017 (BB)


**There are two subdirectories containing different Matlab codes with different purposes:**

## KZ-Picker
The directory [KZ-Picker](KZ-Picker) contains a straight-forward Matlab code to calibrate and select knickpoints as described in Neely et al., 2017 (http://onlinelibrary.wiley.com/doi/10.1002/2017JF004250/full). Only Matlab, Topotoolbox, and this code are required. Only Matlab, Topotoolbox and this code are required. **This code is particularly useful for calibrating individual basins and identifying best-fit parameters.**

See [KZ-Picker/Readme.md](KZ-Picker/Readme.md) for additional information.

A box diagram (PDF) shows how the code operates [KZ-Picker_workflow.pdf](KZ-Picker/KZP-Picker_workflow.pdf).

## KZ-Topo-Picker
The directory [KZ-Topo-Picker](KZ-Topo-Picker) contains a more extensive Matlab code that generates several figures, maps, plots, longitudinal river profiles, shapefiles and additional data. It performs an automatic channel-steepness analysis for all subbasins of a given size (parameter) based an log area-log slope plots and chi plots and calculate channel curvature and steepness and can be adapted to process large DEMs.
However, this code requires installing addons (Topotoolbox, Ghostscript, GDAL, and scripts from Mathworks) to be fully functioning. This is intended for users with some experience in DEM processing and knowledge in Matlab. It generates high-resolution, editable PDF figures and longitudinal river profiles and knickpoint plots and mpas. This directory also contains some function that are still developed and improved, for example wind-gap and stream-capture detection algorithms.

Importantly, this Matlab code generates a set of geotif and shapefiles (if GDAL is installed) and allows to further analysis the results in a GIS. For example, vector files with normalized steepness indices (derived from slope-area and chi plots) as well as knickpoint locations are saved. Standard directories where data are save in are _DEM_geotiff_, _DEM_shapefiles_, and _KZP_shapefiles_.

See [KZ-Topo-Picker/Readme.md](KZ-Topo-Picker/Readme.md) for additional information. 

Both codes run on Linux, Mac OSX, and Windows operating systems and have been tested with various versions of Matlab.

## Example Datasets
We provide several example DEMs. For an example on demonstrating the calibration procedure, we refer the detailed manual available at [Examples/KZ_picker_instructions_ABN_7-8-17.pdf](Examples/KZ_picker_instructions_ABN_7-8-17.pdf). This describes the calbriation and processing steps as performed by the [KZ-Picker](KZ-Picker).

An example data including calibration files is available for the the Smugglers catchment for Santa Cruz Island:

1. _smugglers_1m_dem.7z_ at [Examples/SCI-1m-Smugglers](Examples/SCI-1m-Smugglers). An example dataset at containing a 1-m DEM of the northern part of Santa Cruz Island, southern California. The example GeoTIF with calibration files are available at https://www.dropbox.com/s/kpujj3mfezzd34b/smugglers_1m.7z?dl=0 or https://boxup.uni-potsdam.de/index.php/f/38152339


Additional datasets and example - without detailed description, but with paramter files - are available here:

1. _USGS_SCI_10m_UTM11_NAD83.7z_ at [Examples/SCI-USGS-10m](Examples/SCI-USGS-10m) the USGS 10-m DEM for Santa Cruz Island with parameter files. The GeoTIF DEM is available at: https://www.dropbox.com/s/p0ow2sttn3m7tyi/USGS_SCI_10m_UTM11_NAD83.7z?dl=0 https://boxup.uni-potsdam.de/index.php/f/38195855

2. _SCI-Lidar_noveg_DTM_10m_UTM11_NAD83.7z_ at [Examples/SCI-Lidar-10m](Examples/SCI-Lidar-10m) a Lidar data (DTM, digital terrain model with vegetation removed), processed to 10m with parameter files. These data in its originral form were also used in Neely et al., 2017. The GeoTIF DEM is available at: https://www.dropbox.com/s/luac5mms70n1xwh/SCI-Lidar_noveg_DTM_10m_UTM11_NAD83.7z?dl=0 or https://boxup.uni-potsdam.de/index.php/f/38195854

3. _Pozo_DTM_noveg_UTM11_NAD83_1m.7z_ at [Examples/SCI-1m-Pozo](Examples/SCI-1m-Pozo). An example dataset at containing a 1-m DEM of the southwestern part of Santa Cruz Island, southern California. The GeoTIF DEM is available at: https://www.dropbox.com/s/9judp8fy3fp4mbw/Pozo_DTM_noveg_UTM11_NAD83_1m.7z?dl=0 and https://boxup.uni-potsdam.de/index.php/f/38195891
