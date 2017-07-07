# Processing Digital Elevation Models (DEMs) and generating topometric indices

Automatic Stream Pre-Processing, Selection of Knickzones, and optional calibration of knickzone selection parameters.
*The theoretical background and application of this code has been described in:*
Neely, A., Bookhagen B., Burbank, D.W. (2017): An automated knickzone selection algorithm (KZ-Picker) to analyze transient landscapes: Calibration and validation, JGR Earth Surface, doi:10.1002/2017JF004250, available at:
http://onlinelibrary.wiley.com/doi/10.1002/2017JF004250/full

Code developed by Al Neely (abn5031@psu.edu) and Bodo Bookhagen (bodo.bookhagen@uni-potsdam.de)
Latest update: 7-July-2017

The code was tested on a Ubuntu 16.04 LTS and Windows 7 operating system with Matalb R2015a and Matlab R2016b


## Installation
Before running the code, you will need to install the additional items:

1. _TopoToolbox_ (https://topotoolbox.wordpress.com/download/). I suggest using the github repository (https://github.com/wschwanghart/topotoolbox). For additional information, see: https://topotoolbox.wordpress.com/. The TopoToolbox requires the Image Processing Toolbox and the Statistics Toolbox (which are mostly installed in academic environment). The Mapping Toolbox will come in handy, but is not required. The Curve Fitting Toolbox will be used if available and will produce slightly different fitting parameters, e.g. (for 95% confidence Intervals). Make sure to add the Topotoolbox to the Matlab PATH.

2. Install the _KZ_Picker_ (KZ-Picker) code.  Add the KZ_Picker directory containing processing codes to the Matlab PATH.


## Processing

1. Edit the 'PARAMETERS_INPUTS_KZ_picker' to set input and output filenames. Then use MASTER_RUN_KZ.m to process the DEM. There are additional information in the beginning of the MASTER_RUN_KZ.m file.

## Examples
There exists an example dataset at containing a 1-m DEM of a part of Santa Cruz Island, southern California: smugglers_1m_dem.
This example file is available at https://www.dropbox.com/s/kpujj3mfezzd34b/smugglers_1m.7z?dl=0 or https://boxup.uni-potsdam.de/index.php/f/38152339

## Additional information
BEFORE using your DEM with this script, make sure that you use a projected coordinate system. We strongly suggest to use an UTM or other
coordinate system with equal areas (or nearly equal areas) for all pixels. The geographic coordinate system of most SRTM data has different x and y length of each pixel differ (higher distortion with greater 
distance from the equator) and will need to be reprojected.

This script assumes there is a DEM_MAT_parameters.mat file in the current directory. You can create this file (a list of parameters) with: PARAMETERS_INPUTS_KZ_picker. Edit this file for adjusting input, output and other processing parameters.
We advise to store this parameter file in the local directory where your DEM GeoTIF file is located. All other Matlab code should be in a separate directory and you the paths to this should be included in the Matlab PATH.

This script will create up to 8 directories containing output results.

