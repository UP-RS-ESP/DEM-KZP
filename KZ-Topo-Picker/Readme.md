# DEM-Topo-Picker (DEM-KZP)
Digital Elevation Model (DEM) and KnickZone Picker (KZP) processing

Code developed by Al Neely (abn5031@psu.edu) and Bodo Bookhagen (bodo.bookhagen@uni-potsdam.de)

Latest update: 10-July-2017 (BB)

The code has been tested with Matlab R2012b, R2014b and R2015b, R2016b, R2017a. 
It requires the Statistical Toolbox, the Topotoolbox, and 
export_fig (see below). If the Curve Fitting 
Toolbox is available, it will be used (_smooth_). If no Curve Fitting 
Toolbox is available, _sgolayfilt_ is used (with similar and comparable 
results). It has been tested on Windows 7, Windows 10, and Ubuntu 16.04 LTS.

## Installation

Before running the code and taking advantage of all feautres (i.e.,
high-resolution figures and projected shapefile), you will need to
install the additional items:

1. _TopoToolbox_ (https://topotoolbox.wordpress.com/download/). I suggest
using the github repository
(https://github.com/wschwanghart/topotoolbox): `git clone https://github.com/wschwanghart/topotoolbox`. For additional
information, see: https://topotoolbox.wordpress.com/. Add
the TopoToolbox directory to the Matlab PATH: `addpath(genpath('/home/bodo/Dropbox/Matlab-work/topotoolbox'))`.
The TopoToolbox requires the Image Processing Toolbox and the Statistics
Toolbox (which are mostly installed in academic environment). 
The Mapping Toolbox will come in handy, but is not required. The Curve 
Fitting Toolbox will be used if available and will produce slightly
different fitting parameters (e.g., for 95% confidence Intervals).
Make sure to add the Topotoolbox to the Matlab PATH.

2. _export_fig_ from Mathworks MATLAB Central:
(http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)
for creating the various plots and figures in high-quality PDF format.
Make sure to add the location of the file export_fig.m to the Matlab PATH.

3. Install _ghostscript_ from
http://www.ghostscript.com/download/gsdnld.html . This is required for 
high-res PDF outputs for export_fig.m and highly recommended.

4. _OSGEO4WShell_ or the full GDAL suite (http://www.gdal.org/ or 
http://trac.osgeo.org/osgeo4w/). MAC OS X users should install the GDAL 
complete framwork from http://www.kyngchaos.com/software/frameworks. 
This will allow to generate shapefiles from  
tables that are created during the processing and will allow to
polygonize raster (TIF) files. GDAL is freely available for all operating
systems, and we have tested this on Ubuntu 12.04
(sudo apt-get install gdal), Windows 7 and 8.1 
(http://trac.osgeo.org/osgeo4w/), and
Mac OS X (http://www.kyngchaos.com/software/frameworks). This script will
run well, if you use the default gdal installation options. Otherwise,
you may have to change the directory locations in the parameter file.
We rely on: _ogr2ogr_, _gdalsrsinfo_, _gdal_dem_,
_gdal_polygonize_.

5. Install the _KnickZone-Topo-Picker_ (KZ-Topo-Picker) code via github: `git clone https://github.com/UP-RS-ESP/DEM-KZP/tree/master/KZ-Topo-Picker/Matlab` or download from https://github.com/UP-RS-ESP/DEM-KZP/tree/master/KZ-Topo-Picker/Matla. Add
the KZ-Topo-Picker directory to the Matlab PATH. You can use _addpath_ or the GUI: `addpath('/home/bodo/Dropbox/Matlab-work/DEM-KZP/KZ-Topo-Picker/Matlab')`.

## Processing
The driver or controlling file is the file: [Matlab/run_KZP_example.m](Matlab/run_KZP_example.m). This files sets the path and load the parameter file. You will need to modify to the parameter file with your DEM filename. If you use a 1-m DEM, a good starting point is [Matlab/KZP_parameters_1m_example.m](Matlab/KZP_parameters_1m_example.m). For a 10-m DEM a good starting point is [Matlab/KZP_parameters_10m_example.m](Matlab/KZP_parameters_10m_example.m). These have default parameters for filtering and smoothing to detect knickpoints, but these parameters will need to be modified and optimized for each individual DEM.

Additional driver and parameter files are available at: https://github.com/UP-RS-ESP/DEM-KZP/tree/master/Examples

The codes will write to subdirectories of the directory where den DEM is stored. The subdirectories are:
1. _DEM_geotiff_ This contains several geotif files that are generated from the DEM, including various curvature files, filtered topography, steepness values for channels, specific stream power, relief at various radii, and additional files.

2. _DEM_maps_ contains PDFs of map views with topography, slope, curvature, and steepness.

3. _DEM_plots_ contains log area vs. log slope and chi plots and map views for each basin matching the drainage area threshold (e.g., every basin larger than 10 km2).

4. _DEM_shapefiles_ contains several shapefiles with attributes for steepness, relief, drainage area, and many more parameters.

5. _KZP_plots_ contains knickpoint map views and in longitudinal-river profile, and chi space. This shows some of the filtering steps applied to detect knickpoints in chi plots.

6. _KZP_shapefiles_ contains the shapfiles for knickpoint lips and bases with attributes storing elevation, chi, magnitude, slope, and relief of the knickpoint.

## Examples
For additional parameter examples, see [Examples](Examples).
