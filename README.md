# DEM-KPP
Digital Elevation Model (DEM) and KnickPointPicker (KPP) Analyzer

% Code developed by Al Neely and Bodo Bookhagen 10/12/2015, significantly
% modified March 2016
%
% The code has been tested with Matlab R2012b, R2014b and R2015b. 
% It requires the Statistical Toolbox, the Topotoolbox, and 
% export_fig (see below). If the Curve Fitting 
% Toolbox is available, it will be used (_smooth_). If no Curve Fitting 
% Toolbox is available, _sgolayfilt_ is used (with similar and comparable 
% results).
%
% *Installation*
% 
% Before running the code and taking advantage of all feautres (i.e.,
% high-resolution figures and projected shapefile), you will need to
% install the additional items:
%
% 1. _TopoToolbox_ (https://topotoolbox.wordpress.com/download/). I suggest
% using the github repository
% (https://github.com/wschwanghart/topotoolbox). For additional
% information, see: https://topotoolbox.wordpress.com/.
% The TopoToolbox requires the Image Processing Toolbox and the Statistics
% Toolbox (which are mostly installed in academic environment). 
% The Mapping Toolbox will come in handy, but is not required. The Curve 
% Fitting Toolbox will be used if available and will produce slightly
% different fitting parameters, e.g. (for 95% confidence Intervals).
% Make sure to add the Topotoolbox to the Matlab PATH.
%
% 2. _export_fig_ from Mathworks MATLAB Central:
% (http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)
% for creating the various plots and figures in high-quality PDF format.
% Make sure to add the location of the file export_fig.m to the Matlab PATH.
%
% 3. Install _ghostscript_ from
% http://www.ghostscript.com/download/gsdnld.html . This is required for 
% high-res PDF outputs for export_fig.m and highly recommended.
%
% 4. _OSGEO4WShell_ or the full GDAL suite (http://www.gdal.org/ or 
% http://trac.osgeo.org/osgeo4w/). MAC OS X users should install the GDAL 
% complete framwork from http://www.kyngchaos.com/software/frameworks. 
% This will allow to generate shapefiles from  
% tables that are created during the processing and will allow to
% polygonize raster (TIF) files. GDAL is freely available for all operating
% systems, and we have tested this on Ubuntu 12.04
% (sudo apt-get install gdal), Windows 7 and 8.1 
% (http://trac.osgeo.org/osgeo4w/), and
% Mac OS X (http://www.kyngchaos.com/software/frameworks). This script will
% run well, if you use the default gdal installation options. Otherwise,
% you may have to change the directory locations in the parameter file.
% We rely on: _ogr2ogr_, _gdalsrsinfo_, _gdal_dem_,
% _gdal_polygonize_.
%
% 5. Install the _KnickPointPicker_ (KPP) code and subdirectories and add
% the KPP directory to the Matlab PATH.
%
%
% This Matlab code takes MAT file written by
% _KPP_topometrics_v1.m_ and calculates knickpoints and
% additional attributes from _chi-plot_ analysis. It will write several
% shapefiles and csv files that can be read into any GIS software.
%
% You will have to run _knickpoints_preprocessing_v1.m_ at least once before
% running this script.
%
% Parameters for this script can be adjusted/changed in
% _knickpoints_parameters_X.m_.
%
% This file will output the following four shapefiles:
%
% 1: <filename>_kp_bases.shp (knickpoint bases for all tributaries)
%
% 2: <filename>_kp_lips.shp (knickpoint lips for all tributaries)
%
% 3: <filename>_kp_bases_trunk.shp (knickpoint bases for trunk stream)
%
% 4: <filename>_kp_lips_trunk.shp (knickpoint lips for trunk stream)
%
%
% The code developes a database/csv file/shapefile for the bottom and top
% points of all knickzones. Each database has the following attributes with
% their respective attribute name in brackets [].
%
% Knickpoint Database information [shapefile attribute]
%
% 1: knickpoint # [1kp_id]
%
% 2: stream id # [2stream_id]
%
% 3: tributary id # [3trib_id]
%
% 4: slope of stream, bfl slope: elev/chi [4sl_str]
%
% 5: chi coordinate [5chi]
%
% 6: elevation (m) [6elev_m]
%
% 7: knickpoint magnitude, detrended elevation drop (m) [7kp_magnt]
%
% 8: Knickpoint Relief, height of knickpoint (not detrended (m)) [8kp_Rel]
%
% 9: easting (meters utm) [9Easting_m]
%
% 10: northing (meters utm) [10North_m]
%
% 11: upstream drainage area (m2) [11DA_m2]
%
% 12: dist upstream (m) [12kp_DFM]
%
% 13: knickpoint slope (elev drop/distance upstream) [13kp_slp]
%
% 14: elevation drop (detrended elevation drop/distance upstream, m)
% [14kp_slp_d]
%
% 15: knickpoint slope (detrended elevation drop/chi) [15kp_slp_dt/c]
%
% 16: knickpoint length (distance usptream, m) [16kp_len_m]
%
% 17: sgolay smoothing window size (grid cells) [17sgol_smv]
%
% 18: knickpoint lumping search window size (grid cells) [18kp_lm_ws]
%
% 19: minimum knickkpoint size pre-lumping [19kp_pr_lu]
%
% 20: minimum knickpoint size post-lumping (final minimum knickpoitn size)
% [20kp_pt_lu]
%
% 21: minimum steepness anomaly [21kp_stp_a]
%
% 22: minimum stream size for analysis (cells) [22strea_sz]
%
%
%% (1) Load DEM file and other data from preprocessing
%
KPP_processing_1load_v1

%% (2) Iterate through all basin tributaries and extract knickpoints
%
KPP_processing_2kpp_tribs

%% (3) Iterate through all trunk streams and pull out knickpoints
%
KPP_processing_3kpp_trunks

%% (4) Generate figures for each basin
%
KPP_processing_4kpp_mkfigs
