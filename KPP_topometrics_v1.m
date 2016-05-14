%% Processing Digital Elevation Models (DEMs) and generating topometric indices
%
% Automatic Stream Pre-Processing and extracting of topometric parameters
% By Bodo Bookhagen (bodo.bookhagen@uni-potsdam.de) 11/09/2015
%
% The code has been tested with R2012b, R2014a, R2014b, R2015b. Note that
% older versions of Matlab (e.g., R2012b) use a slightly different set of
% instructions and code segments (see variable MATLABV in parameter file).
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
% *Processing*
%
% BEFORE using your DEM with this script, make sure that you use a
% projected coordinate system. We strongly suggest to use an UTM or other
% coordinate system with equal areas (or nearly equal areas) for all
% pixels. The geographic coordinate system that most SRTM data come in are
% no projections and the x and y length of each pixel differ (higher
% distortion with greater distance from the equator).
%
% This pre-processing script assumes there is a DEM_MAT_parameters.mat file
% in the current directory. You can create this file (a list of parameters)
% with: knickpoints_parameters.m
% We advise to store this parameter file in the local directory where youe
% DEM GeoTIF file is located. All other Matlab code should be in a separate
% directory and you the paths to this should be included in the Matlab
% PATH.
%
% This script will check if files exists (geotiff, shapefile, plots) and
% will not-regenerate them, unless you specify the REGEN flag in the
% parameter file. Thus, if you would like to re-process your
% files and re-generate figures and shapefiles, set REGEN to 1 (0 by
% default).
%
% The script will write the following shapefiles:
%
% 1. <filename>_all_MS_proj.shp
%
% 2. <filename>_1e+06_MS_proj.shp (or other drainage area thresholds)
%
% 3. <filename>_db_stats_CNTR.shp
%
% 4. <filename>_db_stats_OUT.shp
%
% 5. <filename>_ridgecrest_MS.shp
%
% 6. <filename>_ridgecrest_MS_Dy_mean_1std.shp
%
% 7. <filename>_ridgecrest_MS_Dy_mean_2std.shp
%
% 8. <filename>_ridgecrest_MS_Dy_parab_1std.shp
%
% 9. <filename>_ridgecrest_MS_Dy_cosh2_1std.shp
%
% 10. <filename>_ridgecrest_MS_Dy_cosh4_1std.shp
%
% 11. in the directory _trunk_ there will be files
% <filename>_MS_trunk_ksn_db1_proj.shp with consecutive numbers for each
% identified trunk stream in the study area.
%
% 12. in the directory _STO_ (stream order) there will be files
% <filename>_STO_02_db_stats_CNTR.shp and
% <filename>_STO_02_db_stats_OUT.shp for the centroid and outlet
% coordinate. The number after _STO_ indicates the stream order and can be
% defined in the parameter file.
%
%
% The files in the folder shapefiles for the _ridgecrest_ contain
% attributes to decipher ridgecrest anomalies. These are:
%
% 1. 1ID
%
% 2. 2X - X Coordinate
%
% 3. 3Y - Y Coordinate
%
% 4. 4Area - Drainage area in m^2
%
% 5. 5dist_Dy1 or 5dist_Dy1 - Distance along ridgecrest profile where
% \Delta elevation exceeded a threshold value (1 or 2 sigma)
%
% 6. 6elevi_nrm - normalized elevation (between 0 - outlet and 1 - highest
% point)
%
% 7. 7Dy1_mn, 7Dy1_parab, 7Dy1_cosh2, 7Dy1_cosh4, or 7Dy2_mn - \Delta
% elevation values ouytside a threshold with respect to the overall
% 1-std.dev. mean (_mn), parabula, cosh *2, cosh * 4, or 2-std. dev. mean.
%
%
% The shapefiles in the folder shapefiles/STO will have the following
% attributes. These exist for the outlet coordinate (*_OUT) and the basin's
% centroid coordinate (*_CNTR).
%
% 1: Basins ID [1ID]
%
% 2: Basin Centroid X Coordinate [2Centr_X]
%
% 3: Basin Centroid Y Coordinate [3Centr_Y]
%
% 4: Basin Outlet X Coordinate [4BasinO_X]
%
% 5: Basin Outlet Y Coordinate [5BasinO_Y]
%
% 6: Steepness Index (ksn), derived from slope area plot from all
% tributaries and trunk streams with theta=-0.45 [6ksn045]
%
% 7: theta value derived from slope area plot from all
% tributaries and trunk streams [7theta_sa]
%
% 8: Steepness Index (ksn) derived from trunk-stream data only, using
% slopearea [8ks_trsa]
%
% 9: theta value derived from trunk-stream data only, using
% slopearea [9the_trsa]
%
% 10: Steepness Index (ksn), derived from trunk-stream calculated with
% weighted regression (or robust fitting) [10ks_rfit]
%
% 11: theta value, derived from trunk-stream calculated with weighted
% regression (or robust fitting) [11the_rfit]
%
% 12: R2 values, derived from trunk-stream calculated with
% weighted regression (or robust fitting) (columns 11the_rfit and 10ks_rfit)
% [12r2_rfit]
%
% 13: RMSE values, derived from trunk-stream calculated with
% weighted regression (or robust fitting) (columns 10the_rfit and 11ks_rfit)
% [13rmse_rf]
%
% 14: Steepness Index 95% confidence Interval, lower bound, derived from
% trunk-stream calculated with weighted regression (or robust fitting)
% [14ks_f_ci1]
%
% 15: Steepness Index 95% confidence Interval, upper bound, derived from
% trunk-stream calculated with weighted regression (or robust fitting)
% [15ks_f_ci2]
%
% 16: theta values 95% confidence Interval, lower bound, derived from
% trunk-stream calculated with weighted regression (or robust fitting)
% [16the_fci1]
%
% 17: theta value 95% confidence Interval, upper bound, derived from
% trunk-stream calculated with weighted regression (or robust fitting)
% [17the_fci2]
%
% 18: Steepness Index p value, derived from trunk-stream calculated with
% weighted regression (or robust fitting) [18prfit_a]
%
% 19: theta value p value, derived from trunk-stream calculated with
% weighted regression (or robust fitting) [19prfit_b]
%
% 20: theta value, derived from trunk-stream calculated with NO weighted
% regression [20the_nrf]
%
% 21: Steepness Index, derived from trunk-stream calculated with NO weighted
% regression [21ksn_nrf]
%
% 22: R2 value, derived from trunk-stream calculated with NO weighted
% regression [22r2_norf]
%
% 23: RMSE value, derived from trunk-stream calculated with NO weighted
% regression [23rmse_nrf]
%
% 24: theta value, derived from trunk-stream calculated with weighted
% regression (or robust fitting) for log-binned data [24the_lsrf]
%
% 25: Steepness Index, derived from trunk-stream calculated with weighted
% regression (or robust fitting) for log-binned data [25ks_lsrf]
%
% 26: R2 value, derived from trunk-stream calculated with weighted
% regression (or robust fitting) for log-binned data [26r2_lsrf]
%
% 27: RMSE value, derived from trunk-stream calculated with weighted
% regression (or robust fitting) for log-binned data [27rmse_lrf]
%
% 28: Steepness Index, p value, derived from trunk-stream calculated with
% weighted regression (or robust fitting) for log-binned data [28plsf_a]
%
% 29: theta value, p value, derived from trunk-stream calculated with
% weighted regression (or robust fitting) for log-binned data [29plsf_b]
%
% 30: adjusted ks with theta taken from chi analysis from the trunk stream
% [30ks_adj]
%
% 31: ksn taken from chi analysis from trunk stream [31ksn_chi]
%
% 32: theta (m/n) taken from chi analysis from the trunk stream [32the_chi]
%
% 33: R2 value from theta (m/n) chi analysis from the trunk stream
% [33the_cR2]
%
% 34: Drainage Area in m2 [34DA_m2]
%
% 35: Drainage Area in km2 [35DA_km2]
%
%
%% (1) Load parameter file and prepare DEM
%
KPP_topometrics_1load_v1

%% (2) Calculating flow direction, flow accumulation, relief and other topographic metrics
%
KPP_topometrics_2preprocessing_v1

%% (3) Calculate  area-slope, ksn, and chi values for all drainage basins contained in this DEM (e.g., draining to ocean/outlets/edge). Generate plots and maps.
%
KPP_topometrics_3SA_chi_analysis_v1

%% (4) Process individual StreamOrder (STO) basins
%
KPP_topometrics_4STO_analysis_v1