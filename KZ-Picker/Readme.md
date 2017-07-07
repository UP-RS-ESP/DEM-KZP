%% Processing Digital Elevation Models (DEMs) and generating topometric indices
%
% Automatic Stream Pre-Processing, Selection of Knickzones, and optional
% calibration of knickzone selection parameters
%
% Knickzone selection code by Al Neely (abn5031@psu.edu)
%   (latest update: 4/14/2017)
%
% Automatic stream pre-processing and topometric calculation originally 
% coded by Bodo Bookhagen (bodo.bookhagen@uni-potsdam.de) 11/09/2015; and
% adapted by Al Neely (abn5031@psu.edu) - (9/26/2016)
%
% Knickzon selection parameter calibration code by Al Neely (abn5031@psu.edu)
%   (latest update: 4/14/2017)
%
% The code was tested on a windows 7 operating system with matalb R2015a. 
%
% *Installation*
%
% Before running the code, you will need to install the additional items:
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
% 2. Install the _KZ_picker_ (KZp) code. Move the 'KZ_picker_parameters.m'
% script and 'A_KZ_RUN.m' script to the folder containing your DEM. Add
% the KZ_picker directory containing processing codes to the Matlab PATH.
%
%
% *Processing*
%
% BEFORE using your DEM with this script, make sure that you use a
% projected coordinate system. We strongly suggest to use an UTM or other
% coordinate system with equal areas (or nearly equal areas) for all
% pixels. The geographic coordinate system of most SRTM data has different 
% x and y length of each pixel differ (higher distortion with greater 
% distance from the equator).
%
% This pre-processing script assumes there is a DEM_MAT_parameters.mat file
% in the current directory. You can create this file (a list of parameters)
% with: KZ_picker_parameters.m
% We advise to store this parameter file in the local directory where youe
% DEM GeoTIF file is located. All other Matlab code should be in a separate
% directory and you the paths to this should be included in the Matlab
% PATH.
%
% The script will create up to 8 directories containing output results.
%
