%% Run Topometric calculation and KnickZone Picker (KZ-Picker or KZP)
%
% Automatic calculation of landscape topometric values and extracting of
% knickpoint locations by Bodo Bookhagen (bodo.bookhagen@uni-potsdam.de) on
% 06/01/2017 and 07/07/2017.
% KZ Picker (KZP) script by Al Neely (abn@umail.ucsb.edu) and Bodo
% Bookhagen (bodo.bookhagen@uni-potsdam.de).
%
% See _KZP_parameters_1m_example.m_, _KZP_topometrics.m_,
% _KZP_knickzone_calibration.m_ and _KZP_knickzone_processing.m_
% for more information and how to install the required software and run the
% script.
%
%% (0) Initial Setup and preparation
%Change to working directory
cd /home/bodo/Dropbox/Matlab-work/KZP-examples/SCI-1m-Pozo
%cd C:/Users/bodob//home/bodo/Dropbox/Matlab-work/KZP-examples/SCI-10m-Pozo
%remove existing files:
%Linux:
%!rm -fr DEM_geotiff DEM_maps DEM_plots DEM_shapefiles KZP_csv KZP_plots *.mat *.csv *.crt *.shp *.shx *.prj *.out *.dbf

%Windows:
%!rmdir /s /q DEM_geotiff DEM_maps DEM_plots DEM_shapefiles KZP_csv KZP_plots 
%!del *.mat *.csv *.crt *.shp *.shx *.prj *.out *.dbf

%!7z x Pozo_DTM_noveg_UTM11_NAD83_10m_geotif.7z

%getting latest topotoolbox from: https://github.com/wschwanghart/topotoolbox
%git clone https://github.com/wschwanghart/topotoolbox
%!adjust the following to match your location!
addpath(genpath('/home/bodo/Dropbox/Matlab-work/topotoolbox'))

%getting latest KZ-Picker code from https://github.com/UP-RS-ESP/DEM-KZP/tree/master/KZ-Picker/Matlab
%git clone https://github.com/UP-RS-ESP/DEM-KZP/tree/master/KZ-Topo-Picker/Matlab
%!adjust the following to match your location!
addpath('/home/bodo/Dropbox/Matlab-work/DEM-KZP/KZ-Topo-Picker/Matlab')

%turn Matlab warnings off that may occur when fitting data
warning off

%load and set parameters from .m file (modify as you see fit)
KZP_parameters_1m_example

%% Calculate topometrics from digital elevation model and generate figures
KZP_topometrics

%% Perform calibration of KnickZone Picker (will require manual interaction and calibration files to be created - see manual)
%not fully implemented yet, use KZ-Picker at https://github.com/UP-RS-ESP/DEM-KZP/tree/master/KZ-Picker/Matlab
%KZP_knickzone_calibration

%% Run KZP and identify knickzones
%relies on parameters set in KZP_parameters_10m_example, adjust for
%different DEM spatial resolution and geographic areas
KZP_knickzone_processing
