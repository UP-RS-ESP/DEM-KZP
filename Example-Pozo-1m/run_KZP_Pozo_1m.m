%% Run Topometric calculation and KnickZone Picker (KZ-Picker or KZP)
%
% Automatic calculation of landscape topometric values and extracting of
% knickpoint locations by Bodo Bookhagen (bodo.bookhagen@uni-potsdam.de) on
% 06/01/2017.
% KZ Picker (KZP) script by Al Neely (abn@umail.ucsb.edu) and Bodo
% Bookhagen (bodo.bookhagen@uni-potsdam.de).
%
% See _KZP_parameters_1m_example.m_, _KZP_topometrics.m_,
% _KZP_knickzone_calibration.m_ and _KZP_knickzone_processing.m_
% for more information and how to install the required software and run the
% script.
%
%% Initial Setup
%Change to working directory
%cd /home/bodo/Dropbox/KZP/Pozo-example-1m
cd C:\Users\bodob\Dropbox\KZP\Pozo-example-1m

%Add path that contains KZP code
%addpath('/home/bodo/Dropbox/Matlab-work/KZP/')
addpath('C:\Users\bodob\Dropbox\Matlab-work\KZP\')

%turn Matlab warnings off that may occur when fitting data
warning off

%load and set parameters from .m file (modify as you see fit)
KZP_parameters_1m_example

%% Calculate topometrics from digital elevation model and generate figures
KZP_topometrics

%% Perform calibration of KnickZone Picker (will require manual interaction and calibration files to be created - see manual)
KZP_knickzone_calibration

%% Run KZP and identify knickzones
KZP_knickzone_processing

