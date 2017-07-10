%% Parameters file for KZ picker

% Follow brief instructions provided to input parameter values

%% Name Parameters
DEM_basename = 'smugglers_1m_dem'; 
% Filename before '.tif' extension of DEM to analyze

%% Minimum Drainage Area
Min_DA_threshold = 1e5;
% Specify the minimum drainage area for the fluvial network (m^2) 
%(Slope/area break from constant slope to concave-up stream profile)

Min_trib_size = 125; % 10m DEM: 15     30m DEM: 11 (cells)
% Removes reaches that are slightly > than Min_DA_threshold and quickly join with 
% a larger stream order. Make slightly larger than smoothing window size

min_str_gradient = 0.0001;
% set minimum river gradient. This value will be forced on all river
% gradients less than min_str_gradient. For example, all 0-gradient
% sections of a river will be set to this value to force a flow direction.
% Needed if areas of DEM are filled

%% Smoothing and fildering parameter values:

smoothing_option = 1; % set = 1 to smooth long profile, 0 for no smoothing

smoothing_window = 125;  % 10m DEM: 15     30m DEM: 11 (cells)
% Savitzky-Golay smoothing window size (Cells) used to filter geomorphic
% noise and small amplitude convexities from raw chi-elev plot

sgolayfilt_order = 11;  
% polynomial order of Savitzky-Golay filter, must be smaller than smoothing_window

min_kp_size1 = 1; 
% minimum knickpoint mmagnitude (m) before combining closely spaced knickzones

lumping_distance_upstream = 75; % (DEFAULT OPTION) (m)
% combines knickzones spaced closer than this streamwise-distance apart (m)
% value should reflect spacing of waterfalls or steps in individual knickzones

%% Filtering small knickzones by height
kp_magnitude_filter_option = 1; % DEFAULT OPTION (m)
% uses knickzone magnitude (m) as measurement to remove small knickzones
%(measures steepness (Ksn) anomaly: the elevation drop across knickzone 
% after normalizing for the average stream ksn)

kp_relief_filter_option = 0; % (m) Intuitive option, simply filters 
% knickzones by the relief across the knickzone reach but biased by 
% stream concavity moving upstream (knickzones further upstream have larger
% background slopes and therefore should have higher relief). NOT detrended
% relative to background steepness

% NOTE:  knickzone relief systematically increases upstream as background 
% stream gradients increase upstream (why kp_magnitude is prefered for
% filtering)

min_kp_size2_magnitude = 5;  % Default option (m)
% minimum knickpoint magnitude (m)
% used if line 46 = 1

min_kp_size2_relief = 10; % filter by drop in elevation 
% (biased by changes in background channel gradient)
% minimum knickpoint relief (m) 
% used if line 51 = 1

%
min_kp_slope = 3;  % minimum knickzone slope (degrees)
% increase if you want to select only steep waterfalls for example (note
% smoothing profile will strongly affect measurement of knickzone slope)

%% Additional Parameters
theta_bf_option = 1; % set = 1 if you want to fit theta to co-linearize 
% (tributaries for chi plots) Set = 0 if you want to use ref. concavity
theta_ref = -0.45; % ref. concavity (usually -0.45)

% option to generate longitudal profiles with KZ's plotted
% (does this for every tributary analyzed.. so a lot of figures
% potentially)
create_long_prof_4_tributaries = 0; % =1 for yes do this
% this won't occur if only calibrating parameters

% scaling markersize for knickzones on Plots (display handle)
kp_plot_size = 0.8; % (may change depending on size of largest knickzones)
% controls the size of the symbols plotted on the output figures. 
% (change from 0-1 smaller #s = smaller symbols)

%% Calibration Option
Calibration_option = 0; 
% ^^ change to 1 if calibrating smoothing and filtering parameters

do_you_have_calibration_KZ_bases = 1; 
% Set = 0 if you only have positions and heights of knickzone lips in your
% calibration dataset. Set = 1 if you have lips and bases

% If calibrating alogrithm parameters, specify the name in front of the
% .csv calibration files that contain the knickzone positions and reliefs
KZ_lips_calib_fname = 'Smug_lips_1m.csv'; % PAGE 5 (4c. in manual)
KZ_bases_calib_fname = 'Smug_bases_1m.csv';
% ^^ these are tables containing the position of the calibration knickzones
% function csvread will read these files in starting at row 2 (assuming row
% 1 is a header)

KZ_calib_easting_column_num = 3; 
% column number in calibration knickzone .csv file containin easting
% coordinates. 10 is default for tables constructed from profiler toolbar

KZ_calib_northing_column_num = 4; 
% colmun number in calibration knickzone .csv file containin northing
% coordinates. 11 is default for tables constructed from profiler toolbar

KZ_calib_relief_column_num = 5; 
% column number in calibration knickzone .csv file containin relief
% coordinates. 13 is default for tables constructed from profiler toolbar
% and following instruction manual. If you only have KZ lips, you don't
% have to include relief measurements in your calibration file.

calibration_snapping_tolerance = 10; % (m) % 30m for 30m DEMs
% spatial tolerance when referencing calibration knickzones to streamobj.

% usually there is some small misfit between the position of the
% calibration knickzones exported using the profiler toolbar and the
% position of the streamobj network.  calibration knickzones more than this
% distance away from the streamobj network are discarded.  these 
% calibration knickzones may have been selected on a tributary that was not
% analyzed or the calibration knickzones were selected at a stream-wise
% position with a lower contriubing drainage area than Min_DA_threshold 
% for the streamOBJ (line 13 in parameters script).

% error tolerance btwn calibration and algorithm KZ boundaries (m)
error_radius = 50;
% spatial tolerance used to determine true positives, false positives, 
% & false negatives. May increase with coarser DEM resolution.


%% Input parameter values to cycle through. Enter same number of values for each parameter
% Probably use < 8 different values for each parameter or else may take a
% while to run.

SG_smoothing_calib = [51 101 125 175]; % sgolay smoothing window size (cells)
% Input the range of smoothing windows you want to loop through
% ^^ will change with different DEM resolutions 
%(decrease with lower resolution DEMs)

% Input the range of knickzone lumping windows that you want to loop through
lumping_search_distance_calib = [25 50 75 100]; % (m) upstream
% if two knickpoints are within this distance from one another they are 
% combined into one larger knickzone extending from the base of the
% downstream knickzone to the lip of the upstream knickzone

% set a minimum knickpoint size (pre-lumping) NOTE: the algorithm results 
% are usually not very sensitive to this parameter
min_kp_size1_calib = [0.25 0.5 1 1.5]; 
% (before running the closely spaced knickpoint combining function) 
% meters of elev drop not explained by stream gradient

% set a minimum knickpoint size for after the lumping function
min_kp_size2_calib = [3 5 7 10]; 
% minimum knickpoint size 2 (EITHER magntiude or relief depending on inputs
% in lines 36 and 41)
% (after running the closely spaced knickpoint combining function) 
% meters of elev drop not explained by stream gradient



%% I'M NOT SURE WHAT THE 6 VARIABLES BELOW DO??? DO WE NEED THESE?
%% Code records your parameter inputs and saves them for processing
DEM_fname = strcat(DEM_basename, '.tif');  % filename of DEM
TIF_DIR_basename = strcat('geotiff/', DEM_basename); % filename of output tifs

DEM_MAT_fname = strcat(DEM_basename, '.mat');  % filename for DEM in Matlab MAT file (NOT USED)
DEM_FIL_fname = strcat(TIF_DIR_basename, '_FIL.tif');  % filename of filled DEM (NOT USED)
DEM_FAC_fname = strcat(TIF_DIR_basename, '_FAC.tif');  % filename of FlowAccumulation Grid (NOT USED)
DEM_dbasin_fname = strcat(TIF_DIR_basename, '_DBASIN.tif');  % filename of DrainageBasin Grid (NOT USED)


save -v7.3 DEM_MAT_parameters.mat *fname Min_DA_threshold kp_plot_size...
    DEM_basename min_str_gradient theta_bf_option theta_ref... 
    Min_trib_size kp_magnitude_filter_option smoothing_window...
    kp_relief_filter_option min_kp_size2_magnitude...
    min_kp_size2_relief min_kp_size1 min_kp_slope ...
    DEM_FIL_fname DEM_FAC_fname DEM_MAT_fname DEM_dbasin_fname DEM_fname...
    calibration_snapping_tolerance Calibration_option min_kp_size2_calib...
    min_kp_size1_calib lumping_search_distance_calib SG_smoothing_calib...
    error_radius sgolayfilt_order...
    lumping_distance_upstream KZ_lips_calib_fname KZ_bases_calib_fname...
    smoothing_option do_you_have_calibration_KZ_bases KZ_calib_easting_column_num... 
    KZ_calib_northing_column_num KZ_calib_relief_column_num create_long_prof_4_tributaries...