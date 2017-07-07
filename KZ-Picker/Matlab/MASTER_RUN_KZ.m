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

%% 4 directories are generated when calibration_option = 0:
%   option uses inputted parameters to select and meaasure knickzones
%   (DEFAULT OPTION)
%
%   Output directories:
%
% 1. \SA_plots
%       ^ contains slope-area plots for each drainage basin analyzed
%
% 2. \maps_and_chi_plots
%       ^ each basin has a chi plot and longitudinal profile with the 
%       knickzone lips and bases plotted
%
% 3. \csv_KZ_databases
%       ^ contains .csv files (for KZ lips and for KZ bases).
%       These .csv files store northing, easting, elev, distance upstream,
%       chi coordinate, knickzone geometeries and parameters used during
%       knickzone selection.  
%
%       ^^ use these .csv files to plot knickzone positions using arcmap!
%       csv files are created for each individual basin and for all
%       knickzones selected
%
%      SEE top of 'E_KZ_knickzone_selection' script for details on what
%      information is stored in each .csv database
%
% 4. \trib_long_profile_figs
%   ^ optional! but if 'create_long_prof_4_tributaries' parameter = true,
%   then this directory will store output long profiles with knickzones
%   plotted for every single tributary analyzed. THESE long profile figures
%   NEED to be hyperlinked to shapefiles in ArcGIS so users can click on
%   tributary shapefiles IN ARC and have these long-profiles appear
%
%
%% 4 directories are generated when calibration_option = 1:
%   This option calibrates knickzone selection smoothing and filtering
%   parameters (IF DESIRED AND SPECIFIED)
%
% 5. \calib_bases
%       ^ contains .csv files that store knickzone bases attributes for
%       the best fit parameter combinations
%
% 6. \calib_lips
%       ^ contains .csv files that store knickzone lips attributes for
%        the best fit parameter combinations
%
% 7. \calib_database_outputs
%       ^ contains 2 .csv files that store accuracy metrics and parameter
%       values for each parameter combination simulation.  See instruction
%       manual for details on what information is stored in these files
%
% 8. \calib_figure_outputs
%       ^ contains figures from calibration, including:
%       a. Calib_KZ_sOBJ_post-ref : position of calibration knickzones
%       before and after referencing them to the StreamOBJ file (will show
%       if some calibration knickzones were too far from the stream
%       network)

%       b. Calib_KZ_long_prof_post-ref : longitudinal profile plotting the
%       calibration knickzone positions

%       c. SG_window_Accuracy : accuracy score as a function of changing
%       the Savitzky golay window size

%       d. Lumping_distance_accuracy : accuracy score as a function of changing
%       the lumping disnance size

%       e. min_kz_height_pre_lumping : accuracy score as a function of changing
%       the minimum knickzone size before combining closely spaced
%       knickzones
%
%       f. min_kz_height_post_lumping : accuracy score as a function of changing
%       the minimum knickzone size after combining closely spaced knickzones
%
%       g. BF_map_RS_only : map of best fit knickzones and parameter values
%       only weighting the spatial agreement of the calibration and
%       algorithm knickzones
%
%       h. BF_long_profile_RS_only: longitudinal profile of best fit knickzones and parameter values
%       only weighting the spatial agreement of the calibration and
%       algorithm knickzones
%
%       IF CALIBRATION KNICKZONE BASES WERE PROVIDED
%       i. BF_map_RSG_: map of best fit knickzones and parameter values
%       weighting the accuracy of measured geometry of knickzones as well
%
%       j. BF_long_profile_RSG: longitudinal profile of best fit knickzones
%       and parameter values only weighting the accuracy of measured geometry 
%       of knickzones as well as spatial agreement
%
%% (0) setup paths and environment
%getting latest topotoolbox from: https://github.com/wschwanghart/topotoolbox
%git clone https://github.com/wschwanghart/topotoolbox
%!adjust the following to match your location!
addpath(genpath('/home/bodo/Dropbox/Matlab-work/topotoolbox'))

%getting latest KZ-Picker code from https://github.com/UP-RS-ESP/DEM-KZP/tree/master/KZ-Picker/Matlab
%git clone https://github.com/UP-RS-ESP/DEM-KZP/tree/master/KZ-Picker/Matlab
%!adjust the following to match your location!
addpath('/home/bodo/Dropbox/Matlab-work/DEM-KZP/KZ-Picker/Matlab')

%edit PARAMETERS_INPUTS_KZ_picker.m to set input and output filenames and
%parameters. Here, parameters are set for the example from the smugglers_1m_dem
%this example file is available at https://www.dropbox.com/s/kpujj3mfezzd34b/smugglers_1m.7z?dl=0
%or https://boxup.uni-potsdam.de/index.php/f/38152339
PARAMETERS_INPUTS_KZ_picker

%% (1) Load parameter file and prepare DEM
%
%change to directory with data and start processing
cd /home/bodo/Dropbox/Matlab-work/KZP-examples/SCI-1m-smugglers
%extract archive
%!7z x smugglers_1m.7z
orgfolder = pwd;

A_KZ_topometrics_1load_v1

%% (2) Calculate flow direction and flow accumulation and choose basins to analyze
%
B_KZ_topometrics_2preprocessing_v1

%% (3) Mask DEM according to drainage basins selected for processing
% (REQUIRES USER INPUT)
%
C_KZ_drainage_basin_plot_and_select

%% (4) Use selections from data cursor to mask draiange basins of interest 
% Calculate Ksn and Chi coordinate for stream-cell in each drainage basin
% Export slope-area figures
D_KZ_masking_DBs

%% (5-Calibration) If calibration option is selected (NOTE: THIS STEP IS SLOW AND USUALLY ONLY NEEDS TO BE RUN ONCE!)
oldfolder = orgfolder;

if Calibration_option == 1
    % run calibration scripts and do not run the processing scripts
    E_calib_KZ_calibration_I
    
    if do_you_have_calibration_KZ_bases == 1
        F_calib_KZ_calibration_comparison_II
    else
        F_calib_KZ_calibration_comparison_II_only_calib_lips
    end
end

%% if calibration is not selected (default) run processing scripts with
% default/inputted smoothing parameters
%NOTE: The example calibrates the data and DOES NOT run the knickzone
%detection. Set Calibration_option = 0 to run the knickzone detection
if Calibration_option == 0
    % (5) Select Knickzones in each drainage basin
    % Select knickzones from Chi/elev data of each stream and tributary
    E_KZ_knickzone_selection
end

fprintf('\n Finished!\n')
