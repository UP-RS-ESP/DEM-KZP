%% Parameter file for KZP
%First section contains parameters for preprocessing, followed by
%parameters for Knickzone Picker

%% Preprocessing
KZP_parameters = struct;
KZP_parameters.DEM_basename = 'Pozo_DTM_noveg_UTM11_NAD83_1m'; 
% we assume that there exists a TIF file <DEM_basename>.tif
% no need to add .tif here
% Make sure that the TIF file is in UTM or other equal-area projected
% coordinates

KZP_parameters.DEM_basename_no_underscore = strrep(KZP_parameters.DEM_basename,'_','\_');
KZP_parameters.theta = -0.45; % reference concavity
KZP_parameters.manual_select_basin = 0; %set to 1 if basin to be processed will be selected manually

KZP_parameters.area_threshold = 1e5;
% set drainage area threshhold in m^2. Drainage areas smaller than this
% value are not considered and no channels are calculated for drainage
% areas smaller than this value. You can investigate good values for  
% this using a log area/log slope plot to determine the start (drainage
% area) of a negative power law relationship.
% *NOTE*: If you change this parameter after an initial processing, you
% will have to remove the *_HYD.mat and *_STR.mat files

KZP_parameters.min_drainage_area_to_process = 1e5;
% Additional drainage area value in m^2 that identifies value for
% statistical calculations. No statistics and topometrics with drainage
% area smaller than this value will be calculated. This is useful if you
% want to process large DEMs, but would like to obtain topometrics only for
% a reduced dataset. Can be set to KZP_parameters.area_threshold or to any
% value larger than KZP_parameters.area_threshold. For example, stream can
% be calculated for KZP_parameters.area_threshold = 5e3 but you decide to
% only calculate statistics for stream with more than 1e5 with
% KZP_parameters.min_drainage_area_to_process = 1e5. This parameter will
% also affect the basins and basin sizes that will be used to calculate
% knickpoints.

KZP_parameters.min_dbasins_stats_to_process = KZP_parameters.min_drainage_area_to_process;
% Additional drainage area parameter for slope-area and chi plots
% (topometrics) calculation. Set to
% KZP_parameters.min_drainage_area_to_process. This is most useful and
% valid only for separate drainage basins, for example a DEM from an island
% with separate drainage basins draining into the ocean. Useful values are
% 1e7 (10km2) or 1e8 (100 
% km2) or 1e9 (1000km2)

KZP_parameters.stream_order = [3 4];
% Stream order for which to create individual basin plots. All basins equal
% this streamorder will be processed for topometric statistics. A stream
% order of 3 for a 30-m DEM includes catchments between 7-300 km2 (mean 45), 
% a stream order of 4 (30-m DEM): 28-1100 km2 (mean: 200 km2). You can also
% give the stream order in vector format. If you give multiple stream
% orders, each number in the vector will generate a separate Shapefile with
% drainage basins and statistics.

KZP_parameters.relief_values_m = [10 25 50];
% Local relief (or radius relief). This vector lists the search radius for
% relief. A value of 1000 indices a 1-km-radius relief calculation. 3
% calculations are performed. Please note that setting these to higher
% values will take a lot of time to calculate/process.
% If you leave this empty (e.g., []), no relief is being calculated

KZP_parameters.str_area1 = 1e6; % Lower drainage area to calculate streamnetwork (STREAMobj) from and generate shapefile
KZP_parameters.str_area2 = 1e6; % Upper drainage area to calculate streamnetwork (STREAMobj) from and generate shapefile
% These two drainage areas indicate the size of catchments to be considered
% when generating shapefiles. Only catchments larger than str_area1
% (usually 1e6 m2) are used for generating a shapefile. For large DEMs,
% this will generate a large (and time-consuming) shapefile. str_area2
% gives the minumum drainage area of a second shapefile. This can be 1-2
% orders of magnitude larger and usually generates a much smaller (and
% faster-to-read) shapefile.

% In order to generate comparable fits with similar log-bin spacing, we
% will need to define the min. and max. drainage area over which you would
% like to perform the slope-area regression. We force this to be the same
% when generating log-bins to make sure the spacing between various
% catchments remains the same. 3rd number indicates number of bins to be
% created (usually 50 or 100). For the regression, we rely on the median of
% the fitted bins.
KZP_parameters.min_max_DA_fits = [1e5 1e8 100];

KZP_parameters.segL = 10; % segment length in meters for smoothing shapefile output, often spatial resolution * 10

KZP_parameters.min_str_gradient = 0.001;
% set minimum river gradient. This value will be forced on all river
% gradients less than min_str_gradient. For example, all 0-gradient
% sections of a river will be set to this value.

% If this is set to 1, relief, curvature, DEM filtering, and
% hillshade DEMs are being calculated. If set to 0, no additional
% calculations are being performed.
KZP_parameters.RELIEF_CURVATURE = 1; 

%If this is set to 1 ridgecrest maps are being generated. This is
%experimental and off by default (=0).
KZP_parameters.RIDGECREST = 0;

KZP_parameters.REGEN = 0; % regenerate all files (1 = Yes, 0 = No)

KZP_parameters.show_figs = 1; %Generate figures? (1 = Yes, 0 = No, 
% 2 = create figures, but no processing (assumes that all MAT files were
% already created - useful for updating figures)
% This will generate figures as PDFs and JPGs in the subfolder maps
% You will need to have installed export_fig from the Mathworks repository
% (free of charge).

KZP_parameters.PaperType_size = 'A4'; % Sets paper type and size for export_fit: 'A4' or 'letter'

%% Parameters for KnickZone-Picker (KZP)

% (0) General KZP parameters
KZ_databases_csv_dir = 'KZP_database_csv';
KZ_databases_shapefile_dir = 'KZP_database_shapefile';
KZ_figs_png_dir = 'KZP_figs_png';
KZ_figs_pdf_dir = 'KZP_figs_pdf';

% (1) General Smoothing and filtering options:
KZP_parameters.smoothing_option = 1; % set = 1 to smooth long profile, 0 for no smoothing
KZP_parameters.smoothing_window = 151;  %Size of smoothing window for Savitzky-Golay filter in cells. Use 101-151 (default: 151) for 1m DEM, for 10m DEM use 15, for 30m DEM use 11 (cells)
KZP_parameters.sgolayfilt_order = 11;  % polynomial order of Savitzky-Golay filter, must be smaller than smoothing_window
KZP_parameters.min_trib_size = 100;     % number of cells or minimum tributary length considered
% number of grid cells (Savitzky-Golay smoothing, consider spatial resolution.
% Suggested values are 25 to 100 gridcells for 1m DEM and values of 10 to 20 for 30m DEMs
KZP_parameters.create_long_prof_4_tributaries = 1; % =1 for yes do this
KZP_parameters.kp_plot_size = 0.5; % (may change depending on size of largest knickzones)


% (2) Combine closely spaced knickzones option (1 must = 1, other must = 0)
KZP_parameters.chi_lump_option = 1; % use units of chi (essentially steepness)
KZP_parameters.lumping_distance_upstream = 75; % distance upstream (combines closely spaced knickpoints)
%combines knickzones spaced closer than this streamwise distance (m). Value
%should reflect spacing of waterfalls or steps in individual knickzones 
KZP_parameters.lumping_search_distance = 50;  % chi distance (combines closely spaced knickpoints)
KZP_parameters.distance_upstream_lump_option = 0; % use units of distance upstream (m) (more tangible)

% (3) Filter small knickzones option
KZP_parameters.kp_magnitude_filter_option = 1;
% uses knickzone magnitude (m) as measurement to remove small knickzones
%(measures steepness (Ksn) anomaly: the elevation drop across knickzone 
% after normalizing for the average stream ksn)
KZP_parameters.kp_relief_filter_option = 0; % use a minimum knickzone relief (m) (elevation drop across knickzone)

% NOTE:  knickzone relief systematically increases upstream as background 
% stream gradients increase upstream (why kp_magnitude is prefered for
% filtering)

KZP_parameters.min_kp_size1 = 1.1; % minimum knickpoint size (magnitude in m) before combining closely spaced knickpoints
KZP_parameters.min_kp_size2 = 3; % minimum knickpoint size after combining closely spaced knickpoints
KZP_parameters.min_kp_slope = 0.001;  % minimum steepness anomoly of knickpoint (ksn greater than expected ksn)
KZP_parameters.min_kp_size2_magnitude = 5; % minimum knickpoint magnitude before combining closely spaced knickpoints
KZP_parameters.min_kp_size2_relief = 10; % minimum knickpoint magnitude before combining closely spaced knickpoints

% (4) Additional Parameters
KZP_parameters.theta_bf_option = 1; % set = 1 if you want to fit theta to co-linearize 
% (tributaries for chi plots) Set = 0 if you want to use ref. concavity
KZP_parameters.theta_ref = KZP_parameters.theta; % ref. concavity (usually -0.45)

%% KZP: Calibration parameters
KZP_parameters.Calibration_option = 0; 
% set to 1 for calibration (default no calibration)

KZP_parameters.do_you_have_calibration_KZ_bases = 0; 
% Set = 0 if you only have positions and heights of knickzone lips in your
% calibration dataset. Set = 1 if you have lips and bases

% If calibrating alogrithm parameters, specify the name in front of the
% .csv calibration files that contain the knickzone positions and reliefs
KZP_parameters.KZ_lips_calib_fname = 'Smug_lips_1m.csv'; % PAGE 5 (4c. in manual)
KZP_parameters.KZ_bases_calib_fname = 'Smug_bases_1m.csv';
% these are tables containing the position of the calibration knickzones
% function csvread will read these files in starting at row 2 (assuming row
% 1 is a header)

KZP_parameters.KZ_calib_easting_column_num = 10; 
% column number in calibration knickzone .csv file containin easting
% coordinates. 10 is default for tables constructed from profiler toolbar

KZP_parameters.KZ_calib_northing_column_num = 11; 
% colmun number in calibration knickzone .csv file containin northing
% coordinates. 11 is default for tables constructed from profiler toolbar

KZP_parameters.KZ_calib_relief_column_num = 13; 
% column number in calibration knickzone .csv file containin relief
% coordinates. 13 is default for tables constructed from profiler toolbar
% and following instruction manual. If you only have KZ lips, you don't
% have to include relief measurements in your calibration file.

KZP_parameters.calibration_snapping_tolerance = 50; % (m) % 30m for 30m DEMs
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
KZP_parameters.error_radius = 50;
% spatial tolerance used to determine true positives, false positives, 
% & false negatives. May increase with coarser DEM resolution.


% Input parameter values to cycle through. Enter same number of values for each parameter
% Probably use < 8 different values for each parameter or else may take a
% while to run.

KZP_parameters.SG_smoothing_calib = [11 51 101 151]; % sgolay smoothing window size (cells)
% Input the range of smoothing windows you want to loop through
% ^^ will change with different DEM resolutions 
%(decrease with lower resolution DEMs)

% Input the range of knickzone lumping windows that you want to loop through
KZP_parameters.lumping_search_distance_calib = [25 50 75 100]; % (m) upstream
% if two knickpoints are within this distance from one another they are 
% combined into one larger knickzone extending from the base of the
% downstream knickzone to the lip of the upstream knickzone

% set a minimum knickpoint size (pre-lumping) NOTE: the algorithm results 
% are usually not very sensitive to this parameter
KZP_parameters.min_kp_size1_calib = [0.25 0.5 0.75 1]; 
% (before running the closely spaced knickpoint combining function) 
% meters of elev drop not explained by stream gradient

% set a minimum knickpoint size for after the lumping function
KZP_parameters.min_kp_size2_calib = [5 10 15 20]; 
% minimum knickpoint size 2 (EITHER magntiude or relief depending on inputs
% in lines 36 and 41)
% (after running the closely spaced knickpoint combining function) 
% meters of elev drop not explained by stream gradient


%% Parameters and settings for GDAL
%Determinig what system we are running on
KZP_parameters.computer_system = computer;
if strcmp(KZP_parameters.computer_system, 'PCWIN64')
    if exist('C:\OSGeo4W64\bin\gdalsrsinfo.exe', 'file') == 2
        KZP_parameters.gdalsrsinfo_cmd = '!C:\OSGeo4W64\bin\gdalsrsinfo.exe';
    elseif exist('C:\OSGeo4W\bin\gdalsrsinfo.exe', 'file') == 2
        KZP_parameters.gdalsrsinfo_cmd = '!C:\OSGeo4W\bin\gdalsrsinfo.exe';
    else
        fprintf(1,'Can not find ''gdalsrsinfo''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.gdalsrsinfo_cmd = '!C:\OSGeo4W64\bin\gdalsrsinfo.exe';
    end
    
    if exist('C:\OSGeo4W64\bin\ogr2ogr.exe', 'file') == 2
        KZP_parameters.ogr2ogr_cmd = '!C:\OSGeo4W64\bin\ogr2ogr.exe';
    elseif exist('C:\OSGeo4W\bin\ogr2ogr.exe', 'file') == 2
        KZP_parameters.ogr2ogr_cmd = '!C:\OSGeo4W\bin\ogr2ogr.exe';
    else
        fprintf(1,'Can not find ''ogr2ogr''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.ogr2ogr_cmd = '!C:\OSGeo4W64\bin\ogr2ogr.exe';
    end
    
    if exist('C:\OSGeo4W64\bin\gdaldem.exe', 'file') == 2
        KZP_parameters.gdaldem_cmd = '!C:\OSGeo4W64\bin\gdaldem.exe';
    elseif exist('C:\OSGeo4W\bin\gdaldem.exe', 'file') == 2
        KZP_parameters.gdaldem_cmd = '!C:\OSGeo4W\bin\gdaldem.exe';
    else
        fprintf(1,'Can not find ''gdaldem''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.gdaldem_cmd = '!C:\OSGeo4W64\bin\gdaldem.exe';
    end

    if exist('C:\OSGeo4W64\bin\gdal_polygonize.py', 'file') == 2
        KZP_parameters.polygonize_cmd = '!C:\OSGeo4W64\bin\gdal_polygonize.py';
    elseif exist('C:\OSGeo4W\bin\gdal_polygonize.py', 'file') == 2
        KZP_parameters.polygonize_cmd = '!C:\OSGeo4W\bin\gdal_polygonize.py';
    else
        fprintf(1,'Can not find ''gdal_polygonize.py''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.polygonize_cmd = '!C:\OSGeo4W64\bin\gdal_polygonize.py';
    end
    KZP_parameters.remove_cmd = '!del';
    KZP_parameters.mv_cmd = '!move';
	KZP_parameters.dir_sep = '\';
end

if strcmp(KZP_parameters.computer_system, 'GLNXA64')
    if exist('/usr/bin/gdalsrsinfo', 'file') == 2
        KZP_parameters.gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /usr/bin/gdalsrsinfo';
    elseif exist('/usr/local/bin/gdalsrsinfo', 'file') == 2
        KZP_parameters.gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdalsrsinfo';
    else
        fprintf(1,'Can not find ''gdalsrsinfo''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH='''' /usr/bin/gdalsrsinfo';
    end
    
    if exist('/usr/bin/gdaldem', 'file') == 2
        KZP_parameters.gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /usr/bin/gdaldem';
    elseif exist('/usr/local/bin/gdaldem.exe', 'file') == 2
        KZP_parameters.gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdaldem';
    else
        fprintf(1,'Can not find ''gdaldem''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /usr/bin/gdaldem';
    end
    
    if exist('/usr/bin/ogr2ogr', 'file') == 2
        KZP_parameters.ogr2ogr_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /usr/bin/ogr2ogr';
    elseif exist('/usr/local/bin/ogr2ogr', 'file') == 2
        KZP_parameters.ogr2ogr_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/ogr2ogr';
    else
        fprintf(1,'Can not find ''ogr2ogr''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.ogr2ogr_cmd = '!env LD_LIBRARY_PATH='''' /usr/bin/ogr2ogr';
    end
    
    if exist('/usr/bin/gdal_polygonize.py', 'file') == 2
        KZP_parameters.polygonize_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /usr/bin/gdal_polygonize.py';
    elseif exist('/usr/local/bin/gdal_polygonize.py', 'file') == 2
        KZP_parameters.polygonize_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdal_polygonize.py';
    else
        fprintf(1,'Can not find ''gdal_polygonize.py''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.polygonize_cmd = '!env LD_LIBRARY_PATH='''' /usr/bin/gdal_polygonize.py';
    end
    KZP_parameters.remove_cmd = '!rm';
	KZP_parameters.mv_cmd = '!mv';
    KZP_parameters.dir_sep = '/';
end
   
if strcmp(KZP_parameters.computer_system, 'MACI64')
    if exist('/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdalsrsinfo', 'file') == 2
        KZP_parameters.gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdalsrsinfo';
    elseif exist('/usr/local/bin/gdalsrsinfo', 'file') == 2
        KZP_parameters.gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdalsrsinfo';
    else
        fprintf(1,'Can not find ''gdalsrsinfo''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdalsrsinfo';
    end
    
    if exist('/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdaldem', 'file') == 2
        KZP_parameters.gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdaldem';
    elseif exist('/usr/local/bin/gdaldem', 'file') == 2
        KZP_parameters.gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdaldem';
    else
        fprintf(1,'Can not find ''gdaldem''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdaldem';
    end
    
    if exist('/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/ogr2ogr', 'file') == 2
        KZP_parameters.ogr2ogr_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /Library/Frameworks/GDAL.framework/Versions/1.11/Programs/ogr2ogr';
    elseif exist('/usr/local/bin/ogr2ogr', 'file') == 2
        KZP_parameters.ogr2ogr_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/ogr2ogr';
    else
        fprintf(1,'Can not find ''ogr2ogr''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.ogr2ogr_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/ogr2ogr';
    end
    
    if exist('/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdal_polygonize.py', 'file') == 2
        KZP_parameters.polygonize_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdal_polygonize.py';
    elseif exist('/usr/local/bin/gdal_polygonize.py', 'file') == 2
        polygonize_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdal_polygonize.py';
    else
        fprintf(1,'Can not find ''gdal_polygonize.py''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        KZP_parameters.polygonize_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdal_polygonize.py';
    end
    KZP_parameters.remove_cmd = '!rm';
	KZP_parameters.mv_cmd = '!mv';
    KZP_parameters.dir_sep = '/';
end


%% Usually no changes necessary below this line

KZP_parameters.max_knickpointsize2plot = 30;

[v d] = version;
if  str2num(v(end-5:end-2)) < 2014
    %running old Matlab
    KZP_parameters.MATLABV = 0;
elseif str2num(v(end-5:end-2)) >= 2014
    KZP_parameters.MATLABV = 1;
end
clear v d

%strsplit only works on newer Matlab version
if KZP_parameters.MATLABV == 1
    foo = strsplit(KZP_parameters.DEM_basename, '/');
    KZP_parameters.DEM_basename_nodir  = foo{end};
else
    foo = textscan(KZP_parameters.DEM_basename,'%s');
    KZP_parameters.DEM_basename_nodir  = char(foo{end});
end
clear foo
KZP_parameters.TIF_DIR_basename = strcat('DEM_geotiff/', KZP_parameters.DEM_basename);
KZP_parameters.SHP_DIR_basename = strcat('DEM_shapefiles/', KZP_parameters.DEM_basename);

KZP_parameters.map_dirname = 'DEM_maps';
KZP_parameters.plots_dirname = 'DEM_plots';
KZP_parameters.shapefile_dirname = 'DEM_shapefiles';
KZP_parameters.KZP_shapefile_dirname = 'KZP_shapefiles';
KZP_parameters.geotiff_dirname = 'DEM_geotiff';
KZP_parameters.KZP_csv_dirname = 'KZP_csv';
KZP_parameters.KZP_plots_dirname = 'KZP_plots';

%filtering options for diffusion filtering
KZP_parameters.difkernelWidth = 5;
KZP_parameters.difSSquared = 0.05;
KZP_parameters.difFilterI = 10;
KZP_parameters.difMethod = 'PeronaMalik2';
KZP_parameters.difTimeIncrement = 0.02;

KZP_parameters.DEM_fname = strcat(KZP_parameters.DEM_basename, '.tif');  % filename of DEM
KZP_parameters.DEM_MAT_fname = strcat(KZP_parameters.DEM_basename, '.mat');  % filename for DEM in Matlab MAT file
KZP_parameters.DEM_HYD_MAT_fname = strcat(KZP_parameters.DEM_basename, '_HYD.mat');  % filename for hydrologically corrected DEM and derivatives
KZP_parameters.DEM_STR_MAT_fname = strcat(KZP_parameters.DEM_basename, '_STR.mat');  % filename for streamobjs in vector format
KZP_parameters.DEM_FIL_fname = strcat(KZP_parameters.TIF_DIR_basename, '_FIL.tif');  % filename of filled DEM
KZP_parameters.DEM_FAC_fname = strcat(KZP_parameters.TIF_DIR_basename, '_FAC.tif');  % filename of FlowAccumulation Grid
KZP_parameters.DEM_dbasin_fname = strcat(KZP_parameters.TIF_DIR_basename, '_DBASIN.tif');  % filename of DrainageBasin Grid
foo = sprintf('_rel_%d_m.tif', KZP_parameters.relief_values_m(1));
KZP_parameters.DEM_rel_1_fname = strcat(KZP_parameters.TIF_DIR_basename, foo);  % filename of Relief_1 Grid
foo = sprintf('_rel_%d_m.tif', KZP_parameters.relief_values_m(2));
KZP_parameters.DEM_rel_2_fname = strcat(KZP_parameters.TIF_DIR_basename, foo);  % filename of Relief_2 Grid
foo = sprintf('_rel_%d_m.tif', KZP_parameters.relief_values_m(3));
KZP_parameters.DEM_rel_3_fname = strcat(KZP_parameters.TIF_DIR_basename, foo);  % filename of Relief_3 Grid
clear foo
KZP_parameters.AOI_DEM_curv_profc_fname = strcat(KZP_parameters.TIF_DIR_basename, '_curv_profc.tif');  % filename of profile curvature
KZP_parameters.AOI_DEM_curv_planc_fname = strcat(KZP_parameters.TIF_DIR_basename, '_curv_planc.tif');  % filename of planform curvature
KZP_parameters.AOI_DEM_curv_meanc_fname = strcat(KZP_parameters.TIF_DIR_basename, '_curv_meanc.tif');  % filename of planform curvature
KZP_parameters.AOI_DEM_diff_fname = strcat(KZP_parameters.TIF_DIR_basename, '_diffusionf.tif');  % filename of diffusion filtered DEM
KZP_parameters.AOI_DEM_wienerf_fname = strcat(KZP_parameters.TIF_DIR_basename, '_wienerf.tif');  % filename of Wiener filtered DEM
KZP_parameters.DEM_gradient8_fname = strcat(KZP_parameters.TIF_DIR_basename, '_gradient8.tif');  % filename of Gradient of DEM
KZP_parameters.DEM_SSP_fname = strcat(KZP_parameters.TIF_DIR_basename, '_SSP.tif');  % filename of Gradient of DEM
KZP_parameters.DEM_ksn045_fname = strcat(KZP_parameters.TIF_DIR_basename, '_Ksn045.tif');  % filename of Gradient of DEM
KZP_parameters.DEM_ks_adj_fname = strcat(KZP_parameters.TIF_DIR_basename, '_Ks_adj.tif');  % filename of Gradient of DEM
KZP_parameters.AOI_dbasins_stats_fname = strcat(KZP_parameters.TIF_DIR_basename, '_dbasins_stats.tif');  % filename of Gradient of DEM
KZP_parameters.AOI_dbasins_stats_vector_fname = strcat(KZP_parameters.SHP_DIR_basename, '_dbasins_stats_shp.shp');  % filename of Gradient of DEM
KZP_parameters.AOI_dbasins_stats_CNTR_csv_fname = strcat(KZP_parameters.DEM_basename, '_db_stats_CNTR.csv');  % filename of Gradient of DEM
KZP_parameters.AOI_dbasins_stats_OUT_csv_fname = strcat(KZP_parameters.DEM_basename, '_db_stats_OUT.csv');  % filename of Gradient of DEM

%set quality flag for output to PDF/JPG
KZP_parameters.quality_flag = '-q75';

save -v7.3 DEM_MAT_parameters.mat KZP_parameters
