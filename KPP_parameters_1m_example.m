%% Parameter file for knickpoint_picker
% You only will need to change this file and adjust the parameters
DEM_basename = 'Pozo_DTM_noveg_UTM11_NAD83'; 
% we assume that there exists a TIF file <DEM_basename>.tif
% no need to add .tif here
% Make sure that the TIF file is in UTM or other equal-area projected
% coordinates

DEM_basename_no_underscore = strrep(DEM_basename,'_','\_');

area_threshold = 5e3;
% set drainage area threshhold in m^2.  5e3 seems to work well for high-resolution 1m DEMs.
% You can investigate good values for this using a slope/area plot and determining
% on average at what drainage area a negative power law relationship starts

min_drainage_area_to_process = 1e5;
% identifies min. size of drainage area to be processed. No catchments
% smaller than this size will be included in the catchment analysis. Useful
% values are 1e6 or 1e7 (depending on size and resolution of DEM). This
% often can be viewed as the drainage area threshold between hillslope and
% fluvial processes.

min_dbasins_stats_to_process = 1e6;
% minimum drainage basin area for which to calculate statistics. For all
% separate drainage basins larger than this value, slope-area and chi plots
% are being calculated and generated. This is valid only for separate
% drainage basins, for example a DEM from an island with separate drainage
% basins draining into the ocean. Useful values are 1e7 (10km2) or 1e8 (100
% km2) or 1e9 (1000km2)

stream_order = [2 3];
% Stream order for which to create individual basin plots. All basins equal
% this streamorder will be processed for topometric statistics. A stream
% order of 3 for a 30-m DEM includes catchments between 7-300 km2 (mean 45), 
% a stream order of 4 (30-m DEM): 28-1100 km2 (mean: 200 km2). You can also
% give the stream order in vector format. If you give multiple stream
% orders, each number in the vector will generate a separate Shapefile with
% drainage basins and statistics.

relief_values_m = [10 25 50];
% Local relief (or radius relief). This vector lists the search radius for
% relief. A value of 1000 indices a 1-km-radius relief calculation. 3
% calculations are performed. Please note that setting these to higher
% values will take a lot of time to calculate/process.
% If you leave this empty (e.g., []), no relief is being calculated

str_area1 = 1e6; % Lower drainage area to calculate streamnetwork (STREAMobj) from and generate shapefile
str_area2 = 1e6; % Upper drainage area to calculate streamnetwork (STREAMobj) from and generate shapefile
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
min_max_DA_fits = [1e5 1e8 100];

segL = 10; % segment length in meters for smoothing shapefile output, often spatial resolution * 10

min_str_gradient = 0.001;
% set minimum river gradient. This value will be forced on all river
% gradients less than min_str_gradient. For example, all 0-gradient
% sections of a river will be set to this value.

% If this is set to 1, relief, curvature, ridgecrest profiles, and
% hillshade DEMs are being calculated. If set to 0, no additional
% calculations are being performed.
MISC_FILES = 1; 

REGEN = 0; % regenerate all files (1 = Yes, 0 = No)

show_figs = 1; %Generate figures? (1 = Yes, 0 = No, 
% 2 = create figures, but no processing (assumes that all MAT files were
% already created - useful for updating figures)
% This will generate figures as PDFs and JPGs in the subfolder maps
% You will need to have installed export_fig from the Mathworks repository
% (free of charge).

PaperType_size = 'A4'; % Sets paper type and size for export_fit: 'A4' or 'letter'

% Parameters for knickpoint picker
min_trib_size = 100;     % cells (minimum tributary length we want to consider)
% number of grid cells (Savitzky-Golay smoothing, consider spatial resolution.
% Suggested values are 25 to 100 gridcells for 1-m data and values of 10 to
% 20 for 30-m data, must be odd!)
smoothing_window = 201;  
sgolayfilt_order = 11;  % polynomial order of Savitzky-Golay filter, must be smaller than smoothing_window
lumping_search_distance = 125;  % chi distance (combines closely spaced knickpoints
min_kp_size1 = 1.5; % minimum knickpoint size before combining closely spaced knickpoints
min_kp_size2 = 3; % minimum knickpoint size after combining closely spaced knickpoints
min_kp_slope = 0.001;  % minimum steepness anomoly of knickpoint (ksn greater than expected ksn)

% Settings for GDAL:
%Determinig what system we are running on
computer_system = computer;
if strcmp(computer_system, 'PCWIN64')
    if exist('C:\OSGeo4W64\bin\gdalsrsinfo.exe', 'file') == 2
        gdalsrsinfo_cmd = '!C:\OSGeo4W64\bin\gdalsrsinfo.exe';
    elseif exist('C:\OSGeo4W\bin\gdalsrsinfo.exe', 'file') == 2
        gdalsrsinfo_cmd = '!C:\OSGeo4W\bin\gdalsrsinfo.exe';
    else
        fprintf(1,'Can not find ''gdalsrsinfo''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        gdalsrsinfo_cmd = '!C:\OSGeo4W64\bin\gdalsrsinfo.exe';
    end
    
    if exist('C:\OSGeo4W64\bin\ogr2ogr.exe', 'file') == 2
        ogr2ogr_cmd = '!C:\OSGeo4W64\bin\ogr2ogr.exe';
    elseif exist('C:\OSGeo4W\bin\ogr2ogr.exe', 'file') == 2
        ogr2ogr_cmd = '!C:\OSGeo4W\bin\ogr2ogr.exe';
    else
        fprintf(1,'Can not find ''ogr2ogr''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        ogr2ogr_cmd = '!C:\OSGeo4W64\bin\ogr2ogr.exe';
    end
    
    if exist('C:\OSGeo4W64\bin\gdaldem.exe', 'file') == 2
        gdaldem_cmd = '!C:\OSGeo4W64\bin\gdaldem.exe';
    elseif exist('C:\OSGeo4W\bin\gdaldem.exe', 'file') == 2
        gdaldem_cmd = '!C:\OSGeo4W\bin\gdaldem.exe';
    else
        fprintf(1,'Can not find ''gdaldem''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        gdaldem_cmd = '!C:\OSGeo4W64\bin\gdaldem.exe';
    end

    if exist('C:\OSGeo4W64\bin\gdal_polygonize.py', 'file') == 2
        polygonize_cmd = '!C:\OSGeo4W64\bin\gdal_polygonize.py';
    elseif exist('C:\OSGeo4W\bin\gdal_polygonize.py', 'file') == 2
        polygonize_cmd = '!C:\OSGeo4W\bin\gdal_polygonize.py';
    else
        fprintf(1,'Can not find ''gdal_polygonize.py''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        polygonize_cmd = '!C:\OSGeo4W64\bin\gdal_polygonize.py';
    end
    remove_cmd = '!del';
    mv_cmd = '!move';
	dir_sep = '\';
end

if strcmp(computer_system, 'GLNXA64')
    if exist('/usr/bin/gdalsrsinfo', 'file') == 2
        gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /usr/bin/gdalsrsinfo';
    elseif exist('/usr/local/bin/gdalsrsinfo', 'file') == 2
        gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdalsrsinfo';
    else
        fprintf(1,'Can not find ''gdalsrsinfo''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH='''' /usr/bin/gdalsrsinfo';
    end
    
    if exist('/usr/bin/gdaldem', 'file') == 2
        gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /usr/bin/gdaldem';
    elseif exist('/usr/local/bin/gdaldem.exe', 'file') == 2
        gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdaldem';
    else
        fprintf(1,'Can not find ''gdaldem''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /usr/bin/gdaldem';
    end
    
    if exist('/usr/bin/ogr2ogr', 'file') == 2
        ogr2ogr_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /usr/bin/ogr2ogr';
    elseif exist('/usr/local/bin/ogr2ogr', 'file') == 2
        ogr2ogr_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/ogr2ogr';
    else
        fprintf(1,'Can not find ''ogr2ogr''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        ogr2ogr_cmd = '!env LD_LIBRARY_PATH='''' /usr/bin/ogr2ogr';
    end
    
    if exist('/usr/bin/gdal_polygonize.py', 'file') == 2
        polygonize_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /usr/bin/gdal_polygonize.py';
    elseif exist('/usr/local/bin/gdal_polygonize.py', 'file') == 2
        polygonize_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdal_polygonize.py';
    else
        fprintf(1,'Can not find ''gdal_polygonize.py''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        polygonize_cmd = '!env LD_LIBRARY_PATH='''' /usr/bin/gdal_polygonize.py';
    end
    remove_cmd = '!rm';
	mv_cmd = '!mv';
    dir_sep = '/';
end
   
if strcmp(computer_system, 'MACI64')
    if exist('/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdalsrsinfo', 'file') == 2
        gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdalsrsinfo';
    elseif exist('/usr/local/bin/gdalsrsinfo', 'file') == 2
        gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdalsrsinfo';
    else
        fprintf(1,'Can not find ''gdalsrsinfo''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        gdalsrsinfo_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdalsrsinfo';
    end
    
    if exist('/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdaldem', 'file') == 2
        gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdaldem';
    elseif exist('/usr/local/bin/gdaldem', 'file') == 2
        gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdaldem';
    else
        fprintf(1,'Can not find ''gdaldem''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        gdaldem_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdaldem';
    end
    
    if exist('/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/ogr2ogr', 'file') == 2
        ogr2ogr_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /Library/Frameworks/GDAL.framework/Versions/1.11/Programs/ogr2ogr';
    elseif exist('/usr/local/bin/ogr2ogr', 'file') == 2
        ogr2ogr_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/ogr2ogr';
    else
        fprintf(1,'Can not find ''ogr2ogr''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        ogr2ogr_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/ogr2ogr';
    end
    
    if exist('/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdal_polygonize.py', 'file') == 2
        polygonize_cmd = '!env LD_LIBRARY_PATH=''/usr/lib'' /Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdal_polygonize.py';
    elseif exist('/usr/local/bin/gdal_polygonize.py', 'file') == 2
        polygonize_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdal_polygonize.py';
    else
        fprintf(1,'Can not find ''gdal_polygonize.py''. Please add path manually in knickpoints_parameters.m file\n');
        %SET the PATH to the gdalsrsinfo_cmd in the next line
        polygonize_cmd = '!env LD_LIBRARY_PATH='''' /usr/local/bin/gdal_polygonize.py';
    end
    remove_cmd = '!rm';
	mv_cmd = '!mv';
    dir_sep = '/';
end


%% Usually no changes necessary below this line
theta = -0.45; % reference concavity

max_knickpointsize2plot = 40;

[v d] = version;
if  str2num(v(end-5:end-2)) < 2014
    %running old Matlab
    MATLABV = 0;
elseif str2num(v(end-5:end-2)) >= 2014
    MATLABV = 1;
end
clear v d

%strsplit only works on newer Matlab version
if MATLABV == 1
    foo = strsplit(DEM_basename, '/');
    DEM_basename_nodir  = foo{end};
else
    foo = textscan(DEM_basename,'%s');
    DEM_basename_nodir  = char(foo{end});
end
clear foo
TIF_DIR_basename = strcat('geotiff/', DEM_basename);
SHP_DIR_basename = strcat('shapefiles/', DEM_basename);

%filtering options for diffusion filtering
difkernelWidth = 5;
difSSquared = 0.05;
difFilterI = 10;
difMethod = 'PeronaMalik2';
difTimeIncrement = 0.02;

DEM_fname = strcat(DEM_basename, '.tif');  % filename of DEM
DEM_MAT_fname = strcat(DEM_basename, '.mat');  % filename for DEM in Matlab MAT file
DEM_HYD_MAT_fname = strcat(DEM_basename, '_HYD.mat');  % filename for hydrologically corrected DEM and derivatives
DEM_STR_MAT_fname = strcat(DEM_basename, '_STR.mat');  % filename for streamobjs in vector format
DEM_FIL_fname = strcat(TIF_DIR_basename, '_FIL.tif');  % filename of filled DEM
DEM_FAC_fname = strcat(TIF_DIR_basename, '_FAC.tif');  % filename of FlowAccumulation Grid
DEM_dbasin_fname = strcat(TIF_DIR_basename, '_DBASIN.tif');  % filename of DrainageBasin Grid
foo = sprintf('_rel_%d_m.tif', relief_values_m(1));
DEM_rel_1_fname = strcat(TIF_DIR_basename, foo);  % filename of Relief_1 Grid
foo = sprintf('_rel_%d_m.tif', relief_values_m(2));
DEM_rel_2_fname = strcat(TIF_DIR_basename, foo);  % filename of Relief_2 Grid
foo = sprintf('_rel_%d_m.tif', relief_values_m(3));
DEM_rel_3_fname = strcat(TIF_DIR_basename, foo);  % filename of Relief_3 Grid
clear foo
AOI_DEM_curv_profc_fname = strcat(TIF_DIR_basename, '_curv_profc.tif');  % filename of profile curvature
AOI_DEM_curv_planc_fname = strcat(TIF_DIR_basename, '_curv_planc.tif');  % filename of planform curvature
AOI_DEM_curv_meanc_fname = strcat(TIF_DIR_basename, '_curv_meanc.tif');  % filename of planform curvature
AOI_DEM_diff_fname = strcat(TIF_DIR_basename, '_diffusionf.tif');  % filename of diffusion filtered DEM
AOI_DEM_wienerf_fname = strcat(TIF_DIR_basename, '_wienerf.tif');  % filename of Wiener filtered DEM
DEM_gradient8_fname = strcat(TIF_DIR_basename, '_gradient8.tif');  % filename of Gradient of DEM
DEM_SSP_fname = strcat(TIF_DIR_basename, '_SSP.tif');  % filename of Gradient of DEM
DEM_ksn045_fname = strcat(TIF_DIR_basename, '_Ksn045.tif');  % filename of Gradient of DEM
DEM_ks_adj_fname = strcat(TIF_DIR_basename, '_Ks_adj.tif');  % filename of Gradient of DEM
AOI_dbasins_stats_fname = strcat(TIF_DIR_basename, '_dbasins_stats.tif');  % filename of Gradient of DEM
AOI_dbasins_stats_vector_fname = strcat(SHP_DIR_basename, '_dbasins_stats_shp.shp');  % filename of Gradient of DEM
AOI_dbasins_stats_CNTR_csv_fname = strcat(DEM_basename, '_db_stats_CNTR.csv');  % filename of Gradient of DEM
AOI_dbasins_stats_OUT_csv_fname = strcat(DEM_basename, '_db_stats_OUT.csv');  % filename of Gradient of DEM

%set quality flag for output to PDF/JPG
quality_flag = '-q75';

save -v7.3 DEM_MAT_parameters.mat *fname PaperType_size show_figs segL ...
    relief_values_m min_drainage_area_to_process area_threshold ...
    DEM_basename DEM_basename_nodir *cmd min_trib_size smoothing_window ...
    lumping_search_distance min_kp_size1 min_kp_size2 min_kp_slope ...
    quality_flag DEM_basename_no_underscore DEM_FIL_fname DEM_FAC_fname ...
    DEM_dbasin_fname DEM_ks_adj_fname DEM_ksn045_fname DEM_SSP_fname ...
    DEM_gradient8_fname DEM_rel_*_fname str_area1 str_area2 ...
    min_str_gradient theta AOI_dbasins_stats_CNTR_csv_fname polygonize_cmd ...
    min_dbasins_stats_to_process AOI_dbasins_stats_OUT_csv_fname ...
    min_max_DA_fits stream_order AOI_DEM_curv_planc_fname ...
    AOI_DEM_curv_profc_fname remove_cmd AOI_DEM_curv_meanc_fname mv_cmd ...
    dir_sep sgolayfilt_order max_knickpointsize2plot MATLABV REGEN ...
    difkernelWidth difSSquared difFilterI difMethod difTimeIncrement ...
    AOI_DEM_diff_fname AOI_DEM_wienerf_fname TIF_DIR_basename ...
    SHP_DIR_basename MISC_FILES
