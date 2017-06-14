%% GUI for KZP (KnickZone Picker)
% Bodo Bookhagen (bodo.bookhagen@uni-potsdam.de), May 2017
%
% 
% KZ Picker (KZP) Matlab code by Bodo Bookhagen
% (bodo.bookhagen@uni-potsdam.de) and Al Neely (abn5031@psu.edu).
%
% loads standard parameters - modified version will be written/saved
if exist('DEM_MAT_parameters.mat') == 2
    load DEM_MAT_parameters.mat
else
    KZP_parameters_1m_example
end

%Generate figures with parameters read from parameter file
%DEM_MAT_parameters.mat

KZP_parameters = str;  % The name of the variable the user wishes to have in base.
S.CNT = 0;  % The number of times user pressed the pushbutton.
S.CHC = [];  % Holds the strings which represent the operations performed.
S.fh = figure('units','pixels',...
              'position',[400 400 300 130],...
              'menubar','none',...
              'name','GUI_32',...
              'numbertitle','off',...
              'resize','off',...
              'deletefcn',{@fig_del,S});
COL = get(S.fh,'color');          
S.pp = uicontrol('style','pop',...
                  'unit','pix',...
                  'position',[10 20 120 30],...
                  'string',{'Add';'Multiply';'Subtract';'Divide';'Power'});
S.ed(1) = uicontrol('style','edit',...
                    'unit','pix',...
                    'position',[10 90 70 30],...
                    'string','3');
S.tx(1) = uicontrol('style','text',...
                    'unit','pix',...
                    'position',[85 90 20 30],...
                    'string','+',...
                    'fontsize',16,...
                    'backgroundcolor',COL);                  
S.ed(2) = uicontrol('style','edit',...
                    'unit','pix',...
                    'position',[110 90 70 30],...
                    'string','2');  
S.tx(2) = uicontrol('style','text',...
                    'unit','pix',...
                    'position',[185 90 20 30],...
                    'string','=',...
                    'fontsize',16,...
                    'backgroundcolor',COL);                 
S.ed(3) = uicontrol('style','edit',...
                    'unit','pix',...
                    'position',[220 90 70 30],...
                    'string','answer');
S.pb = uicontrol('style','push',...
                  'unit','pix',...
                  'position',[160 20 120 30],...
                  'string','Calculate');
set([S.pp,S.pb],'callback',{@pb_call,S});               


KZP_topometrics
KZP_knickzone_processing

% Standard parameters to set (from KZP_parameters_1m)

[FileName,PathName] = uigetfile({'*.tif';'*.mat'},'Select the DEM (.tif or .mat) file');
DEM_basename = FileName(1:end-4);



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

% If this is set to 1, relief, curvature, DEM filtering, and
% hillshade DEMs are being calculated. If set to 0, no additional
% calculations are being performed.
RELIEF_CURVATURE = 1; 

%If this is set to 1 ridgecrest maps are being generated. This is
%experimental and off by default (=0).
RIDGECREST = 0;

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
