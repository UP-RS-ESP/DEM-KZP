%% (1) Load DEM file and other data from preprocessing
%
if exist('DEM_MAT_parameters.mat', 'file') == 2
    load('DEM_MAT_parameters.mat')
end

if exist(DEM_MAT_fname, 'file') == 2
    if exist('AOI_DEM', 'var') ~= 1
        load(DEM_MAT_fname)
    end
end

if exist(DEM_HYD_MAT_fname, 'file') == 2
    if exist('AOI_STR_MS', 'var') ~= 1
        load(DEM_HYD_MAT_fname, 'AOI_DEM_gradient8', 'AOI_rivers_STR', ...
            'AOI_STR_MS');
    end
end

if exist(DEM_STR_MAT_fname, 'file') == 2
    if exist('AOI_STR_S_chiplot', 'var') ~= 1
        load(DEM_STR_MAT_fname, 'AOI_STR_S_trunk_chiplot', 'AOI_STR_S_chiplot', ...
            'AOI_STR_streams_dbasins_unique', 'AOI_STR_S_slopearea_dbasins', ...
            'AOI_STR_area_subset', 'AOI_STR_slope_subset', ...
            'AOI_STR_area_trunk_subset', ...
            'AOI_STR_slope_trunk_subset', 'AOI_STR_dbasins_unique_subset', ...
            'AOI_STR_all_streams_trunk', 'AOI_STR_all_streams_trunk', 'gof_logspace')
    end
end

% Turn warning off
warning('off');

%number of basins to be analyzed
number_of_basins = length(AOI_STR_S_chiplot);

