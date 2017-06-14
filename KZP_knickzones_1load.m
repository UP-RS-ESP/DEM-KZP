%% (1) Load DEM file and other data from preprocessing
%
if (exist('DEM_MAT_parameters.mat', 'file') == 2)
    load DEM_MAT_parameters.mat
else
    fprintf('\n File: DEM_MAT_parameters.mat does not exist.\nGenerate Parameter file with KZP_parameters.m first!');
    return
end
warning('off');

% Generate DEMs
if (exist(KZP_parameters.KZP_csv_dirname, 'dir') ~= 7)
    mkdir(KZP_parameters.KZP_csv_dirname)
end

if (exist(KZP_parameters.shapefile_dirname, 'dir') ~= 7)
    mkdir(KZP_parameters.shapefile_dirname)
end

if (exist(KZP_parameters.KZP_plots_dirname, 'dir') ~= 7)
    mkdir(KZP_parameters.KZP_plots_dirname)
end

if (exist(KZP_parameters.plots_dirname, 'dir') ~= 7)
    mkdir(KZP_parameters.plots_dirname)
end


if exist(KZP_parameters.DEM_MAT_fname, 'file') == 2
    if exist('AOI_DEM', 'var') ~= 1
        load(KZP_parameters.DEM_MAT_fname)
    end
end

if exist(KZP_parameters.DEM_HYD_MAT_fname, 'file') == 2
    if exist('AOI_STR_MS', 'var') ~= 1
        load(KZP_parameters.DEM_HYD_MAT_fname, 'AOI_DEM_gradient8', 'AOI_rivers_STR', ...
            'AOI_STR_MS');
    end
end

fprintf(1,'KZP identifying knickzones step 1 of 4: loading DEM and stream data for %s\n', KZP_parameters.DEM_fname);
if exist(KZP_parameters.DEM_STR_MAT_fname, 'file') == 2
    if exist('AOI_STR_S_chiplot', 'var') ~= 1
        load(KZP_parameters.DEM_STR_MAT_fname, 'AOI_STR_S_trunk_chiplot', 'AOI_STR_S_chiplot', ...
            'AOI_STR_streams_dbasins_unique', 'AOI_STR_S_slopearea_dbasins', ...
            'AOI_STR_area_subset', 'AOI_STR_slope_subset', 'AOI_STR_S_slopearea_dbasins_trunk_adj',...
            'AOI_STR_area_trunk_subset', 'AOI_STR_S_slopearea_dbasins', 'AOI_STR_S_slopearea_dbasins_trunk', ...
            'AOI_STR_slope_trunk_subset', 'AOI_STR_dbasins_unique_subset', 'fitresult_logspace', ...
            'AOI_STR_all_streams_trunk', 'AOI_STR_all_streams_trunk', 'gof_logspace', 'ci_logspace', 'p_logspace_a', 'p_logspace_b', ...
            'AOI_STR_S_slopearea_dbasins_trunk_adj', 'AOI_STR_S_slopearea_dbasins_trunk_chitheta')
    end
else
    fprintf(1, 'Data file %s does not exist. Run preprocessing steps first.\n', KZP_parameters.DEM_STR_MAT_fname)
end

% Turn warning off
warning('off');

%number of basins to be analyzed. This is determined by parameter
%KZP_parameters.min_drainage_area_to_process and set in the parameter .m
%file.
number_of_basins = length(AOI_STR_S_chiplot);
