%% (1) Load parameter file and prepare DEM

if (exist('DEM_MAT_parameters.mat', 'file') == 2)
    load DEM_MAT_parameters.mat
else
    fprintf('\n File: DEM_MAT_parameters.mat does not exist.\nGenerate Parameter file with KZP_parameters_1m.m first!');
    return
end
warning('off');

% Create Directories:
if (exist(KZP_parameters.map_dirname, 'dir') ~= 7)
    mkdir(KZP_parameters.map_dirname)
end

if (exist(KZP_parameters.plots_dirname, 'dir') ~= 7)
    mkdir(KZP_parameters.plots_dirname)
end

if (exist(KZP_parameters.shapefile_dirname, 'dir') ~= 7)
    mkdir(KZP_parameters.shapefile_dirname)
end
if (exist(KZP_parameters.geotiff_dirname, 'dir') ~= 7)
    mkdir(KZP_parameters.geotiff_dirname)
end
if (exist([KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep, 'trunk'], 'dir') ~= 7)
    mkdir([KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep, 'trunk'])
end

if (exist([KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep, 'STO'], 'dir') ~= 7)
    mkdir([KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep, 'STO'])
end

fprintf(1,'KZP topometrics step 1 of 4: loading DEM: %s\n', KZP_parameters.DEM_fname);
if exist(KZP_parameters.DEM_MAT_fname, 'file') == 0 || KZP_parameters.REGEN == 1
    [AOI_DEM, AOI_x, AOI_y] = DEM_load(KZP_parameters.DEM_fname, KZP_parameters.DEM_MAT_fname);
elseif exist(KZP_parameters.DEM_MAT_fname, 'file') == 2
    if exist('AOI_DEM', 'var') ~= 1
        load(KZP_parameters.DEM_MAT_fname)
    end
end
