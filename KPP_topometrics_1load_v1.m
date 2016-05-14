%% (1) Load parameter file and prepare DEM

if (exist('DEM_MAT_parameters.mat', 'file') == 2)
    load DEM_MAT_parameters.mat
else
    fprintf('\n File: DEM_MAT_parameters.mat does not exist.\n');
end
warning('off');

% Pre-processing starts here:
if (exist('maps', 'dir') ~= 7)
    mkdir maps
end

if (exist('plots', 'dir') ~= 7)
    mkdir plots
end

if (exist('shapefiles', 'dir') ~= 7)
    mkdir shapefiles
end
if (exist('geotiff', 'dir') ~= 7)
    mkdir geotiff
end
if (exist('shapefiles/trunk', 'dir') ~= 7)
    mkdir shapefiles/trunk
end

if (exist('shapefiles/STO', 'dir') ~= 7)
    mkdir shapefiles/STO
end

fprintf(1,'step 1 of 4: loading DEM: %s\n', DEM_fname);
if exist(DEM_MAT_fname, 'file') == 0 || REGEN == 1
    [AOI_DEM, AOI_x, AOI_y] = DEM_load(DEM_fname, DEM_MAT_fname);
elseif exist(DEM_MAT_fname, 'file') == 2
    if exist('AOI_DEM', 'var') ~= 1
        load(DEM_MAT_fname)
    end
end
