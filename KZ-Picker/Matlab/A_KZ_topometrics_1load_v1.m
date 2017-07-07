%% (1) Load parameter file and prepare DEM

if (exist('DEM_MAT_parameters.mat', 'file') == 2)
    load DEM_MAT_parameters.mat
else
    fprintf('\n File: DEM_MAT_parameters.mat does not exist.\n');
end
warning('off');

mkdir SA_plots; % make directory to store output slope-area plots
    
if Calibration_option == 0; % if not performing calibration functions
   mkdir maps_and_chi_plots; % make directory to store maps and chiplots
    mkdir csv_KZ_databases; % make directory to store shapefiles with knickzone positions/information
end

% Pre-processing starts here:
fprintf(1,'step 1 of 4: loading DEM: %s\n', DEM_fname);
[AOI_DEM, AOI_x, AOI_y] = A1_DEM_load(DEM_fname,DEM_MAT_fname);

