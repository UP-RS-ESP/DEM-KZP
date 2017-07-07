%% Plot drainage basins deliniated from input DEM

% user will select the drainage basins they would like to analyze using the
% data cursor tool

% user inputs drainage basin index #s so the code can mask files

% plot drainage basins
figure(1)
imagesc(AOI_dbasins)

datacursormode on 
% user selects basins and inputs them
basin_index =input('Use the data cursor to select each drainage basin: enter index of each basin btwn brackets: ex. [1 2 3] ')

% basin_index used to identify and label different catchments