function basin_index = DEM_select_DB(AOI_dbasins)
%function basin_index = KZP_select_DB(AOI_dbasins)
%select drainage basin to process

figure, clf
imagesc(AOI_dbasins), colorbar

datacursormode on ;
% user selects basins and inputs them
basin_index = input('Use the data cursor to select each drainage basin: enter index of each basin btwn brackets: ex. [40 42 45] ');
