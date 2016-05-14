function [AOI_dbasins_regionprops, AOI_dbasins, AOI_dbasins_outlet] = ...
    BSN_calcstats(AOI_FD, min_drainage_area_to_process, DEM_dbasin_fname, REGEN)
%calculate basin statistics

fprintf(1,'\tcalculating basin statistics\n');

if exist('AOI_dbasins', 'var') ~= 2
    [AOI_dbasins, AOI_dbasins_outlet] = drainagebasins(AOI_FD);
end

AOI_dbasins_regionprops = regionprops(...
    AOI_dbasins.Z,'Area','Centroid', 'PixelIdxList'); %region stored in number of pixels
AOI_dbasins_idx = find([AOI_dbasins_regionprops.Area].*...
    (AOI_dbasins.cellsize.^2) < min_drainage_area_to_process);
AOI_dbasins_2remove = AOI_dbasins_regionprops(AOI_dbasins_idx);
clear AOI_dbasins_2remove_idx
for j = 1:length(AOI_dbasins_2remove)
    if j == 1
        AOI_dbasins_2remove_idx(1:...
            length(AOI_dbasins_2remove(j).PixelIdxList)) = ...
            AOI_dbasins_2remove(j).PixelIdxList;
    else
        AOI_dbasins_2remove_idx(numel(AOI_dbasins_2remove_idx)+1:...
            numel(AOI_dbasins_2remove_idx)+...
            length(AOI_dbasins_2remove(j).PixelIdxList)) = ...
            AOI_dbasins_2remove(j).PixelIdxList;
    end
end
AOI_dbasins_regionprops(AOI_dbasins_idx) = [];
AOI_dbasins_outlet(AOI_dbasins_idx) = [];
%remove small basins from AOI_dbasins
AOI_dbasins.Z(AOI_dbasins_2remove_idx) = 0;
% HERE could add regional statistics for each catchment

if exist(DEM_dbasin_fname, 'file') ~= 2 || REGEN == 1
    GRIDobj2geotiff(AOI_dbasins,DEM_dbasin_fname);
end
