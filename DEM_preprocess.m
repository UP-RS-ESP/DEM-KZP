function [AOI_FIL, AOI_FD, AOI_FAC, AOI_STR_w, AOI_mg, AOI_rivers_STR, ...
    AOI_rivers_STR_area1, AOI_rivers_STR_area2, AOI_resolution, minApix, ...
    AOI_FAC_w, AOI_rivers_slope, AOI_rivers_area, AOI_rivers_w] = ...
    DEM_preprocess(AOI_DEM, REGEN, dir_sep, gdaldem_cmd, DEM_FIL_fname, ...
    DEM_FAC_fname, DEM_basename, area_threshold, min_str_gradient, ...
    str_area1, str_area2, min_drainage_area_to_process, geotiff_dirname, MISC_FILES)
%
% preprocess and hydrologically correct DEM
% fill sinks in DEM:
if exist('AOI_FIL', 'var') ~= 1 || REGEN == 1
    AOI_FIL = fillsinks(AOI_DEM);
end
if MISC_FILES == 1
    if exist(DEM_FIL_fname, 'file') ~= 2 || REGEN == 1
        GRIDobj2geotiff(AOI_FIL,DEM_FIL_fname);
    end
end
if MISC_FILES == 1
    fprintf(1,'\tgenerating hillshade image from DEM using gdaldem\n');
    if exist(sprintf('%s%s%s%s', geotiff_dirname, dir_sep, DEM_basename, '_hs.tif'), 'file') ~= 2 || REGEN == 1
        eval([gdaldem_cmd, ' hillshade ', sprintf('%s%s', DEM_basename, '.tif'), ' ', ...
            sprintf('%s%s%s%s', geotiff_dirname, dir_sep, DEM_basename, '_hs.tif')]);
    end
end

% create flow direction (carve lets streams cut down through e.g.
% bridges):
fprintf(1,'\tgenerating flow direction grid\n');
if exist('AOI_FD', 'var') ~= 1 || REGEN == 1
    AOI_FD = FLOWobj(AOI_FIL,'preprocess','carve');
end

% Generate flow accumulation array (may take time):
fprintf(1,'\tgenerating flow accumulation grid\n');
if exist('AOI_FAC', 'var') ~= 1 || REGEN == 1
    AOI_FAC = flowacc(AOI_FD);
end
if MISC_FILES == 1
    if exist(DEM_FAC_fname, 'file') ~= 2 || REGEN == 1
        GRIDobj2geotiff(AOI_FAC,DEM_FAC_fname);
    end
end
AOI_resolution = AOI_DEM.refmat(2,1);
% get the DEM resolution, used to convert area threshold to meters squared

minApix = area_threshold/(AOI_resolution*AOI_resolution);
minApix = ceil(minApix);  % convert area threshold to meters specified above

AOI_FAC_w = AOI_FAC > minApix;  % masks flow accum grid (above threshold)

AOI_rivers_w = AOI_FAC.*(AOI_FAC.cellsize.^2) > ...
    min_drainage_area_to_process;
AOI_STR_w = STREAMobj(AOI_FD,AOI_FAC_w);
AOI_rivers_STR = STREAMobj(AOI_FD,AOI_rivers_w);
% calculate slope and area from DEM

% set a minimum gradient (no place in DEM has a slope of 0)
AOI_mg = imposemin(AOI_FD,AOI_DEM,min_str_gradient);

AOI_rivers_slope = gradient(AOI_rivers_STR,AOI_mg,'unit','tangent');
AOI_rivers_area = AOI_FAC.Z(AOI_rivers_STR.IXgrid).*...
    (AOI_FAC.cellsize).^2;

AOI_rivers_STR_area1 = STREAMobj(AOI_FD,'minarea',str_area1, ...
    'unit', 'mapunits');
AOI_rivers_STR_area2 = STREAMobj(AOI_FD,'minarea',str_area2, ...
    'unit', 'mapunits');

