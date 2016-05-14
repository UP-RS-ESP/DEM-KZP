function [AOI_rivers_STR_area1, AOI_rivers_STR_area2, AOI_rivers_STR, ...
    AOI_rivers_slope, AOI_rivers_area] = DEM_riverslopearea(AOI_FD, str_area1, str_area2, 
AOI_rivers_STR_area1 = STREAMobj(AOI_FD,'minarea',str_area1, ...
        'unit', 'mapunits');
    AOI_rivers_STR_area2 = STREAMobj(AOI_FD,'minarea',str_area2, ...
        'unit', 'mapunits');
    
    if exist('AOI_rivers_STR', 'var') ~= 2 || REGEN == 1
        AOI_rivers_w = AOI_FAC.*(AOI_FAC.cellsize.^2) > ...
            min_drainage_area_to_process;
        AOI_rivers_STR = STREAMobj(AOI_FD,AOI_rivers_w);
    end
    if exist('AOI_rivers_area', 'var') ~= 2 || REGEN == 1 
        % calculate slope and area from DEM
        AOI_rivers_slope = gradient(AOI_rivers_STR,AOI_mg,'unit','tangent');
        AOI_rivers_area = AOI_FAC.Z(AOI_rivers_STR.IXgrid).*...
            (AOI_FAC.cellsize).^2;
    end
