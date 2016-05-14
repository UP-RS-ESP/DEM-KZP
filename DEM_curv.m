function [AOI_DEM_curv_profc, AOI_DEM_curv_planc, AOI_DEM_curv_meanc] = ...
    DEM_curv(AOI_DEM, AOI_DEM_curv_profc_fname, AOI_DEM_curv_planc_fname, AOI_DEM_curv_meanc_fname, REGEN)
%using topotoolbox' curvature

fprintf(1,'\tcalculating planform, profile, and mean curvature\n');
    if exist('AOI_DEM_curv_profc', 'var') ~= 1
        [AOI_DEM_curv_profc] = curvature(AOI_DEM, 'profc');
        [AOI_DEM_curv_planc] = curvature(AOI_DEM, 'planc');
        [AOI_DEM_curv_meanc] = curvature(AOI_DEM, 'meanc');
    end
    if exist(AOI_DEM_curv_profc_fname, 'file') ~= 2 || REGEN == 1
        GRIDobj2geotiff(AOI_DEM_curv_profc,AOI_DEM_curv_profc_fname);
        GRIDobj2geotiff(AOI_DEM_curv_planc,AOI_DEM_curv_planc_fname);
        GRIDobj2geotiff(AOI_DEM_curv_meanc,AOI_DEM_curv_meanc_fname);
    end
