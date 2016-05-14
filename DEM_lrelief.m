function [AOI_DEM_rel_1, AOI_DEM_rel_2, AOI_DEM_rel_3] = ...
    DEM_lrelief(AOI_DEM, relief_values_m, DEM_rel_1_fname, DEM_rel_2_fname, DEM_rel_3_fname, REGEN)

if ~isempty(relief_values_m)
    % calculate relief and gradient
    fprintf(1,'\tcalculating relief 1: %d m\n', relief_values_m(1));
    AOI_DEM_rel_1 = localtopography(AOI_DEM, relief_values_m(1));
    if exist(DEM_rel_1_fname, 'file') ~= 2 || REGEN == 1
        GRIDobj2geotiff(AOI_DEM_rel_1,DEM_rel_1_fname);
    end
    
    fprintf(1,'\tcalculating relief 2: %d m\n', relief_values_m(2));
    AOI_DEM_rel_2 = localtopography(AOI_DEM, relief_values_m(2));
    if exist(DEM_rel_2_fname, 'file') ~= 2 || REGEN == 1
        GRIDobj2geotiff(AOI_DEM_rel_2,DEM_rel_2_fname);
    end
    
    fprintf(1,'\tcalculating relief 3: %d m\n', relief_values_m(3));
    AOI_DEM_rel_3 = localtopography(AOI_DEM, relief_values_m(3));
    if exist(DEM_rel_3_fname, 'file') ~= 2 || REGEN == 1
        GRIDobj2geotiff(AOI_DEM_rel_3,DEM_rel_3_fname);
    end
end
