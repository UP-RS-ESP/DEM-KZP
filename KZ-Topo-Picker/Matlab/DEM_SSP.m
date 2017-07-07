function [AOI_DEM_gradient8, AOI_Q, AOI_SSP] = ...
    DEM_SSP(AOI_mg, DEM_gradient8_fname, AOI_FAC, DEM_SSP_fname, REGEN, MISC_FILES)
% Calculuate gradient 8 and specific stream power (SSP) from DEM

fprintf(1,'\tcalculating gradient8\n');
if exist('AOI_DEM_gradient8', 'var') ~= 1 || REGEN == 1
    AOI_DEM_gradient8 = gradient8(AOI_mg);
end
if exist(DEM_gradient8_fname, 'file') ~= 2 || REGEN == 1
    GRIDobj2geotiff(AOI_DEM_gradient8, DEM_gradient8_fname);
end

fprintf(1,'\tcalculating SSP\n');
% calculate Specific Stream Power, SSP
AOI_Q = AOI_FAC;
AOI_Q.Z = (AOI_FAC.Z .* ((AOI_FAC.cellsize^2) / (60*60*24*365)));
AOI_SSP = AOI_Q;
AOI_SSP = 9810 .* AOI_DEM_gradient8 .* AOI_Q ./ (AOI_Q .^ 0.5);
AOI_SSP.name = 'SSP';
AOI_SSP.zunit = 'W/m^2';

if MISC_FILES == 1
    if exist(DEM_SSP_fname, 'file') ~= 2 || REGEN == 1
        GRIDobj2geotiff(AOI_SSP, DEM_SSP_fname);
    end
end
