function [AOI_ksn045, AOI_ksnadj, AOI_slopearea] = ...
    DEM_ksn(theta, AOI_FAC, AOI_DEM_gradient8, AOI_rivers_STR, AOI_mg, DEM_ksn045_fname, DEM_ksnadj_fname)
% Calcuate steepness index for entire DEM using theta = 0.45 and calculate
% steepness index with adjusted theta from region using log area-log slope
% plot

% calculate ksn with theta = 0.45
AOI_ksn045 = AOI_FAC;
AOI_ksn045 = AOI_DEM_gradient8 ./ (AOI_FAC.*...
    (AOI_FAC.cellsize^2)).^theta;
AOI_ksn045.name = 'ksn045';
AOI_ksn045.zunit = 'm^0.9';
%calculate area-slope relation for entire area (may not be useful)
fprintf(1,'\tcalculating slope-area relation for all streams\n');
AOI_slopearea = slopearea(AOI_rivers_STR,AOI_mg,...
    AOI_FAC, 'areabinlocs', 'median', 'gradaggfun', 'median', ...
    'streamgradient', 'robust', 'plot', false);

% calculate ksn with adjusted theta for this region
AOI_ksnadj = AOI_FAC;
AOI_ksnadj = AOI_DEM_gradient8 ./ (AOI_FAC.*...
    (AOI_FAC.cellsize^2)).^AOI_slopearea.theta;
AOI_ksnadj.name = 'ksn_adjusted';
AOI_ksnadj.zunit = 'm^0.9';

if exist(DEM_ksn045_fname, 'file') ~= 2
    GRIDobj2geotiff(AOI_ksn045, DEM_ksn045_fname);
end
if exist(DEM_ksnadj_fname, 'file') ~= 2
    GRIDobj2geotiff(AOI_ksnadj, DEM_ksnadj_fname);
end
