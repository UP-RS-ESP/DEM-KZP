function [AOI_FIL, AOI_FD, AOI_mg] = ...
    B1_DEM_preprocess(AOI_DEM, min_str_gradient)
%
% preprocess and hydrologically correct DEM
% fill sinks in DEM:

AOI_FIL = fillsinks(AOI_DEM);

% create flow direction (carve lets streams cut down through e.g.
% bridges):
AOI_FD = FLOWobj(AOI_FIL,'preprocess','carve');

% set a minimum gradient (no place in DEM has a slope of 0)
AOI_mg = imposemin(AOI_FD,AOI_DEM,min_str_gradient);





