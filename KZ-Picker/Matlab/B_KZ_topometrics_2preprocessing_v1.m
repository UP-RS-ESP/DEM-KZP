%% (2) Calculating flow direction, flow accumulation, and deliniate drainage basins
%
fprintf(1,'step 2 of 4: pre-processing DEM: %s\n', DEM_fname);

    % Pre-process DEM (e.g., calculate flow direction, flow accumulation,
    % and stream networks)
[AOI_FIL, AOI_FD, AOI_mg] = ...
    B1_DEM_preprocess(AOI_DEM, min_str_gradient);

% Deliniate Drainage basins
    AOI_dbasins = drainagebasins(AOI_FD);
