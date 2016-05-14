%% (2) Calculating flow direction, flow accumulation, relief and other topographic metrics
%
fprintf(1,'step 2 of 4: pre-processing DEM: %s\n', DEM_fname);
if exist(DEM_HYD_MAT_fname, 'file') == 0 || REGEN == 1
    % Pre-process DEM (e.g., calculate flow direction, flow accumulation,
    % slope, and stream networks)
    if exist('AOI_FAC', 'var') ~= 1 || REGEN == 1
        [AOI_FIL, AOI_FD, AOI_FAC, AOI_STR_w, AOI_mg, AOI_rivers_STR, ...
            AOI_rivers_STR_area1, AOI_rivers_STR_area2, AOI_resolution, ...
            minApix, AOI_FAC_w, AOI_rivers_slope, AOI_rivers_area, AOI_rivers_w] = ...
            DEM_preprocess(AOI_DEM, REGEN, dir_sep, gdaldem_cmd, DEM_FIL_fname, ...
            DEM_FAC_fname, DEM_basename, area_threshold, min_str_gradient, ...
            str_area1, str_area2, min_drainage_area_to_process, MISC_FILES);
    end
    
    % delineate drainage basins
    if exist('AOI_dbasins_regionprops', 'var') ~= 1 || REGEN == 1
        [AOI_dbasins_regionprops, AOI_dbasins, AOI_dbasins_outlet] = ...
            BSN_calcstats(AOI_FD, min_dbasins_stats_to_process, DEM_dbasin_fname, REGEN);
    end
    
    % Calculate relief
    if MISC_FILES == 1
        if exist('AOI_DEM_rel_3', 'var') ~= 1 || REGEN == 1
            [AOI_DEM_rel_1, AOI_DEM_rel_2, AOI_DEM_rel_3] = ...
                DEM_lrelief(AOI_DEM, relief_values_m, DEM_rel_1_fname, DEM_rel_2_fname, DEM_rel_3_fname, REGEN);
        end
    end
    
    % calculatute Curvature from unfiltered DEM
    if MISC_FILES == 1
        if exist('AOI_DEM_curv_meanc', 'var') ~= 1 || REGEN == 1
            [AOI_DEM_curv_profc, AOI_DEM_curv_planc, AOI_DEM_curv_meanc] = ...
                DEM_curv(AOI_DEM, AOI_DEM_curv_profc_fname, AOI_DEM_curv_planc_fname, AOI_DEM_curv_meanc_fname, REGEN);
        end
    end
    
    % diffusion-filtering of DEM as described in Passalacqua et al. (2010a,
    % 2010b) and Perona and Malik (1990) and calculate curvature
    if MISC_FILES == 1
        if exist('AOI_DEM_diffusionf', 'var') ~= 1 || REGEN == 1
            AOI_DEM_diffusionf = ...
                DEM_diffusionf(AOI_DEM, difkernelWidth, difSSquared, ...
                difTimeIncrement, difFilterI, difMethod, AOI_DEM_diff_fname, REGEN);
            AOI_DEM_diffusionf_curv = DEM_curv2(AOI_DEM_diffusionf, strcat(TIF_DIR_basename, '_diffusionf_curv.tif'));
        end
    end
    
    % Wiener filtering of DEM and calculate curvature
    if MISC_FILES == 1
        if exist('AOI_DEM_wiener_curv', 'var') ~= 1 || REGEN == 1
            AOI_DEM_wiener = filter(AOI_DEM, 'wiener', [5 5]);
            AOI_DEM_wiener_curv = DEM_curv2(AOI_DEM_wiener, strcat(TIF_DIR_basename, '_wienerf_curv.tif'));
        end
    end
    
    %Calculate stream gradient and SSP
    if exist('AOI_SSP', 'var') ~= 1 || REGEN == 1
        [AOI_DEM_gradient8, AOI_Q, AOI_SSP] = ...
            DEM_SSP(AOI_mg, DEM_gradient8_fname, AOI_FAC, DEM_SSP_fname, REGEN, MISC_FILES);
    end
    
    % extract ridgecrests
    if MISC_FILES == 1
        if exist('AOI_ridgecrest_MS', 'var') ~= 1 || REGEN == 1
            [AOI_ridgecrest_MS] = ...
                DEM_ridgecrest(AOI_dbasins_regionprops, AOI_dbasins, AOI_dbasins_outlet);
        end
    end
    
    %Get DEM and curvature values for ridgecrest
    if MISC_FILES == 1
        if exist('AOI_ridgecrest_DEM', 'var') ~= 1 || REGEN == 1
            [AOI_ridgecrest_DEM] = ...
                DEM_ridgevalues(AOI_ridgecrest_MS, AOI_DEM);
            [AOI_ridgecrest_diffusionf_curv] = ...
                DEM_ridgevalues(AOI_ridgecrest_MS, AOI_DEM_diffusionf_curv);
            [AOI_ridgecrest_wiener_curv] = ...
                DEM_ridgevalues(AOI_ridgecrest_MS, AOI_DEM_wiener_curv);
        end
    end
    
    % Analyze ridgecrests and get normalized profiles
    if MISC_FILES == 1
        if exist('AOI_ridgecrest_MS_deltay_all', 'var') ~= 1 || REGEN == 1
            [AOI_ridgecrest_MS, AOI_ridgecrest_MS_Dy_all, AOI_ridgecrest_MS_yi_all, ridgecrest_stepsize] = ...
                DEM_ridgecrest_analyze(AOI_ridgecrest_MS, AOI_ridgecrest_DEM, ...
                DEM_HYD_MAT_fname);
        end
    end
    
    if exist('AOI_ks_adj', 'var') ~= 1 || REGEN == 1
        [AOI_ksn045, AOI_ks_adj, AOI_slopearea] = ...
            DEM_ksn(theta, AOI_FAC, AOI_DEM_gradient8, AOI_rivers_STR, AOI_mg, DEM_ksn045_fname, DEM_ks_adj_fname);
    end
    
    %Generate Stream ORDER files for adding IDs to catchments
    if exist('AOI_STO_N', 'var') ~= 1 || REGEN == 1
        AOI_STO = streamorder(AOI_FD, AOI_rivers_w);
        for k = 1:length(stream_order)
            stream_order2use = stream_order(k);
            [AOI_STO_N{k}] = drainagebasins(AOI_FD, AOI_STO, stream_order2use);
            stream_order_label(k,:) = sprintf('STO%d',stream_order(k));
        end
    end
    
    %extract trunk stream
    if exist('AOI_STR_streams_dbasins_trunk_grid', 'var') ~= 1 || REGEN == 1
        AOI_STR_streams_dbasins = klargestconncomps(AOI_rivers_STR);
        AOI_STR_streams_dbasins_trunk = trunk(AOI_STR_streams_dbasins);
        AOI_STR_streams_dbasins_trunk_grid = STREAMobj2GRIDobj(AOI_STR_streams_dbasins_trunk);
        %idxnan = find(AOI_STR_streams_dbasins_trunk_grid.Z == 0);
        %AOI_STR_streams_dbasins_trunk_grid.Z(idxnan) = NaN;
    end
    
    if MISC_FILES == 1
        AOI_STR_MS = STREAMobj2mapstruct(AOI_rivers_STR,...
            'seglength', segL, 'attributes',...
            {'ksn045' AOI_ksn045 @mean ...
            'ks_adj' AOI_ks_adj @mean ...
            'SSP' AOI_SSP @mean ...
            'rel1' AOI_DEM_rel_1 @mean ...
            'rel2' AOI_DEM_rel_2 @mean ...
            'rel3' AOI_DEM_rel_3 @mean ...
            'uparea_m2' (AOI_FAC.*(AOI_FAC.cellsize^2)) @mean ...
            'gradient' AOI_DEM_gradient8 @mean ...
            'elev' AOI_DEM @mean ...
            'elev_mg' AOI_mg @mean ...
            'prof_curv' AOI_DEM_curv_profc @mean ...
            'plan_curv' AOI_DEM_curv_planc @mean ...
            'mean_curv' AOI_DEM_curv_meanc @mean ...
            'diff_curv' AOI_DEM_diffusionf_curv @mean ...
            'wiene_curv' AOI_DEM_wiener_curv @mean ...
            stream_order_label(1,:) AOI_STO_N{1} @max ...
            stream_order_label(2,:) AOI_STO_N{2} @max ...
            'TRUNK_ID' AOI_STR_streams_dbasins_trunk_grid @max
            });
    else
        AOI_STR_MS = STREAMobj2mapstruct(AOI_rivers_STR,...
            'seglength', segL, 'attributes',...
            {'ksn045' AOI_ksn045 @mean ...
            'ks_adj' AOI_ks_adj @mean ...
            'SSP' AOI_SSP @mean ...
            'uparea_m2' (AOI_FAC.*(AOI_FAC.cellsize^2)) @mean ...
            'gradient' AOI_DEM_gradient8 @mean ...
            'elev' AOI_DEM @mean ...
            'elev_mg' AOI_mg @mean ...
            });
    end
    
    if length(AOI_rivers_STR_area1.x) > 5
        if MISC_FILES == 1
            AOI_STR_MS_area1 = STREAMobj2mapstruct(AOI_rivers_STR_area1,...
                'seglength', segL, 'attributes',...
                {'ksn045' AOI_ksn045 @mean ...
                'ks_adj' AOI_ks_adj @mean ...
                'SSP' AOI_SSP @mean ...
                'rel1' AOI_DEM_rel_1 @mean ...
                'rel2' AOI_DEM_rel_2 @mean ...
                'rel3' AOI_DEM_rel_3 @mean ...
                'uparea_m2' (AOI_FAC.*(AOI_FAC.cellsize^2)) @mean ...
                'gradient' AOI_DEM_gradient8 @mean ...
                'elev' AOI_DEM @mean ...
                'elev_mg' AOI_mg @mean ...
                'prof_curv' AOI_DEM_curv_profc @mean ...
                'plan_curv' AOI_DEM_curv_planc @mean ...
                'mean_curv' AOI_DEM_curv_meanc @mean ...
                'diff_curv' AOI_DEM_diffusionf_curv @mean ...
                'wiene_curv' AOI_DEM_wiener_curv @mean ...
                stream_order_label(1,:) AOI_STO_N{1} @max ...
                stream_order_label(2,:) AOI_STO_N{2} @max ...
                'TRUNK_ID' AOI_STR_streams_dbasins_trunk_grid @max
                });
        else
            AOI_STR_MS_area1 = STREAMobj2mapstruct(AOI_rivers_STR_area1,...
                'seglength', segL, 'attributes',...
                {'ksn045' AOI_ksn045 @mean ...
                'ks_adj' AOI_ks_adj @mean ...
                'SSP' AOI_SSP @mean ...
                'uparea_m2' (AOI_FAC.*(AOI_FAC.cellsize^2)) @mean ...
                'gradient' AOI_DEM_gradient8 @mean ...
                'elev' AOI_DEM @mean ...
                'elev_mg' AOI_mg @mean ...
                });
        end
    else
        AOI_STR_MS_area1 = NaN;
    end
    
    if length(AOI_rivers_STR_area2.x) > 5
        if MISC_FILES == 1
            AOI_STR_MS_area2 = STREAMobj2mapstruct(AOI_rivers_STR_area2,...
                'seglength', segL, 'attributes',...
                {'ksn045' AOI_ksn045 @mean ...
                'ks_adj' AOI_ks_adj @mean ...
                'SSP' AOI_SSP @mean ...
                'rel1' AOI_DEM_rel_1 @mean ...
                'rel2' AOI_DEM_rel_2 @mean ...
                'rel3' AOI_DEM_rel_3 @mean ...
                'uparea_m2' (AOI_FAC.*(AOI_FAC.cellsize^2)) @mean ...
                'gradient' AOI_DEM_gradient8 @mean ...
                'elev' AOI_DEM @mean ...
                'elev_mg' AOI_mg @mean ...
                'prof_curv' AOI_DEM_curv_profc @mean ...
                'plan_curv' AOI_DEM_curv_planc @mean ...
                'mean_curv' AOI_DEM_curv_meanc @mean ...
                'diff_curv' AOI_DEM_diffusionf_curv @mean ...
                'wiene_curv' AOI_DEM_wiener_curv @mean ...
                stream_order_label(1,:) AOI_STO_N{1} @max ...
                stream_order_label(2,:) AOI_STO_N{2} @max ...
                'TRUNK_ID' AOI_STR_streams_dbasins_trunk_grid @max                
                });
            AOI_STR_MS_areaonly = STREAMobj2mapstruct(AOI_rivers_STR_area2,...
                'attributes',...
                {'uparea_m2' (AOI_FAC.*(AOI_FAC.cellsize^2)) @mean ...
                'gradient' AOI_DEM_gradient8 @mean ...
                'elev' AOI_DEM @mean ...
                'elev_mg' AOI_mg @mean ...
                });
        else
            AOI_STR_MS_area2 = STREAMobj2mapstruct(AOI_rivers_STR_area2,...
                'seglength', segL, 'attributes',...
                {'ksn045' AOI_ksn045 @mean ...
                'ks_adj' AOI_ks_adj @mean ...
                'SSP' AOI_SSP @mean ...
                'uparea_m2' (AOI_FAC.*(AOI_FAC.cellsize^2)) @mean ...
                'gradient' AOI_DEM_gradient8 @mean ...
                'elev' AOI_DEM @mean ...
                'elev_mg' AOI_mg @mean ...
                });
            AOI_STR_MS_areaonly = STREAMobj2mapstruct(AOI_rivers_STR_area2,...
                'attributes',...
                {'uparea_m2' (AOI_FAC.*(AOI_FAC.cellsize^2)) @mean ...
                'gradient' AOI_DEM_gradient8 @mean ...
                'elev' AOI_DEM @mean ...
                'elev_mg' AOI_mg @mean ...
                });
            
        end
    else
        AOI_STR_MS_area2 = NaN;
    end
    
    fprintf(1,'\tsaving variables\n');
    % saves all variables so you don't need to repeat the time-intensive steps
    if MISC_FILES == 1
        save(DEM_HYD_MAT_fname, 'AOI_FIL', 'AOI_FD', 'AOI_FAC', 'AOI_resolution', ...
            'minApix', 'AOI_FAC_w', 'AOI_STR_w', 'AOI_mg', 'AOI_STR_MS', ...
            'AOI_DEM_rel*', 'AOI_Q', 'AOI_SSP', 'AOI_DEM_gradient8', 'AOI_dbasins', ...
            'AOI_rivers_slope', 'AOI_rivers_area', 'AOI_ksn045', 'AOI_ks_adj', ...
            'AOI_slopearea', 'AOI_rivers_STR', 'AOI_rivers_STR_area1', ...
            'AOI_rivers_STR_area2', 'AOI_STR_MS_area1', 'AOI_STR_MS_area2', ...
            'AOI_STR_MS_areaonly', 'AOI_rivers_w', 'AOI_DEM_curv_*', 'AOI_DEM_wiene*', ...
            'AOI_DEM_diffusion*', 'AOI_ridgecres*', 'AOI_ridgecrest_MS', ...,
            'ridgecrest_stepsize', '-v7.3');
    else
        save(DEM_HYD_MAT_fname, 'AOI_FIL', 'AOI_FD', 'AOI_FAC', 'AOI_resolution', ...
            'minApix', 'AOI_FAC_w', 'AOI_STR_w', 'AOI_mg', 'AOI_STR_MS', ...
            'AOI_Q', 'AOI_SSP', 'AOI_DEM_gradient8', 'AOI_dbasins', ...
            'AOI_rivers_slope', 'AOI_rivers_area', 'AOI_ksn045', 'AOI_ks_adj', ...
            'AOI_slopearea', 'AOI_rivers_STR', 'AOI_rivers_STR_area1', ...
            'AOI_STR_MS_areaonly', 'AOI_rivers_STR_area2', 'AOI_STR_MS_area1', ...
            'AOI_STR_MS_area2', 'AOI_rivers_w', '-v7.3');
    end
elseif exist(DEM_HYD_MAT_fname, 'file') == 2
    if exist('AOI_STR_MS', 'var') ~= 1 && show_figs ~= 2
        load(DEM_HYD_MAT_fname)
    end
end

if MISC_FILES == 1
    load(DEM_HYD_MAT_fname, 'AOI_DEM_gradient8', 'AOI_rivers_STR', ...
        'AOI_STR_MS*', 'AOI_DEM_curv_*', 'AOI_DEM_rel_2', 'AOI_STR_MS_areaonly', ...
        'AOI_DEM_wiener_curv', 'AOI_DEM_diffusionf_curv', ...
        'AOI_slopearea', 'AOI_ridgecres*', 'ridgecrest_stepsize');
    BSN_writeshapefiles_long(dir_sep, DEM_basename, DEM_basename_nodir, DEM_fname, AOI_STR_MS, ...
        AOI_STR_MS_area1, AOI_STR_MS_area2, AOI_STR_MS_areaonly, AOI_ridgecrest_MS, AOI_ridgecrest_DEM, ...
        str_area1, str_area2, gdalsrsinfo_cmd, ogr2ogr_cmd, remove_cmd, mv_cmd);
    
    shapeout_fn_prj = sprintf('shapefiles%s%s_%s_PT_proj.shp', ...
        dir_sep, DEM_basename, num2str(area_threshold,'%1.0e'));
    if exist(shapeout_fn_prj, 'file') ~= 2 || REGEN == 1
        DEM_FAC_pt_export(AOI_DEM, AOI_FAC, area_threshold, AOI_FD, ...
            AOI_mg, dir_sep, DEM_basename, gdalsrsinfo_cmd, ogr2ogr_cmd, remove_cmd, mv_cmd);
    end
else
    load(DEM_HYD_MAT_fname, 'AOI_DEM_gradient8', 'AOI_rivers_STR', ...
        'AOI_STR_MS*', 'AOI_slopearea');
    BSN_writeshapefiles_short(dir_sep, DEM_basename, DEM_basename_nodir, DEM_fname, AOI_STR_MS, ...
        AOI_STR_MS_area1, AOI_STR_MS_area2, AOI_STR_MS_areaonly, ...
        str_area1, str_area2, gdalsrsinfo_cmd, ogr2ogr_cmd, remove_cmd, mv_cmd);
end

% Create figures showing DEM, hillshade, and river network, save as PDF:
if show_figs == 1 || show_figs == 2
    if show_figs == 2
        load(DEM_MAT_fname)
        load(DEM_HYD_MAT_fname, 'AOI_DEM_gradient8', 'AOI_rivers_STR', ...
            'AOI_STR_MS', 'AOI_DEM_curv_*', 'AOI_DEM_rel_2', ...
            'AOI_DEM_wiener_curv', 'AOI_DEM_diffusionf_curv', ...
            'AOI_slopearea', 'AOI_ridgecres*', 'ridgecrest_stepsize');
    end
    if show_figs == 1
        load(DEM_HYD_MAT_fname, 'AOI_DEM_gradient8', 'AOI_rivers_STR*', ...
            'AOI_STR_MS', 'AOI_DEM_curv_*', 'AOI_DEM_rel_2', ...
            'AOI_DEM_wiener_curv', 'AOI_DEM_diffusionf_curv', ...
            'AOI_slopearea', 'AOI_ridgecres*', 'ridgecrest_stepsize');
    end
    if MISC_FILES == 1
        DEM_mkfigures_long(MATLABV, PaperType_size, quality_flag, DEM_basename_nodir, ...
            AOI_DEM, DEM_basename_no_underscore, AOI_rivers_STR_area2, AOI_DEM_gradient8, ...
            AOI_DEM_curv_profc, AOI_DEM_curv_planc, AOI_DEM_curv_meanc, ...
            AOI_DEM_wiener_curv, AOI_DEM_diffusionf_curv, AOI_STR_MS, ...
            AOI_DEM_rel_2, relief_values_m, AOI_slopearea, AOI_ridgecrest_MS, ...
            AOI_ridgecrest_DEM, AOI_ridgecrest_MS_Dy_all, ...
            AOI_ridgecrest_MS_yi_all, ridgecrest_stepsize)
    else
        DEM_mkfigures_short(MATLABV, PaperType_size, quality_flag, DEM_basename_nodir, ...
            AOI_DEM, DEM_basename_no_underscore, AOI_rivers_STR_area2, AOI_DEM_gradient8, ...
            AOI_STR_MS, AOI_slopearea)
    end
end

clear AOI_FIL AOI_resolution AOI_FAC_w AOI_STR_w ...
    AOI_DEM_rel* AOI_Q AOI_SSP AOI_DEM_curv_* ...
    AOI_rivers_slope AOI_rivers_area AOI_ksn045 AOI_ks_adj ...
    AOI_slopearea AOI_rivers_STR AOI_rivers_STR_area1 ...
    AOI_rivers_STR_area2 AOI_STR_MS_area1AOI_STR_MS_area2 AOI_DEM_wiene* ...
    AOI_DEM_diffusion AOI_ridgecres* AOI_ridgecrest_MS
