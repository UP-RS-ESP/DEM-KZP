%% (3) Calculate  area-slope, ksn, and chi values for all drainage basins contained in this DEM (e.g., draining to ocean/outlets/edge)
%
% If there is only one drainage basin, only one file will be processed
%
% Ths part performs the following steps:
% (1) find dbasins with drainage areas > min_dbasins_stats_to_process
% (often 1e6 m^2)
%
% (2) use AOI_rivers to mask out stream networks (i.e., create
% a new streamnetwork mask: AOI_STR_streams_dbasins_unique
% with trunk streams in: AOI_STR_all_streams_trunk)
%
% (3) calculate slope/area and chiplots for each catchment using trunk
% stream
%
% (4) Write to shapefile/csv

fprintf('KZP topometrics step 3 of 4: area-slope and chi-plot processing for all basins draining to edge\n');
if exist('AOI_DEM', 'var') ~= 1
    load(KZP_parameters.DEM_MAT_fname)
end
if exist('AOI_dbasins', 'var') ~= 1
    load(KZP_parameters.DEM_HYD_MAT_fname, 'AOI_dbasins', 'AOI_FAC', 'AOI_FD', 'minApix', ...
        'AOI_DEM_gradient8', 'AOI_mg', 'AOI_STR_MS')
end
if exist(KZP_parameters.DEM_STR_MAT_fname, 'file') == 2
    load(KZP_parameters.DEM_STR_MAT_fname)
elseif KZP_parameters.show_figs == 1 || KZP_parameters.show_figs == 2
    if exist('AOI_dbasins_unique', 'var') ~= 1 && exist('AOI_dbasins_unique_nr', 'var') ~= 1
        AOI_dbasins_unique = unique(AOI_dbasins.Z);
        if AOI_dbasins_unique(1) == 0
            AOI_dbasins_unique(1) = [];
        end
        AOI_dbasins_unique_stats = regionprops(AOI_dbasins.Z, ...
            'PixelIdxList','Area','Centroid');
        idx_area = find([AOI_dbasins_unique_stats.Area] > ...
            ((KZP_parameters.min_dbasins_stats_to_process)/(AOI_DEM.cellsize.^2)));
        AOI_dbasins_unique_stats = AOI_dbasins_unique_stats(idx_area);
        if length(AOI_dbasins_unique_stats) > 0
            for i = 1:length(AOI_dbasins_unique_stats)
                idx_basinid(i) = ...
                    AOI_dbasins.Z([AOI_dbasins_unique_stats(i).PixelIdxList(1)]);
            end
            % idx2remove = find(idx_basinid == AOI_dbasins_unique);
            AOI_dbasins_unique = idx_basinid;
            AOI_dbasins_unique_nr = numel(AOI_dbasins_unique);
        else
            AOI_dbasins_unique = NaN;
            AOI_dbasins_unique_nr = NaN;
            fprintf('\tCould not find any basins\n');
        end
    else
        if exist('AOI_dbasins_unique', 'var')
            AOI_dbasins_unique_stats = regionprops(AOI_dbasins.Z, ...
                'PixelIdxList','Area','Centroid');
            idx_area = find([AOI_dbasins_unique_stats.Area] > ...
                ((KZP_parameters.min_dbasins_stats_to_process)/(AOI_DEM.cellsize.^2)));
            AOI_dbasins_unique_stats = AOI_dbasins_unique_stats(idx_area);
            if length(AOI_dbasins_unique_stats) > 0
                for i = 1:length(AOI_dbasins_unique_stats)
                    idx_basinid(i) = ...
                        AOI_dbasins.Z([AOI_dbasins_unique_stats(i).PixelIdxList(1)]);
                end
                % idx2remove = find(idx_basinid == AOI_dbasins_unique);
                AOI_dbasins_unique = idx_basinid;
                AOI_dbasins_unique_nr = numel(AOI_dbasins_unique);
            else
                AOI_dbasins_unique = NaN;
                AOI_dbasins_unique_nr = NaN;
                fprintf('\tCould not find any basins\n');
            end
        end
    end
end
if ~isnan(AOI_dbasins_unique_nr)
    AOI_dbasins_stats_min_area = NaN(AOI_dbasins_unique_nr,35);
    % create empty mask to indicate basins
    AOI_dbasins_stats = AOI_dbasins; AOI_dbasins_stats.Z = ...
        zeros(size(AOI_dbasins.Z), 'uint32');
    AOI_dbasins_stats.name='drainage basins w/ stats';
    
    fprintf(['\tCalculating basin statistics for basins >%G m^2 (>%G km^2 or %2.2E m^2).\n',...
        'At basin # of %d: '], KZP_parameters.min_dbasins_stats_to_process, ...
        KZP_parameters.min_dbasins_stats_to_process/1e6, KZP_parameters.min_dbasins_stats_to_process, AOI_dbasins_unique_nr);
    for i = 1:AOI_dbasins_unique_nr
        fprintf('%d, ', i);
        no_plot = 0;
        if mod(i,10) == 0
            fprintf('\n\tAt basin # of %d: ', AOI_dbasins_unique_nr);
        end
        idx_basin = AOI_dbasins_unique(i);
        AOI_dbasins_current = AOI_dbasins;
        idx_ncurrent_basin = find(AOI_dbasins_current.Z ~= idx_basin);
        AOI_dbasins_current.Z(idx_ncurrent_basin) = 0;
        AOI_FAC_dbasins_w = (AOI_FAC.*(AOI_FAC.cellsize.^2) > ...
            KZP_parameters.min_drainage_area_to_process) & (AOI_dbasins_current > 0);  % masks flow accum grid (above threshold)
        AOI_STR_dbasins_unique_subset{i} = STREAMobj(AOI_FD, AOI_FAC_dbasins_w);
        
        % Mask basin in AOI_dbasins_stats
        idx_current_basin = find(AOI_dbasins_current.Z == idx_basin);
        AOI_dbasins_stats.Z(idx_current_basin) = AOI_dbasins_unique(i);
        
        %only calculate log area-log slope plot, if there are more than 20 elements in vector
        if length(AOI_STR_dbasins_unique_subset{i}.x) > 20
            AOI_STR_streams_dbasins_unique{i} = ...
                klargestconncomps(AOI_STR_dbasins_unique_subset{i});
            AOI_STR_all_streams_trunk{i} = trunk(AOI_STR_streams_dbasins_unique{i});
            idx_outlet = find(AOI_STR_all_streams_trunk{i}.distance == ...
                min(AOI_STR_all_streams_trunk{i}.distance));
            BasinO_X = AOI_STR_all_streams_trunk{i}.x(idx_outlet);
            BasinO_Y = AOI_STR_all_streams_trunk{i}.y(idx_outlet);
            
            % calculate slope and area from DEM: entire dataset
            AOI_STR_slope_subset{i} = ...
                gradient(AOI_STR_streams_dbasins_unique{i}, AOI_mg, 'unit','tangent');
            AOI_STR_area_subset{i} = ...
                AOI_FAC.Z(AOI_STR_streams_dbasins_unique{i}.IXgrid).*(AOI_FAC.cellsize).^2;
            % remove gradient less than minimum gradient defined in min_str_gradient
            idx0 = find(AOI_STR_slope_subset{i} <= KZP_parameters.min_str_gradient);
            AOI_STR_slope_subset{i}(idx0) = []; AOI_STR_area_subset{i}(idx0) = [];
            
            % calculate slope and area from DEM: only trunk stream
            AOI_STR_slope_trunk_subset{i} = ...
                gradient(AOI_STR_all_streams_trunk{i}, AOI_mg, 'unit','tangent');
            AOI_STR_area_trunk_subset{i} = ...
                AOI_FAC.Z(AOI_STR_all_streams_trunk{i}.IXgrid).*(AOI_FAC.cellsize).^2;
            % remove gradient less than minimum gradient defined in min_str_gradient
            idx0 = find(AOI_STR_slope_trunk_subset{i} <= KZP_parameters.min_str_gradient);
            AOI_STR_slope_trunk_subset{i}(idx0) = [];
            AOI_STR_area_trunk_subset{i}(idx0) = [];
            
            %fprintf(1,'\tcalculating slope-area relation for all drainages.\n');
            AOI_STR_S_slopearea_dbasins{i} = slopearea(AOI_STR_streams_dbasins_unique{i},...
                AOI_mg, AOI_FAC, 'theta', abs(KZP_parameters.theta), ...
                'areabins', 100, 'areabinlocs', 'median', 'gradaggfun', 'median', ...
                'fitmethod', 'ls', 'streamgradient', 'robust', 'plot', false);
            
            %fprintf(1,'\tcalculating slope-area relation for trunk stream.\n');
            AOI_STR_S_slopearea_dbasins_trunk{i} = slopearea(AOI_STR_all_streams_trunk{i},...
                AOI_mg, AOI_FAC, 'theta', abs(KZP_parameters.theta), ...
                'areabins', 100, 'areabinlocs', 'median', 'gradaggfun', 'median', ...
                'fitmethod', 'ls', 'streamgradient', 'robust', 'plot', false);
            
            %fprintf(1,'\tcalculating slope-area relation for trunk stream with free theta.\n');
            AOI_STR_S_slopearea_dbasins_trunk_adj{i} = slopearea(AOI_STR_all_streams_trunk{i},...
                AOI_mg, AOI_FAC, ...
                'areabins', 100, 'areabinlocs', 'median', 'gradaggfun', 'median', ...
                'fitmethod', 'ls', 'streamgradient', 'robust', 'plot', false);
            
            % calculate ks with adjusted theta for this region
            AOI_ks_adj_crop = AOI_DEM_gradient8 ./ ...
                (AOI_FAC.*(AOI_FAC.cellsize^2)).^...
                AOI_STR_S_slopearea_dbasins_trunk{i}.theta;
            AOI_ks_adj_crop.name = 'ks_adjusted';
            AOI_ks_adj_crop.zunit = 'm^0.9';
            
            AOI_ks045_crop = AOI_DEM_gradient8 ./ ...
                (AOI_FAC.*(AOI_FAC.cellsize^2)).^...
                KZP_parameters.theta;
            AOI_ks045_crop.name = 'ks045';
            AOI_ks045_crop.zunit = 'm^0.9';
            
            
            % Create shapefile with relevant data for area > min_drainage_area_to_process
            % this will take a while
            % fprintf(1,'\tcreating mapstruct.\n');
            AOI_STR_MS_crop{i} = STREAMobj2mapstruct(...
                AOI_STR_streams_dbasins_unique{i},...
                'seglength', KZP_parameters.segL, 'attributes', ...
                {'ks045' AOI_ks045_crop @mean ...
                'ks_adj' AOI_ks_adj_crop @mean });
            shapeout_fn = sprintf('%s%strunk%s%s_MS_trunk_ksn_db%d.shp', KZP_parameters.shapefile_dirname, ...
                KZP_parameters.dir_sep, KZP_parameters.dir_sep, KZP_parameters.DEM_basename, i);
            shapeout_fn_out = sprintf('%s%strunk%s%s_MS_trunk_ksn_db%d.out', KZP_parameters.shapefile_dirname, ...
                KZP_parameters.dir_sep, KZP_parameters.dir_sep, KZP_parameters.DEM_basename, i);
            shapeout_all_fn = sprintf('%s%strunk%s%s_MS_trunk_ksn_db%d.*', KZP_parameters.shapefile_dirname, ...
                KZP_parameters.dir_sep, KZP_parameters.dir_sep, KZP_parameters.DEM_basename, i);
            if exist(shapeout_fn, 'file') ~= 2
                %fprintf(1,'\twriting shapefile and adding projection.\n');
                shapewrite(AOI_STR_MS_crop{i},shapeout_fn);
                %Because shapewrite doesn't add projection information, we have to add
                %these manually via ogr2ogr or ArcMAP (or something similar)
                shapeout_fn_prj = sprintf('%s%strunk%s%s_MS_trunk_ksn_db%d_proj.shp', KZP_parameters.shapefile_dirname, ...
                    KZP_parameters.dir_sep, KZP_parameters.dir_sep, KZP_parameters.DEM_basename, i);
                eval([KZP_parameters.gdalsrsinfo_cmd, ' -o wkt ', KZP_parameters.DEM_fname, '> projection.prj']);
                eval([KZP_parameters.ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj ', ...
                    shapeout_fn_prj, ' ', shapeout_fn, ' 2> ', shapeout_fn_out]);
                eval([KZP_parameters.remove_cmd, ' ', shapeout_all_fn]);
            end
        end
        
        %fprintf(1,'\tcalculate chiplot for all streams.\n');
        %check  if stream network has more than 1 outlet
        %only calculate log area-log slope plot, if there are more than 20 elements in vector
        if length(AOI_STR_dbasins_unique_subset{i}.x) > 20
            nrc = numel(AOI_STR_streams_dbasins_unique{i}.x);
            M   = sparse(double(AOI_STR_streams_dbasins_unique{i}.ix),...
                double(AOI_STR_streams_dbasins_unique{i}.ixc),true,nrc,nrc);
            outlet = sum(M,2) == 0 & sum(M,1)'~=0;
            if nnz(outlet)==1, AOI_STR_S_chiplot{i} = ...
                    chiplot(AOI_STR_streams_dbasins_unique{i}, AOI_mg, ...
                    AOI_FAC, 'a0', KZP_parameters.min_drainage_area_to_process, ...
                    'fitto', 'all', 'plot', false);
            elseif nnz(outlet) > 1, AOI_STR_S_chiplot_outlet(i) = nnz(outlet); end
            clear nrc M outlet
            
            %fprintf(1,'\tcalculate chiplot for trunk stream.\n');
            % Get chiplot info for just the trunkstreams for each stream
            nrc = numel(AOI_STR_all_streams_trunk{i}.x);
            M   = sparse(double(AOI_STR_all_streams_trunk{i}.ix),...
                double(AOI_STR_all_streams_trunk{i}.ixc),true,nrc,nrc);
            outlet = sum(M,2) == 0 & sum(M,1)'~=0;
            if nnz(outlet)==1, AOI_STR_S_trunk_chiplot{i} = ...
                    chiplot(AOI_STR_all_streams_trunk{i}, AOI_mg,AOI_FAC, ...
                    'a0', KZP_parameters.min_drainage_area_to_process, 'trunkstream', ...
                    AOI_STR_all_streams_trunk{i}, 'fitto', 'ts', 'plot', false);
            elseif nnz(outlet) > 1, AOI_STR_S_trunk_chiplot_outlet(i) = nnz(outlet); end
            clear nrc M outlet
            
            %fprintf(1,'\tcalculating adusted ks for trunk stream with theta from chi analysis.\n');
            AOI_STR_S_slopearea_dbasins_trunk_chitheta{i} = slopearea(AOI_STR_all_streams_trunk{i},...
                AOI_mg, AOI_FAC, 'theta', abs(AOI_STR_S_trunk_chiplot{i}.mn), ...
                'areabins', 100, 'areabinlocs', 'median', 'gradaggfun', 'median', ...
                'fitmethod', 'ls', 'streamgradient', 'robust', 'plot', false);
            
            %fprintf(1,'\tgenerate regression statistics based on area and slope.\n');
            %generate fit with Curve Fitting toolbox
            area = AOI_STR_area_trunk_subset{i}; grad = AOI_STR_slope_trunk_subset{i};
            idx0 = find(grad <= KZP_parameters.min_str_gradient*10); grad(idx0) = []; area(idx0) = [];
            if length(area) >= 1000/KZP_parameters.segL
                %enough data points
                xData = area; yData = grad;
                if license('test', 'curve_fitting_toolbox') == 1
                    [xData, yData] = prepareCurveData(xData, yData);
                else
                    idx = find(~isfinite(xData)); xData(idx) = []; yData(idx) = [];
                    xData = xData'; yData = yData';
                end
                % Set up fittype and options:
                % First power law with robust weighting of all data points
                if license('test', 'curve_fitting_toolbox') == 1
                    try
                        ft = fittype( 'power1' ); opts = fitoptions( ft );
                        opts.Display = 'Off'; opts.Lower = [-Inf -Inf];
                        opts.Robust = 'Bisquare'; opts.Upper = [Inf Inf];
                        [fitresult, gof] = fit(xData, yData, ft, opts );
                        ci = confint(fitresult);
                    catch err
                        fitresult.a = NaN;  fitresult.b = NaN;
                        gof.adjrsquare = NaN; gof.rmse = NaN; ci(1:4) = NaN;
                    end
                else
                    [fitr, fitr_stats] = robustfit(log10(xData), log10(yData));
                    fitresult.a = exp(fitr(1));  fitresult.b = fitr(2);
                    sse = fitr_stats.dfe * fitr_stats.robust_s^2;
                    phat = fitr(1) + fitr(2)*log10(xData);
                    ssr = norm(phat-mean(phat))^2;
                    adjrsquare = 1 - sse / (sse + ssr);
                    gof.adjrsquare = adjrsquare; gof.rmse = fitr_stats.robust_s;
                    fitr_ci = nlparci(fitr,fitr_stats.resid,'cov',fitr_stats.covb);
                    fitr_ci = fitr_ci';
                    fitr_ci(1) = fitresult.a - fitr_ci(1);
                    fitr_ci(2) = fitresult.a + fitr_ci(2);
                    fitr_ci(3) = fitresult.b - fitr_ci(3);
                    fitr_ci(4) = fitresult.b + fitr_ci(4);
                    ci(1:4) = fitr_ci;
                end
                % Get p and t statistics using regstats:
                % stats_rweight = regstats(log10(xData), log10(yData), 'linear');
                % p_rweight_a = stats_rweight.tstat.pval(1);
                % p_rweight_b = stats_rweight.tstat.pval(2);
                % Better to use robustfit if looking for p and t-value statistics:
                [fitr, fitr_stats] = robustfit(log10(xData), log10(yData));
                p_rweight_a = fitr_stats.p(1); p_rweight_b = fitr_stats.p(2);
                
                % another option to calculate p value:
                % [R_rweight, p_rweight_corr] = corr(log10(xData), log10(yData));
                
                % Power Law with no weighting
                if license('test', 'curve_fitting_toolbox') == 1
                    try
                        ft = fittype( 'power1' ); opts = fitoptions( ft ); opts.Display = 'Off'; opts.Lower = [-Inf -Inf];
                        opts.Robust = 'Off'; opts.Upper = [Inf Inf];
                        [fitresult2, gof2] = fit(xData, yData, ft, opts ); ci2 = confint(fitresult2);
                    catch err
                        fitresult2.a = NaN; fitresult2.b = NaN; gof2.adjrsquare = NaN; gof2.rmse = NaN; ci2(1:4) = NaN;
                    end
                else
                    [fitr, fitr_stats] = robustfit(log10(xData), log10(yData), 'ols');
                    fitresult2.a = exp(fitr(1)); fitresult2.b = fitr(2);
                    sse = fitr_stats.dfe * fitr_stats.robust_s^2;
                    phat = fitr(1) + fitr(2)*log10(xData);
                    ssr = norm(phat-mean(phat))^2;
                    adjrsquare = 1 - sse / (sse + ssr);
                    gof2.adjrsquare = adjrsquare; gof2.rmse = fitr_stats.robust_s;
                    fitr_ci = nlparci(fitr,fitr_stats.resid,'cov',fitr_stats.covb);
                    fitr_ci = fitr_ci';
                    fitr_ci(1) = fitresult2.a - fitr_ci(1);
                    fitr_ci(2) = fitresult2.a + fitr_ci(2);
                    fitr_ci(3) = fitresult2.b - fitr_ci(3);
                    fitr_ci(4) = fitresult2.b + fitr_ci(4);
                    ci2(1:4) = fitr_ci;
                end
                
                % Power Law of binned data
                min_logspace = floor(log10(KZP_parameters.min_max_DA_fits(1)));
                max_logspace = ceil(log10(KZP_parameters.min_max_DA_fits(2)));
                % create min_max_DA_fits(3) log-spaced points between 10^min_logspace and 10^max_logspace:
                edges_log = logspace(min_logspace, max_logspace, KZP_parameters.min_max_DA_fits(3));
                for j = 1:length(edges_log)-1
                    idx = find((xData) >= edges_log(j) & ...
                        (xData) < edges_log(j+1));
                    xData_median(j) = nanmedian(xData(idx));
                    yData_median(j) = nanmedian(yData(idx));
                    yData_std(j) = nanstd(yData(idx));
                    clear idx
                end
                if sum(isnan(xData_median)) < 20
                    % Not more than 20 values missing
                    if license('test', 'curve_fitting_toolbox') == 1
                        [xData_logspace, yData_logspace, wData_logspace] = prepareCurveData(xData_median, yData_median, yData_std);
                    else
                        idx = find(~isfinite(xData_median)); xData_median(idx) = []; yData_median(idx) = []; yData_std(idx) = [];
                        xData_logspace = xData_median'; yData_logspace = yData_median'; wData_logspace = yData_std';
                    end
                    
                    if license('test', 'curve_fitting_toolbox') == 1
                        try
                            ft = fittype( 'power1' ); opts = fitoptions( ft ); opts.Display = 'Off'; opts.Lower = [-Inf -Inf]; opts.Weights = wData_logspace;
                            opts.Robust = 'Bisquare'; opts.Upper = [Inf Inf];
                            [fitresult_logspace, gof_logspace{i}] = fit(xData_logspace, yData_logspace, ft, opts ); ci_logspace = confint(fitresult_logspace);
                        catch err
                            fitresult_logspace.a = NaN; fitresult_logspace.b = NaN;
                            gof_logspace{i}.adjrsquare = NaN; gof_logspace{i}.rmse = NaN; ci_logspace(1:4) = NaN;
                        end
                    else
                        [fitr, fitr_stats] = robustfit(log10(xData_logspace), log10(yData_logspace));
                        fitresult_logspace.a = exp(fitr(1));  fitresult_logspace.b = fitr(2);
                        sse = fitr_stats.dfe * fitr_stats.robust_s^2;
                        phat = fitr(1) + fitr(2)*log10(xData_logspace);
                        ssr = norm(phat-mean(phat))^2;
                        adjrsquare = 1 - sse / (sse + ssr);
                        gof_logspace{i}.adjrsquare = adjrsquare; gof_logspace{i}.rmse = fitr_stats.robust_s;
                        fitr_ci = nlparci(fitr,fitr_stats.resid,'cov',fitr_stats.covb);
                        fitr_ci = fitr_ci';
                        fitr_ci(1) = fitresult_logspace.a - fitr_ci(1);
                        fitr_ci(2) = fitresult_logspace.a + fitr_ci(2);
                        fitr_ci(3) = fitresult_logspace.b - fitr_ci(3);
                        fitr_ci(4) = fitresult_logspace.b + fitr_ci(4);
                        ci_logspace(1:4) = fitr_ci;
                    end
                    stats_logspace = regstats(log10(xData_logspace), log10(yData_logspace), 'linear');
                    p_logspace_a = stats_logspace.tstat.pval(1); p_logspace_b = stats_logspace.tstat.pval(2);
                else
                    fitresult_logspace.a = NaN; fitresult_logspace.b = NaN;
                    gof_logspace{i}.adjrsquare = NaN; gof_logspace{i}.rmse = NaN; ci_logspace(1:4) = NaN;
                    p_logspace_a = NaN; p_logspace_b = NaN;
                end
                
                
                
                %pixelvalue, Centroid_X, Centroid_Y, BasinO_X, BasinO_Y, ks_slopearea, theta_slopearea, ks_trunk_slopearea, theta_trunk_slopearea, theta (concavity), y interception (ksn), r2, rmse, confint: a-, confint: a+, confint:b-, confint: b+, pval_a, pval_b, theta_robust_off, r2_robust_off, rmse_robust_off
                utm_x_centroid = AOI_x(round(AOI_dbasins_unique_stats(i).Centroid(2)), round(AOI_dbasins_unique_stats(i).Centroid(1)));
                utm_y_centroid = AOI_y(round(AOI_dbasins_unique_stats(i).Centroid(2)), round(AOI_dbasins_unique_stats(i).Centroid(1)));
                AOI_dbasins_stats_min_area(i,:) = [double(AOI_dbasins_unique(i)), double(utm_x_centroid), double(utm_y_centroid), BasinO_X, BasinO_Y, double(AOI_STR_S_slopearea_dbasins{i}.ks), ...
                    double(AOI_STR_S_slopearea_dbasins{i}.theta), double(AOI_STR_S_slopearea_dbasins_trunk{i}.ks), double(AOI_STR_S_slopearea_dbasins_trunk{i}.theta),...
                    fitresult.b, fitresult.a, gof.adjrsquare, gof.rmse, ci(1), ci(2), ci(3), ci(4), p_rweight_a, p_rweight_b, ...
                    fitresult2.b, fitresult2.a, gof2.adjrsquare, gof2.rmse, ...
                    fitresult_logspace.b, fitresult_logspace.a, gof_logspace{i}.adjrsquare, gof_logspace{i}.rmse, p_logspace_a, p_logspace_b, ...
                    AOI_STR_S_slopearea_dbasins_trunk_adj{i}.ks, AOI_STR_S_trunk_chiplot{i}.ks, AOI_STR_S_trunk_chiplot{i}.mn, AOI_STR_S_trunk_chiplot{i}.R2, ...
                    max(AOI_STR_S_slopearea_dbasins{i}.a), max(AOI_STR_S_slopearea_dbasins{i}.a)/1e6];
            else
                gof_logspace{i}.adjrsquare = NaN; gof_logspace{i}.rmse = NaN;
                utm_x_centroid = AOI_x(round(AOI_dbasins_unique_stats(i).Centroid(2)), round(AOI_dbasins_unique_stats(i).Centroid(1)));
                utm_y_centroid = AOI_y(round(AOI_dbasins_unique_stats(i).Centroid(2)), round(AOI_dbasins_unique_stats(i).Centroid(1)));
                AOI_dbasins_stats_min_area(i,:) = [AOI_dbasins_unique(i), utm_x_centroid, utm_y_centroid, BasinO_X, BasinO_Y, double(AOI_STR_S_slopearea_dbasins{i}.ks), ...
                    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
            end
        else
            fprintf(1,'\tBasin does not contain enough data points, skipping to next\n');
            gof_logspace{i}.adjrsquare = NaN; gof_logspace{i}.rmse = NaN;
            utm_x_centroid = AOI_x(round(AOI_dbasins_unique_stats(i).Centroid(2)), round(AOI_dbasins_unique_stats(i).Centroid(1)));
            utm_y_centroid = AOI_y(round(AOI_dbasins_unique_stats(i).Centroid(2)), round(AOI_dbasins_unique_stats(i).Centroid(1)));
            AOI_dbasins_stats_min_area(i,:) = [AOI_dbasins_unique(i), utm_x_centroid, utm_y_centroid, BasinO_X, BasinO_Y, NaN, ...
                NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
        end
        if length(AOI_STR_dbasins_unique_subset{i}.x) > 20
            
            if KZP_parameters.show_figs == 1 || KZP_parameters.show_figs == 2
                figure; set(gcf,'Visible', 'off');
                slope_area_plot_fname = sprintf('%s%s%s_basin_%d_sa.pdf', ...
                    KZP_parameters.plots_dirname, KZP_parameters.dir_sep, KZP_parameters.DEM_basename_nodir, i);
                if exist(slope_area_plot_fname, 'file') ~= 2
                    subplot(2,1,1,'align')
                    set(gcf,'units','normalized','position',[0 0 1 1]);
                    set(gcf, 'PaperOrientation', 'landscape');
                    set(gcf, 'PaperType', KZP_parameters.PaperType_size);
                    try
                        AOI_STR_S_slopearea_dbasins_plot{i} = slopearea(...
                            AOI_STR_all_streams_trunk{i},...
                            AOI_mg, AOI_FAC, 'theta', abs(KZP_parameters.theta), ...
                            'areabins', 100, 'areabinlocs', 'median', ...
                            'gradaggfun', 'median', ...
                            'fitmethod', 'ls', 'streamgradient', 'robust', 'plot', true);
                    catch err
                        fprintf(1, '\tCan not create plot for id: %d\n', i);
                        no_plot = 1;
                    end
                    if no_plot == 0
                        hold on
                        set(AOI_STR_S_slopearea_dbasins_plot{i}.hLine, 'LineWidth', 2);
                        set(AOI_STR_S_slopearea_dbasins_plot{i}.hLine, ...
                            'Color', [0 0.4470 0.7410]);
                        set(AOI_STR_S_slopearea_dbasins_plot{i}.hPoints, 'LineWidth', 2);
                        set(AOI_STR_S_slopearea_dbasins_plot{i}.hPoints, 'MarkerSize', 10);
                        loglog(AOI_STR_area_trunk_subset{i}, ...
                            AOI_STR_slope_trunk_subset{i}, '+', 'Markersize', 2, ...
                            'color', [0.5 0.5 0.5])
                        set(gca,'Children',flipud(get(gca,'Children')));
                        xlabel('Drainage Area (m^2)', 'Fontsize', 12);
                        ylabel('Slope (m/m)', 'Fontsize', 12);
                        title_string = sprintf('DEM: %s, Basin #: %d, Trunk Stream (theta=0.45)', ...
                            KZP_parameters.DEM_basename_no_underscore, i);
                        title(title_string, 'Fontsize', 16), grid;
                        grid on
                        %adding text: ksn, ks, theta, theta_ref, r2 and regression boundaries from SA plots
                        s = sprintf('From slope-area analyses (binned) (DA=%3.2f km^2): k_{sn} trunk = %3.1f with theta = %0.2f, k_{s} trunk = %3.1f with theta = %0.2f\nk_{s} trunk = %3.1f with theta from chi plot = %0.2f, k_{s} trunk = %3.1f +/- %3.1f with theta from robust regression = %0.2f +/- %0.2f\nr^2 = %0.2f, rmse=%2.2e, p(k_{s}) = %1.2e, p(theta) = %1.2e', ...
                            max(AOI_STR_S_slopearea_dbasins{i}.a)/1e6, AOI_STR_S_slopearea_dbasins_trunk{i}.ks, AOI_STR_S_slopearea_dbasins_trunk{i}.theta, ...
                            AOI_STR_S_slopearea_dbasins_trunk_adj{i}.ks, AOI_STR_S_slopearea_dbasins_trunk_adj{i}.theta, ...
                            AOI_STR_S_slopearea_dbasins_trunk_chitheta{i}.ks, AOI_STR_S_slopearea_dbasins_trunk_chitheta{i}.theta, ...
                            fitresult_logspace.a, abs(fitresult_logspace.a-ci_logspace(1)), fitresult_logspace.b, abs(fitresult_logspace.b-ci_logspace(3)), ...
                            gof_logspace{i}.adjrsquare, gof_logspace{i}.rmse, p_logspace_a, p_logspace_b);
                        foo=axis;
                        text(0.25, 0.85, s, 'Fontsize', 14, 'Units', 'normalized')
                        hold off
                        
                        subplot(2,1,2,'align')
                        hold on
                        plot(AOI_STR_S_chiplot{i}.chi, AOI_STR_S_chiplot{i}.elev, 'color',[.5 .5 .5], 'Linewidth', 1)
                        plot(AOI_STR_S_trunk_chiplot{i}.chi, AOI_STR_S_trunk_chiplot{i}.elev, 'color',[0 0 0], 'Linewidth', 2)
                        set(gca,'Children',flipud(get(gca,'Children')));
                        xlabel('\chi (m)', 'Fontsize', 12);
                        ylabel('Elevation (m)', 'Fontsize', 12);
                        title_string = sprintf('DEM: %s, Basin #: %d, Trunk Stream (m/n = %0.2f)', ...
                            KZP_parameters.DEM_basename_no_underscore, i, AOI_STR_S_trunk_chiplot{i}.mn);
                        title(title_string, 'Fontsize', 16), grid;
                        %adding text: ksn, ks, theta, theta_ref, r2 and regression boundaries from SA plots
                        s = sprintf('From Chi analysis (DA=%3.2f km^2):\nk_{s} trunk = %3.1f with m/n= %0.2f\nbeta trunk = %3.1f, betaSE = %1.2e, r^2 = %0.2f\nk_{s} (all streams)= %3.1f with m/n= %0.2f, r^2 = %0.2f', ...
                            max(AOI_STR_S_slopearea_dbasins{i}.a)/1e6, AOI_STR_S_trunk_chiplot{i}.ks, AOI_STR_S_trunk_chiplot{i}.mn, ...
                            AOI_STR_S_trunk_chiplot{i}.beta, AOI_STR_S_trunk_chiplot{i}.betase, AOI_STR_S_trunk_chiplot{i}.R2, ...
                            AOI_STR_S_chiplot{i}.ks, AOI_STR_S_chiplot{i}.mn, AOI_STR_S_chiplot{i}.R2);
                        foo=axis;
                        text(0.65, 0.2, s, 'Fontsize', 16, 'Units', 'normalized');
                        grid on
                        hold off
                        if exist('export_fig') == 2
                            export_fig(slope_area_plot_fname,KZP_parameters.quality_flag,'-pdf');
                        else
                            saveas(gcf,slope_area_plot_fname, 'pdf');
                        end
                    end
                end
                
                slope_area_map_fname = sprintf('%s%s%s_basin_%d_ksn.pdf', ...
                    KZP_parameters.plots_dirname, KZP_parameters.dir_sep, KZP_parameters.DEM_basename_nodir, i);
                if exist(slope_area_map_fname, 'file') ~= 2
                    figure; set(gcf,'Visible', 'off');
                    symbolspec_ksn045 = makesymbolspec('line',...
                        {'ksn045' [prctile([AOI_STR_MS.ksn045], 5) ...
                        prctile([AOI_STR_MS.ksn045], 95)] 'color' jet(6)});
                    symbolspec_ks_adj = makesymbolspec('line',...
                        {'ks_adj' [prctile([AOI_STR_MS_crop{i}.ks_adj], 5) ...
                        prctile([AOI_STR_MS_crop{i}.ks_adj], 95)] 'color' jet(6)});
                    clf
                    set(gcf,'units','normalized','position',[0 0 1 1]);
                    set(gcf, 'PaperOrientation', 'landscape');
                    set(gcf, 'PaperType', KZP_parameters.PaperType_size);
                    subplot(2,1,1,'align')
                    imageschs(AOI_DEM, AOI_DEM, 'caxis', ...
                        [floor(min(AOI_DEM(:))) ceil(max(AOI_DEM(:)))], ...
                        'colormap',gray,'colorbar',false)
                    hold on
                    mapshow(AOI_STR_MS,'SymbolSpec',symbolspec_ksn045);
                    ylabel('UTM-Northing (m)', 'Fontsize', 12);
                    xlabel('UTM-Easting (m)', 'Fontsize', 12);
                    title_string = sprintf('%s: K_{sn} with theta -0.45 ', ...
                        KZP_parameters.DEM_basename_no_underscore);
                    title(title_string, 'Fontsize', 14), grid;
                    colorbar
                    caxis([prctile([AOI_STR_MS.ksn045], 5) ...
                        prctile([AOI_STR_MS.ksn045], 95)])
                    xrange = max(AOI_STR_dbasins_unique_subset{i}.x) - min(AOI_STR_dbasins_unique_subset{i}.x);
                    yrange = max(AOI_STR_dbasins_unique_subset{i}.y) - min(AOI_STR_dbasins_unique_subset{i}.y);
                    axis([min(AOI_STR_dbasins_unique_subset{i}.x)-(xrange/15) ...
                        max(AOI_STR_dbasins_unique_subset{i}.x)+(xrange/15) ...
                        min(AOI_STR_dbasins_unique_subset{i}.y)-(yrange/15) ...
                        max(AOI_STR_dbasins_unique_subset{i}.y)+(yrange/15)])
                    hold off
                    
                    subplot(2,1,2,'align')
                    imageschs(AOI_DEM, AOI_DEM, 'caxis', ...
                        [floor(min(AOI_DEM(:))) ceil(max(AOI_DEM(:)))], ...
                        'colormap',gray,'colorbar',false)
                    hold on
                    mapshow(AOI_STR_MS_crop{i},'SymbolSpec',symbolspec_ks_adj);
                    ylabel('UTM-Northing (m)', 'Fontsize', 12);
                    xlabel('UTM-Easting (m)', 'Fontsize', 12);
                    title_string = sprintf('%s: K_{s} adjusted with theta = %0.3f ', ...
                        KZP_parameters.DEM_basename_no_underscore, AOI_STR_S_slopearea_dbasins{i}.theta);
                    title(title_string, 'Fontsize', 14), grid;
                    colorbar
                    caxis([prctile([AOI_STR_MS_crop{i}.ks_adj], 5) ...
                        prctile([AOI_STR_MS_crop{i}.ks_adj], 95)])
                    xrange = max(AOI_STR_dbasins_unique_subset{i}.x) - min(AOI_STR_dbasins_unique_subset{i}.x);
                    yrange = max(AOI_STR_dbasins_unique_subset{i}.y) - min(AOI_STR_dbasins_unique_subset{i}.y);
                    axis([min(AOI_STR_dbasins_unique_subset{i}.x)-(xrange/15) ...
                        max(AOI_STR_dbasins_unique_subset{i}.x)+(xrange/15) ...
                        min(AOI_STR_dbasins_unique_subset{i}.y)-(yrange/15) ...
                        max(AOI_STR_dbasins_unique_subset{i}.y)+(yrange/15)])
                    hold off
                    if exist('export_fig') == 2
                        export_fig(slope_area_map_fname,KZP_parameters.quality_flag,'-pdf');
                    else
                        saveas(gcf,slope_area_map_fname, 'pdf');
                    end
                end
            end
        end
        
        clear AOI_DEM_gradient8_crop AOI_FAC_crop AOI_ks_adj_crop
        close all
        
        clear area grad xData yData xData_logspace yData_logspace wData_logspace
        % Save basin statistics into table
    end
    fprintf('\nSaving Matlab file: %s\n', KZP_parameters.DEM_STR_MAT_fname);
    save (KZP_parameters.DEM_STR_MAT_fname,'AOI_STR_S_trunk_chiplot','AOI_STR_S_chiplot', 'AOI_STR_S_slopearea_dbasins',...
        'AOI_STR_all_streams_trunk', 'AOI_STR_streams_dbasins_unique', 'AOI_STR_slope_subset', ...
        'AOI_STR_area_subset', 'AOI_STR_dbasins_unique_subset', 'AOI_STR_MS_crop', 'AOI_STR_slope_trunk_subset', ...
        'AOI_STR_area_trunk_subset', 'AOI_STR_S_slopearea_dbasins_trunk', 'AOI_dbasins_stats_min_area', ...
        'AOI_dbasins_unique_nr', 'AOI_dbasins_unique', 'AOI_dbasins_stats', 'fitresult_logspace', 'ci_logspace', ...
        'AOI_STR_S_slopearea_dbasins_trunk', 'AOI_ks045_crop', 'AOI_STR_S_slopearea_dbasins_trunk_adj', 'p_logspace_a', 'p_logspace_b', ...
        'gof_logspace', 'AOI_STR_S_slopearea_dbasins_trunk_chitheta', 'AOI_dbasins_unique_stats', '-v7.3');
end

if exist(KZP_parameters.DEM_STR_MAT_fname, 'file') == 2
    if exist('AOI_dbasins_unique_nr', 'var') ~= 1
        load(KZP_parameters.DEM_STR_MAT_fname)
    end
    for i = 1:AOI_dbasins_unique_nr
        no_plot = 0;
        %only calculate log area-log slope plot, if there are more than 20 elements in vector
        if length(AOI_STR_dbasins_unique_subset{i}.x) > 20
            % Create shapefile with relevant data for area > min_drainage_area_to_process
            % this will take a while
            %fprintf(1,'\tcreating mapstruct.\n');
            shapeout_fn = sprintf('%s%strunk%s%s_MS_trunk_ksn_db%d.shp', KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep, KZP_parameters.dir_sep, KZP_parameters.DEM_basename, i);
            shapeout_fn_out = sprintf('%s%strunk%s%s_MS_trunk_ksn_db%d.out', KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep, KZP_parameters.dir_sep, KZP_parameters.DEM_basename, i);
            shapeout_all_fn = sprintf('%s%strunk%s%s_MS_trunk_ksn_db%d.*', KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep, KZP_parameters.dir_sep, KZP_parameters.DEM_basename, i);
            shapeout_fn_prj = sprintf('%s%strunk%s%s_MS_trunk_ksn_db%d_proj.shp', KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep, KZP_parameters.dir_sep, KZP_parameters.DEM_basename, i);
            if exist(shapeout_fn_prj, 'file') ~= 2
                fprintf('\tAt basin %d of %d, ', i, AOI_dbasins_unique_nr);
                fprintf(1,'\twriting shapefile: %s\n', shapeout_fn);
                shapewrite(AOI_STR_MS_crop{i},shapeout_fn);
                %Because shapewrite doesn't add projection information, we have to add
                %these manually via ogr2ogr or ArcMAP (or something similar)
                eval([KZP_parameters.gdalsrsinfo_cmd, ' -o wkt ', DEM_fname, '> projection.prj']);
                eval([KZP_parameters.ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj ', ...
                    shapeout_fn_prj, ' ', shapeout_fn, ' 2>', shapeout_fn_out]);
                eval([KZP_parameters.remove_cmd, ' ', shapeout_all_fn]);
            end
        end
    end
end
% Now could add code to combine shapefiles into one larger shapefile.

if exist('AOI_dbasins_stats_min_area', 'var') ~= 0
    % write grid file containing dbasin mask and link it to
    % statistics file
    headers = {'1ID', '2Centr_X', '3Centr_Y', '4BasinO_X', '5BasinO_Y', '6ksn045', '7theta_sa', ...
        '8ks_trsa', '9the_trsa', '10ks_rfit', '11the_rfit', '12r2_rfit', '13rmse_rf', '14ks_f_ci1', '15ks_f_ci2', '16the_fci1', '17the_fci2', ...
        '18prfit_a', '19prfit_b', '20the_norf', '21ksn_norf', '22r2_nrf', '23rmse_nrf', '24the_lsrf', '25ks_lsrf', '26r2_lsrf', '27rmse_lrf', ...
        '28plsf_a', '29plsf_b', '30ks_adj', '31ksn_chi', '32the_chi', '33the_cR2', '34DA_m2', '35DA_km2'};
    precision = '%8.5f';
    idxnan = find(isnan(AOI_dbasins_stats_min_area));
    AOI_dbasins_stats_min_areab = AOI_dbasins_stats_min_area; AOI_dbasins_stats_min_areab(idxnan) = -9999;
    if exist(KZP_parameters.AOI_dbasins_stats_CNTR_csv_fname, 'file') ~= 2
        fprintf('\tWriting basin csv files\n');
        csvwrite_with_headers(KZP_parameters.AOI_dbasins_stats_CNTR_csv_fname,AOI_dbasins_stats_min_areab,headers,precision);
        csvwrite_with_headers(KZP_parameters.AOI_dbasins_stats_OUT_csv_fname,AOI_dbasins_stats_min_areab,headers,precision);
    end
    if exist(KZP_parameters.AOI_dbasins_stats_fname, 'file') ~= 2
        GRIDobj2geotiff(AOI_dbasins_stats, KZP_parameters.AOI_dbasins_stats_fname);
    end
    % The following command is very time consuming:
    % eval([polygonize_cmd, ' -8 -mask ', AOI_dbasins_stats_fname, ' ', AOI_dbasins_stats_fname, ' -b 1 -f "ESRI Shapefile" ', AOI_dbasins_stats_vector_fname, ' DBasin DBasin-ID']);
    
    % Creating point shapefile that show Basin properties for each basin as
    % point from Centroid_X, Centroid_Y coordinate
    db_stats_shapeout_fn = strcat(KZP_parameters.DEM_basename_nodir, '_db_stats_CNTR.shp');
    db_stats_shapeout_fn_out = strcat(KZP_parameters.DEM_basename_nodir, '_db_stats_CNTR.out');
    db_stats_shapeout_all_fn = strcat(KZP_parameters.DEM_basename_nodir, '_db_stats_CNTR.*');
    db_stats_crt_fn = strcat(KZP_parameters.DEM_basename_nodir, '_db_stats_CNTR.crt');
    if exist(strcat(sprintf('%s%s', KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep),db_stats_shapeout_fn), 'file') ~= 2
        db_stats_FID = fopen(db_stats_crt_fn, 'w+');
        string2write = sprintf(['<OGRVRTDataSource>\n  <OGRVRTLayer name=\"%s\">\n', ...
            '    <SrcDataSource relativeToVRT=\"1\">%s</SrcDataSource>\n', ...
            '    <GeometryType>wkbPoint</GeometryType>\n', ...
            '    <LayerSRS>WGS84</LayerSRS>\n',...
            '    <GeometryField encoding="PointFromColumns" x="2Centr_X" y="3Centr_Y"/>\n',...
            '    <Field name="1ID" type="Integer" width="8"/>\n', ...
            '    <Field name="2Centr_X" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="3Centr_Y" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="4BasinO_X" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="5BasinO_Y" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="6ksn045" type="Real" width="8"/>\n', ...
            '    <Field name="7theta_sa" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="8ks_trsa" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="9the_trsa" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="10ks_rfit" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="11the_rfit" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="12r2_rfit" type="Real" width="8" precision="7"/>\n',  ...
            '    <Field name="13rmse_rf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="14ks_f_ci1" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="15ks_f_ci2" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="16the_fci1" type="Real" width="8"/>\n', ...
            '    <Field name="17the_fci2" type="Real" width="8" precision="7"/>\n',  ...
            '    <Field name="18prfit_a" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="19prfit_b" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="20the_nrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="21ksn_nrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="22r2_norf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="23rmse_nrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="22the_lsrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="23ks_lsrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="24the_lsrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="25ks_lsrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="26r2_lsrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="27rmse_lrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="28plsf_a" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="29plsf_b" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="30ks_adj" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="31ksn_chi" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="32the_chi" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="33the_cR2" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="34DA_m2" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="35DA_km2" type="Real" width="8" precision="7"/>\n', ...
            '  </OGRVRTLayer>\n', '</OGRVRTDataSource>\n'], ...
            strcat(KZP_parameters.DEM_basename_nodir, '_db_stats_CNTR'), KZP_parameters.AOI_dbasins_stats_CNTR_csv_fname);
        fwrite(db_stats_FID, string2write);
        fclose(db_stats_FID);
        eval([KZP_parameters.gdalsrsinfo_cmd, ' -o wkt ', KZP_parameters.DEM_fname, '> projection.prj']);
        eval([KZP_parameters.ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj -f "ESRI Shapefile" ', ...
            db_stats_shapeout_fn, ' ', db_stats_crt_fn, ' 2> ', db_stats_shapeout_fn_out]);
    end
    eval([KZP_parameters.mv_cmd, ' ', db_stats_shapeout_all_fn, ' ', sprintf('%s%s', KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep)]);
    
    % Write shapefile with outlet coordinates
    db_stats_shapeout_fn  = strcat(KZP_parameters.DEM_basename_nodir, '_db_stats_OUT.shp');
    db_stats_shapeout_fn_out  = strcat(KZP_parameters.DEM_basename_nodir, '_db_stats_OUT.out');
    db_stats_shapeout_all_fn  = strcat(KZP_parameters.DEM_basename_nodir, '_db_stats_OUT.*');
    db_stats_crt_fn = strcat(KZP_parameters.DEM_basename_nodir, '_db_stats_OUT.crt');
    if exist(strcat(sprintf('%s%s', KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep),db_stats_shapeout_fn), 'file') ~= 2
        db_stats_FID = fopen(db_stats_crt_fn, 'w+');
        string2write = sprintf(['<OGRVRTDataSource>\n  <OGRVRTLayer name=\"%s\">\n', ...
            '    <SrcDataSource relativeToVRT=\"1\">%s</SrcDataSource>\n', ...
            '    <GeometryType>wkbPoint</GeometryType>\n', ...
            '    <LayerSRS>WGS84</LayerSRS>\n',...
            '    <GeometryField encoding="PointFromColumns" x="4BasinO_X" y="5BasinO_Y"/>\n',...
            '    <Field name="1ID" type="Integer" width="8"/>\n', ...
            '    <Field name="2Centr_X" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="3Centr_Y" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="4BasinO_X" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="5BasinO_Y" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="6ksn045" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="7theta_sa" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="8ks_trsa" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="9the_trsa" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="10ks_rfit" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="11the_rfit" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="12r2_rfit" type="Real" width="8" precision="7"/>\n',  ...
            '    <Field name="13rmse_rf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="14ks_f_ci1" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="15ks_f_ci2" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="16the_fci1" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="17the_fci2" type="Real" width="8" precision="7"/>\n',  ...
            '    <Field name="18prfit_a" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="19prfit_b" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="20the_nrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="21ksn_nrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="22r2_norf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="23rmse_nrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="22the_lsrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="23ks_lsrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="24the_lsrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="25ks_lsrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="26r2_lsrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="27rmse_lrf" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="28plsf_a" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="29plsf_b" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="30ks_adj" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="31ksn_chi" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="32the_chi" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="33the_cR2" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="34DA_m2" type="Real" width="8" precision="7"/>\n', ...
            '    <Field name="35DA_km2" type="Real" width="8" precision="7"/>\n', ...
            '  </OGRVRTLayer>\n', '</OGRVRTDataSource>\n'], ...
            strcat(KZP_parameters.DEM_basename_nodir, '_db_stats_OUT'), KZP_parameters.AOI_dbasins_stats_OUT_csv_fname);
        fwrite(db_stats_FID, string2write);
        fclose(db_stats_FID);
        eval([KZP_parameters.gdalsrsinfo_cmd, ' -o wkt ', KZP_parameters.DEM_fname, '> projection.prj']);
        eval([KZP_parameters.ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj -f "ESRI Shapefile" ', ...
            db_stats_shapeout_fn, ' ', db_stats_crt_fn, ' 2> ', db_stats_shapeout_fn_out]);
    end
    eval([KZP_parameters.mv_cmd, ' ', db_stats_shapeout_all_fn, ' ', sprintf('%s%s', KZP_parameters.shapefile_dirname, KZP_parameters.dir_sep)]);
end
