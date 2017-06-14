% (4) Process individual StreamOrder (STO) basins
% Calculate stream order for basin delineation
if exist('AOI_DEM', 'var') ~= 1
    load(KZP_parameters.DEM_MAT_fname)
end
if exist('AOI_rivers_w', 'var') ~= 1
    load(KZP_parameters.DEM_HYD_MAT_fname, 'AOI_dbasins', 'AOI_FAC', 'AOI_FD', 'minApix', ...
        'AOI_DEM_gradient8', 'AOI_mg', 'AOI_STR_MS', 'AOI_STO_N', 'AOI_rivers_w')
end
clear AOI_STO_N AOI_STO_N_unique
AOI_STO = streamorder(AOI_FD, AOI_rivers_w);

for k = 1:length(KZP_parameters.stream_order)
    stream_order2use = KZP_parameters.stream_order(k);
    fprintf('KZP topometrics step 4 of 4: area-slope and chi-plot processing for stream order basins (STO = %d)\n', stream_order2use);
    AOI_STO_dbasins_stats_CNTR_csv_fname = sprintf('%s_STO_%02d_db_stats_CNTR.csv', KZP_parameters.DEM_basename_nodir, stream_order2use);
    AOI_STO_dbasins_stats_OUT_csv_fname = sprintf('%s_STO_%02d_db_stats_OUT.csv', KZP_parameters.DEM_basename_nodir, stream_order2use);
    if exist(strcat(sprintf('shapefiles%sSTO%s', KZP_parameters.dir_sep, KZP_parameters.dir_sep), AOI_STO_dbasins_stats_CNTR_csv_fname), 'file') == 2 && ...
            exist(strcat(sprintf('shapefiles%sSTO%s', KZP_parameters.dir_sep, KZP_parameters.dir_sep), AOI_STO_dbasins_stats_OUT_csv_fname), 'file') == 2
        continue
    end
    
    [AOI_STO_N, AOI_STO_N_outlet] = drainagebasins(AOI_FD, AOI_STO, stream_order2use);
    AOI_STO_N_fname = sprintf('%s%s%s_STO_%02d.tif', KZP_parameters.geotiff_dirname, KZP_parameters.dir_sep, KZP_parameters.DEM_basename_nodir, stream_order2use);
    if exist(AOI_STO_N_fname, 'file') ~= 2
        GRIDobj2geotiff(AOI_STO_N,AOI_STO_N_fname);
    end
    
    %     AOI_STO_N_vector_fname = sprintf('%s_STO_%02d_polygon.shp', DEM_basename_nodir, stream_order2use);
    %     if exist(AOI_STO_N_vector_fname, 'file') ~= 2
    %         eval([polygonize_cmd, ' -8 -mask ', AOI_STO_N_fname, ' ', AOI_STO_N_fname, ' -b 1 -f "ESRI Shapefile" ', AOI_STO_N_vector_fname, ' DBasin DBasin-ID']);
    %     end
    
    if exist('AOI_STO_N_unique', 'var') ~= 1
        AOI_STO_N_unique = unique(AOI_STO_N.Z); % unique drainage basins
        if (AOI_STO_N_unique(1)) == 0
            AOI_STO_N_unique(1) = [];
        end
        AOI_STO_N_nr = numel(AOI_STO_N_unique); % nr of drainage basins
        if AOI_STO_N_nr == 0
            fprintf(1, '\tNo drainage basins for this stream order (%d).\nSkipping and continue to next stream order.\n', stream_order2use);
            clear AOI_STO_N_unique
            continue;
        end
        
        AOI_STO_N_stats = regionprops(AOI_STO_N.Z,'PixelIdxList','Area','Centroid');
        counter = size(AOI_STO_N_stats,1);
        if KZP_parameters.RIDGECREST == 1
            % extract ridgecrests
            [AOI_ridgecrest_STO_MS{k}] = ...
                DEM_ridgecrest(AOI_STO_N_stats, AOI_STO_N, AOI_STO_N_outlet);
            [AOI_ridgecrest_STO_DEM{k}] = ...
                DEM_ridgevalues(AOI_ridgecrest_STO_MS{k}, AOI_DEM);
        end
        
        AOI_STO_dbasins_stats = NaN(counter,35);
        AOI_STO_pixelvalue_N = NaN(counter,1);
        AOI_trunk_N_outlet = NaN(counter,1);
        fprintf('\tat STO %d basin x of %d: ', stream_order2use, counter);
        for i = 1:counter
            if i == 1 || mod(i, 10)==0, fprintf('%d, ', i);end
            if mod(i,150)==0, fprintf('\n\tat STO %d basin x of %d:', stream_order2use, counter);end
            
            %clip out area of interest (catchment with given streamorder)
            AOI_STO_msk = AOI_STO_N;
            pixelvalue =  double(AOI_STO_N.Z(AOI_STO_N_stats(i).PixelIdxList(1)));
            AOI_STO_pixelvalue_N(i) = double(pixelvalue);
            idx2remove = find(AOI_STO_N.Z ~= pixelvalue);
            AOI_STO_msk.Z(idx2remove) = 0;
            idx2remove = find(AOI_STO_N.Z ~= pixelvalue);
            AOI_FAC_msk = AOI_FAC;
            AOI_FAC_msk.Z(idx2remove) = 1;
            AOI_STO_msk_N = STREAMobj(AOI_FD, AOI_FAC_msk > (KZP_parameters.min_drainage_area_to_process./(AOI_FAC_msk.cellsize^2)));
            AOI_STO_trunk_N(i) = trunk(AOI_STO_msk_N);
            
            if length(AOI_STO_msk_N.x) > 20
                AOI_STO_sa_N(i) = slopearea(AOI_STO_msk_N, AOI_DEM, AOI_FAC_msk, 'areabinlocs', 'median', 'gradaggfun', 'median', 'streamgradient', 'robust', 'plot', false);
                idx_outlet = find(AOI_STO_trunk_N(i).distance == min(AOI_STO_trunk_N(i).distance));
                BasinO_X = AOI_STO_trunk_N(i).x(idx_outlet);
                BasinO_Y = AOI_STO_trunk_N(i).y(idx_outlet);
                
                % calculate slope and area from DEM: entire dataset
                AOI_STO_STR_slope_subset{i} = gradient(AOI_STO_msk_N, AOI_mg, 'unit','tangent');
                AOI_STO_STR_area_subset{i} = AOI_FAC.Z(AOI_STO_msk_N.IXgrid).*(AOI_FAC.cellsize).^2;
                % remove gradient less than minimum gradient defined in min_str_gradient
                idx0 = find(AOI_STO_STR_slope_subset{i} <= KZP_parameters.min_str_gradient);
                AOI_STO_STR_slope_subset{i}(idx0) = []; AOI_STO_STR_area_subset{i}(idx0) = [];
                
                % calculate slope and area from DEM: only trunk stream
                AOI_STO_STR_slope_trunk_subset{i} = gradient(AOI_STO_trunk_N(i), AOI_mg, 'unit','tangent');
                AOI_STO_STR_area_trunk_subset{i} = AOI_FAC.Z(AOI_STO_trunk_N(i).IXgrid).*(AOI_FAC.cellsize).^2;
                % remove gradient less than minimum gradient defined in min_str_gradient
                idx0 = find(AOI_STO_STR_slope_trunk_subset{i} <= KZP_parameters.min_str_gradient);
                AOI_STO_STR_slope_trunk_subset{i}(idx0) = []; AOI_STO_STR_area_trunk_subset{i}(idx0) = [];
                
                if length(AOI_STO_msk_N.ix) > 5
                    AOI_STO_STR_S_slopearea{i} = slopearea(AOI_STO_msk_N,...
                        AOI_mg, AOI_FAC, ...
                        'areabins', 100, 'areabinlocs', 'median', 'gradaggfun', 'median', ...
                        'fitmethod', 'ls', 'streamgradient', 'robust', 'plot', false);
                    if length(AOI_STO_STR_S_slopearea{i}.a) < 2
                        AOI_STO_STR_S_slopearea{i}.a = NaN; AOI_STO_STR_S_slopearea{i}.g = NaN; AOI_STO_STR_S_slopearea{i}.ks = NaN; AOI_STO_STR_S_slopearea{i}.theta = NaN;
                    end
                else
                    AOI_STO_STR_S_slopearea{i}.a = NaN; AOI_STO_STR_S_slopearea{i}.g = NaN; AOI_STO_STR_S_slopearea{i}.ks = NaN; AOI_STO_STR_S_slopearea{i}.theta = NaN;
                end
                
                AOI_STO_STR_S_slopearea_trunk{i} = slopearea(AOI_STO_trunk_N(i),...
                    AOI_mg, AOI_FAC, ...
                    'areabins', 100, 'areabinlocs', 'median', 'gradaggfun', 'median', ...
                    'fitmethod', 'ls', 'streamgradient', 'robust', 'plot', false);
                if length(AOI_STO_STR_S_slopearea_trunk{i}.a) < 2
                    AOI_STO_STR_S_slopearea_trunk{i}.a = NaN; AOI_STO_STR_S_slopearea_trunk{i}.g = NaN; AOI_STO_STR_S_slopearea_trunk{i}.ks = NaN; AOI_STO_STR_S_slopearea_trunk{i}.theta = NaN;
                end
                
                %check  if stream network has more than 1 outlet
                nrc = numel(AOI_STO_msk_N.x);
                M   = sparse(double(AOI_STO_msk_N.ix),double(AOI_STO_msk_N.ixc),true,nrc,nrc);
                outlet = sum(M,2) == 0 & sum(M,1)'~=0;
                if nnz(outlet)==1, AOI_STO_STR_S_chiplot{i} = chiplot(AOI_STO_msk_N, AOI_mg, AOI_FAC, 'a0', KZP_parameters.min_drainage_area_to_process, 'fitto', 'all', 'plot', false);
                elseif nnz(outlet) > 1, AOI_STO_STR_S_chiplot_outlet(i) = nnz(outlet); end
                clear nrc M outlet
                
                % Get chiplot info for just the trunkstreams for each stream
                nrc = numel(AOI_STO_trunk_N(i).x);
                M   = sparse(double(AOI_STO_trunk_N(i).ix),double(AOI_STO_trunk_N(i).ixc),true,nrc,nrc);
                outlet = sum(M,2) == 0 & sum(M,1)'~=0;
                if nnz(outlet)==1
                    AOI_STO_STR_S_trunk_chiplot{i} = chiplot(AOI_STO_trunk_N(i), AOI_mg, AOI_FAC, 'a0', KZP_parameters.min_drainage_area_to_process, 'trunkstream', AOI_STO_trunk_N(i), 'fitto', 'ts', 'plot', false);
                elseif nnz(outlet) > 1, AOI_STO_STR_S_trunk_chiplot_outlet(i) = nnz(outlet); end
                clear nrc M outlet
                
                %fprintf(1,'\tcalculating adusted ks for trunk stream with theta from chi analysis.\n');
                AOI_STO_STR_S_slopearea_trunk_ks_adj{i} = slopearea(AOI_STO_trunk_N(i),...
                    AOI_mg, AOI_FAC, 'theta', abs(AOI_STO_STR_S_trunk_chiplot{i}.mn), ...
                    'areabins', 100, 'areabinlocs', 'median', 'gradaggfun', 'median', ...
                    'fitmethod', 'ls', 'streamgradient', 'robust', 'plot', false);
                
                %generate fit with Curve Fitting toolbox - use only
                %mainstem for fitting
                area = AOI_STO_STR_area_trunk_subset{i};
                grad = AOI_STO_STR_slope_trunk_subset{i};
                idx0 = find(grad <= KZP_parameters.min_str_gradient);
                grad(idx0) = []; area(idx0) = [];
                if length(area) >= 10
                    %enough data points
                    xData = area; yData = grad;
                    if license('test', 'curve_fitting_toolbox') == 1
                        [xData, yData] = prepareCurveData(xData, yData);
                    else
                        idx = find(~isfinite(xData));
                        xData(idx) = []; yData(idx) = [];
                        xData = xData'; yData = yData';
                    end
                    
                    % Set up fittype and options:
                    % First power law with robust weighting of all data points
                    ft = fittype( 'power1' ); opts = fitoptions( ft );
                    opts.Display = 'Off'; opts.Lower = [-Inf -Inf];
                    opts.Robust = 'Bisquare'; opts.Upper = [Inf Inf];
                    opts.MaxFunEvals=1000;
                    if license('test', 'curve_fitting_toolbox') == 1
                        try
                            [fitresult, gof] = fit(xData, yData, ft, opts ); ci = confint(fitresult);
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
                        clear sse phat ssr adjrsquare
                    end
                    % Get p and t statistics using regstats:
                    %stats_rweight = regstats(log10(xData), log10(yData), 'linear');
                    %p_rweight_a = stats_rweight.tstat.pval(1); p_rweight_b = stats_rweight.tstat.pval(2);
                    % Better to use robustfit if looking for p and t-value statistics:
                    [fitr, fitr_stats] = robustfit(log10(xData), log10(yData));
                    p_rweight_a = fitr_stats.p(1); p_rweight_b = fitr_stats.p(2);
                    
                    % Power Law with no weighting
                    ft = fittype( 'power1' ); opts = fitoptions( ft ); opts.Display = 'Off'; opts.Lower = [-Inf -Inf];
                    opts.Robust = 'Off'; opts.Upper = [Inf Inf]; opts.MaxFunEvals=1000;
                    if license('test', 'curve_fitting_toolbox') == 1
                        try
                            [fitresult2, gof2] = fit(xData, yData, ft, opts ); ci2 = confint(fitresult2);
                        catch err
                            fitresult2.a = NaN; fitresult2.b = NaN; gof2.adjrsquare = NaN; gof2.rmse = NaN; ci2(1:4) = NaN;
                        end
                    else
                        [fitr, fitr_stats] = robustfit(log10(xData), log10(yData));
                        fitresult2.a = exp(fitr(1));  fitresult2.b = fitr(2);
                        sse = fitr_stats.dfe * fitr_stats.robust_s^2;
                        phat = fitr(1) + fitr(2)*log10(xData);
                        ssr = norm(phat-mean(phat))^2;
                        adjrsquare = 1 - sse / (sse + ssr);
                        gof2.adjrsquare = adjrsquare; gof2.rmse = fitr_stats.robust_s;
                        fitr_ci = nlparci(fitr,fitr_stats.resid,'cov',fitr_stats.covb);
                        fitr_ci = fitr_ci';
                        fitr_ci(1) = fitresult.a - fitr_ci(1);
                        fitr_ci(2) = fitresult.a + fitr_ci(2);
                        fitr_ci(3) = fitresult.b - fitr_ci(3);
                        fitr_ci(4) = fitresult.b + fitr_ci(4);
                        ci2(1:4) = fitr_ci;
                        clear sse phat ssr adjrsquare
                    end
                    
                    % Power Law of binned data
                    min_logspace = floor(log10(KZP_parameters.min_max_DA_fits(1)));
                    max_logspace = ceil(log10(KZP_parameters.min_max_DA_fits(2)));
                    % create n numbers of log-spaced points (stored in min_max_DA_fits(3))
                    % between 10^min_logspace and 10^max_logspace:
                    edges_log = logspace(min_logspace, max_logspace, KZP_parameters.min_max_DA_fits(3));
                    for j = 1:length(edges_log)-1
                        idx = find((xData) >= edges_log(j) & ...
                            (xData) < edges_log(j+1));
                        xData_median(j) = nanmedian(xData(idx));
                        yData_median(j) = nanmedian(yData(idx));
                        yData_std(j) = nanstd(yData(idx));
                        clear idx
                    end
                    if license('test', 'curve_fitting_toolbox') == 1
                        [xData_logspace, yData_logspace, wData_logspace] = prepareCurveData(xData_median, yData_median, yData_std);
                    else
                        idx = find(~isfinite(xData_median)); xData_median(idx) = []; yData_median(idx) = []; yData_std(idx) = [];
                        xData_logspace = xData_median'; yData_logspace = yData_median'; wData_logspace = yData_std';
                    end
                    if length(xData_logspace) < 5
                        fprintf(1, '\n\tat STO %d: %d. Number of log-binned values: %d, no log-bin fitting performed\n', stream_order2use, i, length(xData_logspace));
                        fitresult_logspace.a = NaN; fitresult_logspace.b = NaN; gof_logspace2.adjrsquare = NaN; gof_logspace2.rmse = NaN; ci_logspace(1:4) = NaN;
                        p_logspace_a = NaN; p_logspace_b = NaN;
                    else
                        if license('test', 'curve_fitting_toolbox') == 1
                            ft = fittype( 'power1' ); opts = fitoptions( ft ); opts.Display = 'Off'; opts.Lower = [-Inf -Inf]; opts.Weights = wData_logspace;
                            opts.Robust = 'Bisquare'; opts.Upper = [Inf Inf]; opts.MaxFunEvals=1000;
                            try
                                [fitresult_logspace, gof_logspace2] = fit(xData_logspace, yData_logspace, ft, opts ); ci_logspace = confint(fitresult);
                            catch err
                                fitresult_logspace.a = NaN; fitresult_logspace.b = NaN; gof_logspace2.adjrsquare = NaN; gof_logspace2.rmse = NaN; ci_logspace(1:4) = NaN;
                            end
                        else
                            [fitr, fitr_stats] = robustfit(log10(xData), log10(yData));
                            fitresult_logspace.a = exp(fitr(1));  fitresult_logspace.b = fitr(2);
                            sse = fitr_stats.dfe * fitr_stats.robust_s^2;
                            phat = fitr(1) + fitr(2)*log10(xData);
                            ssr = norm(phat-mean(phat))^2;
                            adjrsquare = 1 - sse / (sse + ssr);
                            gof_logspace2.adjrsquare = adjrsquare;
                            gof_logspace2.rmse = fitr_stats.robust_s;
                            fitr_ci = nlparci(fitr,fitr_stats.resid,'cov',fitr_stats.covb);
                            fitr_ci = fitr_ci';
                            fitr_ci(1) = fitresult.a - fitr_ci(1);
                            fitr_ci(2) = fitresult.a + fitr_ci(2);
                            fitr_ci(3) = fitresult.b - fitr_ci(3);
                            fitr_ci(4) = fitresult.b + fitr_ci(4);
                            ci_logspace(1:4) = fitr_ci;
                            clear sse phat ssr adjrsquare
                        end
                        stats_logspace = regstats(log10(xData_logspace), log10(yData_logspace), 'linear');
                        p_logspace_a = stats_logspace.tstat.pval(1); p_logspace_b = stats_logspace.tstat.pval(2);
                    end
                    
                    %pixelvalue, Centroid_X, Centroid_Y, BasinO_X, BasinO_Y, ks_slopearea, theta_slopearea, ks_trunk_slopearea, theta_trunk_slopearea, theta (concavity), y interception (ksn), r2, rmse, confint: a-, confint: a+, confint:b-, confint: b+, theta_robust_off, r2_robust_off, rmse_robust_off
                    utm_x_centroid = AOI_x(round(AOI_STO_N_stats(i).Centroid(2)), round(AOI_STO_N_stats(i).Centroid(1)));
                    utm_y_centroid = AOI_y(round(AOI_STO_N_stats(i).Centroid(2)), round(AOI_STO_N_stats(i).Centroid(1)));
                    % Save basin statistics into table
                    
                    AOI_STO_dbasins_stats(i,:) = [double(AOI_STO_N_unique(i)), double(utm_x_centroid), double(utm_y_centroid), BasinO_X, BasinO_Y, double(AOI_STO_STR_S_slopearea{i}.ks), ...
                        double(AOI_STO_STR_S_slopearea{i}.theta), double(AOI_STO_STR_S_slopearea_trunk{i}.ks), double(AOI_STO_STR_S_slopearea_trunk{i}.theta),...
                        fitresult.b, fitresult.a, gof.adjrsquare, gof.rmse, ci(1), ci(2), ci(3), ci(4), p_rweight_a, p_rweight_b, ...
                        fitresult2.b, fitresult2.a, gof2.adjrsquare, gof2.rmse, ...
                        fitresult_logspace.b, fitresult_logspace.a, gof_logspace2.adjrsquare, gof_logspace2.rmse, p_logspace_a, p_logspace_b, ...
                        AOI_STO_STR_S_slopearea_trunk_ks_adj{i}.ks, AOI_STO_STR_S_trunk_chiplot{i}.ks, AOI_STO_STR_S_trunk_chiplot{i}.mn, AOI_STO_STR_S_trunk_chiplot{i}.R2, ...
                        max(AOI_STO_STR_S_slopearea{i}.a), max(AOI_STO_STR_S_slopearea{i}.a)/1e6];
                elseif length(area) < 10
                    AOI_STO_dbasins_stats(i,:) = [double(AOI_STO_N_unique(i)), utm_x_centroid, utm_y_centroid, BasinO_X, BasinO_Y, double(AOI_STO_STR_S_slopearea{i}.ks), ...
                        double(AOI_STO_STR_S_slopearea{i}.theta), double(AOI_STO_STR_S_slopearea_trunk{i}.ks), double(AOI_STO_STR_S_slopearea_trunk{i}.theta), ...
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, max(AOI_STO_STR_S_slopearea{i}.a), max(AOI_STO_STR_S_slopearea{i}.a)/1e6];
                end
                clear area grad xData yData xData_logspace yData_logspace wData_logspace  p_rweight_a p_rweight_b xData_median yData_median yData_std
            else
                utm_x_centroid = AOI_x(round(AOI_STO_N_stats(i).Centroid(2)), round(AOI_STO_N_stats(i).Centroid(1)));
                utm_y_centroid = AOI_y(round(AOI_STO_N_stats(i).Centroid(2)), round(AOI_STO_N_stats(i).Centroid(1)));
                AOI_STO_dbasins_stats(i,:) = [double(AOI_STO_N_unique(i)), utm_x_centroid, utm_y_centroid, NaN, NaN, NaN, ...
                    NaN, NaN, NaN, ...
                    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
            end
            %calculate characteristic channel length for tributaries/STO
        end
        if KZP_parameters.RIDGECREST == 1
            %save MapStructure of ridgecrest values
            fname2save = sprintf('%s_STO_%02d_MS.mat', KZP_parameters.DEM_basename_nodir, stream_order2use);
            save(fname2save, 'AOI_ridgecrest_STO_MS', 'AOI_ridgecrest_STO_DEM', '-v7.3')
        end
        fprintf('\n')
        %write grid file containing dbasin mask and link it to
        %statistics file
        headers = {'1ID', '2Centr_X', '3Centr_Y', '4BasinO_X', '5BasinO_Y', '6ksn045', '7theta_sa', ...
            '8ks_trsa', '9the_trsa', '10ks_rfit', '11the_rfit', '12r2_rfit', '13rmse_rf', '14ks_f_ci1', '15ks_f_ci2', '16the_fci1', '17the_fci2', ...
            '18prfit_a', '19prfit_b', '20the_norf', '21ksn_norf', '22r2_nrf', '23rmse_nrf', '24the_lsrf', '25ks_lsrf', '26r2_lsrf', '27rmse_lrf', ...
            '28plsf_a', '29plsf_b', '30ks_adj', '31ksn_chi', '32the_chi', '33the_cR2', '34DA_m2', '35DA_km2'};
        precision = '%8.5f';
        idxnan = find(isnan(AOI_STO_dbasins_stats)); AOI_STO_dbasins_statsb = AOI_STO_dbasins_stats; AOI_STO_dbasins_statsb(idxnan) = -9999;
        if exist(strcat(sprintf('shapefiles%sSTO%s', KZP_parameters.dir_sep, KZP_parameters.dir_sep), AOI_STO_dbasins_stats_CNTR_csv_fname), 'file') ~= 2
            csvwrite_with_headers(AOI_STO_dbasins_stats_CNTR_csv_fname, AOI_STO_dbasins_statsb, headers,precision)
        end
        if exist(strcat(sprintf('shapefiles%sSTO%s', KZP_parameters.dir_sep, KZP_parameters.dir_sep), AOI_STO_dbasins_stats_OUT_csv_fname), 'file') ~= 2
            csvwrite_with_headers(AOI_STO_dbasins_stats_OUT_csv_fname, AOI_STO_dbasins_statsb, headers,precision)
        end
        
        % Creating point shapefile that show Basin properties for each basin as
        % point from Centroid_X, Centroid_Y coordinate
        db_stats_shapeout_fn = sprintf('%s_STO_%02d_db_stats_CNTR.shp', KZP_parameters.DEM_basename_nodir, stream_order2use);
        db_stats_shapeout_fn_out = sprintf('%s_STO_%02d_db_stats_CNTR.out', KZP_parameters.DEM_basename_nodir, stream_order2use);
        db_stats_shapeout_all_fn  = sprintf('%s_STO_%02d_db_stats_CNTR.*', KZP_parameters.DEM_basename_nodir, stream_order2use);
        if exist(strcat(sprintf('shapefiles%sSTO%s', KZP_parameters.dir_sep, KZP_parameters.dir_sep), db_stats_shapeout_fn), 'file') ~= 2
            db_stats_crt_fn = sprintf('%s_STO_%02d_db_stats_CNTR.crt', KZP_parameters.DEM_basename_nodir, stream_order2use);
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
                sprintf('%s_STO_%02d_db_stats_CNTR', KZP_parameters.DEM_basename_nodir, stream_order2use), AOI_STO_dbasins_stats_CNTR_csv_fname);
            fwrite(db_stats_FID, string2write);
            fclose(db_stats_FID);
            eval([KZP_parameters.gdalsrsinfo_cmd, ' -o wkt ', KZP_parameters.DEM_fname, '> projection.prj']);
            eval([KZP_parameters.ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj -f "ESRI Shapefile" ', ...
                db_stats_shapeout_fn, ' ', db_stats_crt_fn, ' 2> ', db_stats_shapeout_fn_out]);
            eval([KZP_parameters.mv_cmd, ' ', db_stats_shapeout_all_fn, ' ', sprintf('shapefiles%sSTO%s', KZP_parameters.dir_sep, KZP_parameters.dir_sep)]);
            eval([KZP_parameters.remove_cmd, ' ', sprintf('shapefiles%sSTO%s%s', KZP_parameters.dir_sep, KZP_parameters.dir_sep, db_stats_crt_fn)]);
        end
        
        % Write shapefile with outlet coordinates
        db_stats_shapeout_fn  = sprintf('%s_STO_%02d_db_stats_OUT.shp', KZP_parameters.DEM_basename_nodir, stream_order2use);
        db_stats_shapeout_fn_out  = sprintf('%s_STO_%02d_db_stats_OUT.out', KZP_parameters.DEM_basename_nodir, stream_order2use);
        db_stats_shapeout_all_fn  = sprintf('%s_STO_%02d_db_stats_OUT.*', KZP_parameters.DEM_basename_nodir, stream_order2use);
        if exist(strcat(sprintf('shapefiles%sSTO%s', KZP_parameters.dir_sep, KZP_parameters.dir_sep), db_stats_shapeout_fn), 'file') ~= 2
            db_stats_crt_fn = sprintf('%s_STO_%02d_db_stats_OUT.crt', KZP_parameters.DEM_basename_nodir, stream_order2use);
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
                sprintf('%s_STO_%02d_db_stats_OUT', KZP_parameters.DEM_basename_nodir, stream_order2use), AOI_STO_dbasins_stats_OUT_csv_fname);
            fwrite(db_stats_FID, string2write);
            fclose(db_stats_FID);
            eval([KZP_parameters.gdalsrsinfo_cmd, ' -o wkt ', KZP_parameters.DEM_fname, '> projection.prj']);
            eval([KZP_parameters.ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj -f "ESRI Shapefile" ', ...
                db_stats_shapeout_fn, ' ', db_stats_crt_fn, ' 2> ', db_stats_shapeout_fn_out]);
            eval([KZP_parameters.mv_cmd, ' ', db_stats_shapeout_all_fn, ' ', sprintf('shapefiles%sSTO%s', KZP_parameters.dir_sep, KZP_parameters.dir_sep)]);
            eval([KZP_parameters.remove_cmd, ' ', sprintf('shapefiles%sSTO%s%s', KZP_parameters.dir_sep, KZP_parameters.dir_sep, db_stats_crt_fn)]);
        end
    end
    clear AOI_STO_N_unique AOI_STO_N
    % merge files for this streamorder
    %         string2write = sprintf(['<OGRVRTDataSource>\n  <OGRVRTLayer name=\"%s\">\n', ...
    %                 '    <SrcDataSource>%s</SrcDataSource>\n', ...
    %                 '</OGRVRTLayer>\n', ...
    %                 '  </OGRVRTLayer>\n', '</OGRVRTDataSource>\n'],...
    %                 Layername1, Srcname_shp);
    % <OGRVRTDataSource>
    %     <OGRVRTUnionLayer name="unionLayer">
    %         <OGRVRTLayer name="source1">
    %             <SrcDataSource>source1.shp</SrcDataSource>
    %         </OGRVRTLayer>
    %         <OGRVRTLayer name="source2">
    %             <SrcDataSource>source2.shp</SrcDataSource>
    %         </OGRVRTLayer>
    %     </OGRVRTUnionLayer>
    % </OGRVRTDataSource>
    
end
close all
fprintf('\nfinished with KZP_topometrics.\n');
