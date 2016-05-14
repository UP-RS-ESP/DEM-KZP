function DEM_mkfigures_long(MATLABV, PaperType_size, quality_flag, DEM_basename_nodir, ...
    AOI_DEM, DEM_basename_no_underscore, AOI_rivers_STR_area2, AOI_DEM_gradient8, ...
    AOI_DEM_curv_profc, AOI_DEM_curv_planc, AOI_DEM_curv_meanc, ...
    AOI_DEM_wiener_curv, AOI_DEM_diffusionf_curv, AOI_STR_MS, ...
    AOI_DEM_rel_2, relief_values_m, AOI_slopearea, AOI_ridgecrest_MS, ...
    AOI_ridgecrest_DEM, AOI_ridgecrest_MS_Dy_all, AOI_ridgecrest_MS_yi_all, ridgecrest_stepsize)

figure
dem_plot_fname = sprintf('maps/%s_DEM_Slope.pdf', DEM_basename_nodir);
if exist(dem_plot_fname, 'file') ~= 2
    fprintf(1,'\tgenerating slope figures\n');
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', PaperType_size);
    subplot(2,1,1,'align')
    imageschs(AOI_DEM, AOI_DEM, 'caxis', [floor(min(AOI_DEM(:))) ...
        ceil(max(AOI_DEM(:)))])
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: DEM ', DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14), grid;
    hold on
    plot(AOI_rivers_STR_area2, 'k', 'Linewidth', 2); hold off
    
    subplot(2,1,2,'align')
    imageschs(AOI_DEM, AOI_DEM_gradient8, 'caxis', [0 1])
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: Slope ', DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14), grid;
    export_fig(dem_plot_fname,quality_flag,'-pdf');
end
close all

figure
dem_plot_fname = sprintf('maps/%s_Curv_Relief.pdf', DEM_basename_nodir);
if exist(dem_plot_fname, 'file') ~= 2
    fprintf(1,'\tgenerating curvature figures\n');
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', PaperType_size);
    subplot(2,2,1,'align')
    if MATLABV == 1
        min_curv = round(prctile(AOI_DEM_curv_profc.Z(:),5),2);
        max_curv = round(prctile(AOI_DEM_curv_profc.Z(:),95),2);
    else
        min_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_curv_profc.Z(:),5)));
        max_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_curv_profc.Z(:),95)));
    end
    imageschs(AOI_DEM, AOI_DEM_curv_profc, 'caxis', ...
        [min_curv max_curv])
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: Profile Curvature ', ...
        DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14), grid;
    
    subplot(2,2,2,'align')
    if MATLABV == 1
        min_curv = round(prctile(AOI_DEM_curv_planc.Z(:),5),2);
        max_curv = round(prctile(AOI_DEM_curv_planc.Z(:),95),2);
    else
        min_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_curv_planc.Z(:),5)));
        max_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_curv_planc.Z(:),95)));
    end
    imageschs(AOI_DEM, AOI_DEM_curv_planc, 'caxis', ...
        [min_curv max_curv])
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: Planform Curvature ', ...
        DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14), grid;
    
    subplot(2,2,3,'align')
    if MATLABV == 1
        min_curv = round(prctile(AOI_DEM_curv_meanc.Z(:),5),2);
        max_curv = round(prctile(AOI_DEM_curv_meanc.Z(:),95),2);
    else
        min_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_curv_meanc.Z(:),5)));
        max_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_curv_meanc.Z(:),95)));
    end
    imageschs(AOI_DEM, AOI_DEM_curv_meanc, 'caxis', ...
        [min_curv max_curv])
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: Mean Curvature ', ...
        DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14), grid;
    
    subplot(2,2,4,'align')
    if MATLABV == 1
        min_curv = round(prctile(AOI_DEM_rel_2.Z(:),5),2);
        max_curv = round(prctile(AOI_DEM_rel_2.Z(:),95),2);
    else
        min_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_rel_2.Z(:),5)));
        max_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_rel_2.Z(:),95)));
    end
    imageschs(AOI_DEM, AOI_DEM_rel_2, 'caxis', ...
        [min_curv max_curv])
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: %d-m radius local relief ', ...
        DEM_basename_no_underscore, relief_values_m(2));
    title(title_string, 'Fontsize', 14), grid;
    export_fig(dem_plot_fname,quality_flag,'-pdf');
end
close all

figure
dem_plot_fname = sprintf('maps/%s_Curvature_comp.pdf', DEM_basename_nodir);
if exist(dem_plot_fname, 'file') ~= 2
    fprintf(1,'\tgenerating curvature-comparison figure\n');
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', PaperType_size);
    subplot(2,2,1,'align')
    if MATLABV == 1
        min_curv = round(prctile(AOI_DEM_wiener_curv.Z(:),5),2);
        max_curv = round(prctile(AOI_DEM_wiener_curv.Z(:),95),2);
    else
        min_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_wiener_curv.Z(:),5)));
        max_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_wiener_curv.Z(:),95)));
    end
    imageschs(AOI_DEM, AOI_DEM_wiener_curv, 'caxis', ...
        [min_curv max_curv])
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: Wiener-filtered Curvature ', ...
        DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14), grid;
    
    subplot(2,2,2,'align')
    if MATLABV == 1
        min_curv = round(prctile(AOI_DEM_diffusionf_curv.Z(:),5),2);
        max_curv = round(prctile(AOI_DEM_diffusionf_curv.Z(:),95),2);
    else
        min_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_diffusionf_curv.Z(:),5)));
        max_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_diffusionf_curv.Z(:),95)));
    end
    imageschs(AOI_DEM, AOI_DEM_diffusionf_curv, 'caxis', ...
        [min_curv max_curv])
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: Diffusion-filtered Curvature ', ...
        DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14), grid;
    
    subplot(2,2,3,'align')
    data2plot= AOI_DEM_diffusionf_curv - AOI_DEM_wiener_curv;
    if MATLABV == 1
        min_curv = round(prctile(data2plot.Z(:),5),2);
        max_curv = round(prctile(data2plot.Z(:),95),2);
    else
        min_curv = str2num(sprintf('%.2f',prctile(data2plot.Z(:),5)));
        max_curv = str2num(sprintf('%.2f',prctile(data2plot.Z(:),95)));
    end
    imageschs(AOI_DEM, data2plot, 'caxis', ...
        [min_curv max_curv])
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: Diffusion-filtered minus Wiener Curvature ', ...
        DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14), grid;
    
    subplot(2,2,4,'align')
    if MATLABV == 1
        min_curv = round(prctile(AOI_DEM_rel_2.Z(:),5),2);
        max_curv = round(prctile(AOI_DEM_rel_2.Z(:),95),2);
    else
        min_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_rel_2.Z(:),5)));
        max_curv = str2num(sprintf('%.2f',prctile(AOI_DEM_rel_2.Z(:),95)));
    end
    imageschs(AOI_DEM, AOI_DEM_rel_2, 'caxis', ...
        [min_curv max_curv])
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: %d-m radius local relief ', ...
        DEM_basename_no_underscore, relief_values_m(2));
    title(title_string, 'Fontsize', 14), grid;
    export_fig(dem_plot_fname,quality_flag,'-pdf');
end
close all

figure
dem_plot_fname = sprintf('maps/%s_DEM_ksn.pdf', DEM_basename_nodir);
if exist(dem_plot_fname, 'file') ~= 2
    fprintf(1,'\tgenerating k_sn figures\n');
    symbolspec_ksn045 = makesymbolspec('line',...
        {'ksn045' [prctile([AOI_STR_MS.ksn045], 5) ...
        prctile([AOI_STR_MS.ksn045], 95)] 'color' jet(6)});
    symbolspec_ks_adj = makesymbolspec('line',...
        {'ks_adj' [prctile([AOI_STR_MS.ks_adj], 5) ...
        prctile([AOI_STR_MS.ks_adj], 95)] 'color' jet(6)});
    clf
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', PaperType_size);
    subplot(2,1,1,'align')
    imageschs(AOI_DEM, AOI_DEM, 'caxis', ...
        [floor(min(AOI_DEM(:))) ceil(max(AOI_DEM(:)))], ...
        'colormap',gray,'colorbar',false)
    hold on
    mapshow(AOI_STR_MS,'SymbolSpec',symbolspec_ksn045, 'Linewidth', 2);
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: K_{sn} -0.45 ', DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14), grid;
    colorbar
    caxis([prctile([AOI_STR_MS.ksn045], 5) prctile([AOI_STR_MS.ksn045], 95)])
    hold off
    
    subplot(2,1,2,'align')
    imageschs(AOI_DEM, AOI_DEM, 'caxis', [floor(min(AOI_DEM(:))) ...
        ceil(max(AOI_DEM(:)))], 'colormap',gray,'colorbar',false)
    hold on
    mapshow(AOI_STR_MS,'SymbolSpec',symbolspec_ks_adj, 'Linewidth', 2);
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: K_{sn} adjusted, k_{sn}: %0.3f ', ...
        DEM_basename_no_underscore, AOI_slopearea.theta);
    title(title_string, 'Fontsize', 14), grid;
    colorbar
    caxis([prctile([AOI_STR_MS.ks_adj], 5) prctile([AOI_STR_MS.ks_adj], 95)])
    hold off
    export_fig(dem_plot_fname,quality_flag,'-pdf');
end
close all

figure
dem_plot_fname = sprintf('plots/%s_ridgecrest_DEM_curv.pdf', DEM_basename_nodir);
if exist(dem_plot_fname, 'file') ~= 2
    fprintf(1,'\tgenerating ridgecrest plots\n');
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', PaperType_size);
    subplot(2,1,1,'align')
    hold on
    for i = 1:length(AOI_ridgecrest_MS)
        if i == 1
            h1 = plot(AOI_ridgecrest_MS(i).distance_norm, AOI_ridgecrest_DEM(i).values_norm, 'Color', [0.5, 0.5, 0.5]);
        else
            plot(AOI_ridgecrest_MS(i).distance_norm, AOI_ridgecrest_DEM(i).values_norm, 'Color', [0.5, 0.5, 0.5]);
        end
    end
    grid on
    x1 = [-1:ridgecrest_stepsize:1];
    x1_norm = (x1 - min(x1)) ./ (max(x1) - min(x1));
    %normalize
    y1_parabola = -x1.^2; y1_parabola_norm = y1_parabola - min(y1_parabola) ./ (max(y1_parabola) - min(y1_parabola));
    y1_cosh = cosh(x1.*2); y1_cosh_norm = (y1_cosh - min(y1_cosh)) ./ (max(y1_cosh) - min(y1_cosh));
    y1_cosh_norm = -y1_cosh_norm; y1_cosh2_norm = y1_cosh_norm + 1;
    y1_cosh = cosh(x1.*4); y1_cosh_norm = (y1_cosh - min(y1_cosh)) ./ (max(y1_cosh) - min(y1_cosh));
    y1_cosh_norm = -y1_cosh_norm; y1_cosh4_norm = y1_cosh_norm + 1;
    h2 = plot(x1_norm, y1_parabola_norm, 'Linewidth', 2, 'Color', 'k');
    h3 = plot(x1_norm, y1_cosh2_norm, 'Linewidth', 2, 'Color', 'm');
    h4 = plot(x1_norm, y1_cosh4_norm, 'Linewidth', 2, 'Color', 'c');
    if size(nanmean(AOI_ridgecrest_MS_yi_all'),1) > 1 || size(nanmean(AOI_ridgecrest_MS_yi_all'),2) > 1
        h5 = plot(x1_norm, nanmean(AOI_ridgecrest_MS_yi_all'), 'r-', 'Linewidth', 2);
    else
        h5 = plot(x1_norm, AOI_ridgecrest_MS_yi_all, 'r-', 'Linewidth', 2);
    end        
    legend([h1 h2 h3 h4 h5],{'ridgecrest profiles', 'parabola x^2', 'cos hyperb * 2', 'cos hyperb * 4', 'mean of all ridgecrest profiles'});
    xlabel('normalized distance (by max. distance)', 'Fontsize', 12);
    ylabel('normalized elevation (by max. elevation)', 'Fontsize', 12);
    title_string = sprintf('%s: Ridgecrest profiles, normalized distance and elevation', ...
    DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14);

    %filter curvature values
    windowWidth = 21; % filter length - 11-point moving average, increase for smoother output
    kernel = ones(windowWidth,1) / windowWidth;
    subplot(2,1,2,'align')
    hold on
    %     for i = 1:length(AOI_ridgecrest_MS)
    %         plot(AOI_ridgecrest_MS(i).distance_norm, filter(kernel, 1, AOI_ridgecrest_diffusionf_curv(i).values), 'Color', [0.5, 0.5, 0.5]);
    %     end
    %         plot(repmat(x1_norm, [size(AOI_ridgecrest_MS_deltay_all,1) 1]), AOI_ridgecrest_MS_deltay_all, 'Color', [0.5, 0.5, 0.5]);
    for i = 1:size(AOI_ridgecrest_MS_Dy_all,2)
        if i == 1
            h1 = plot(x1_norm, AOI_ridgecrest_MS_Dy_all(:,i), 'Color', [0.5, 0.5, 0.5]);
        else
            plot(x1_norm, AOI_ridgecrest_MS_Dy_all(:,i), 'Color', [0.5, 0.5, 0.5]);
        end
    end
    if size(nanmean(AOI_ridgecrest_MS_Dy_all'), 1) > 1 || size(nanmean(AOI_ridgecrest_MS_Dy_all'), 2) > 1
        h2 = plot(x1_norm, nanmean(AOI_ridgecrest_MS_Dy_all'), 'k-', 'Linewidth', 2);
        h3 = plot([min(x1_norm) max(x1_norm)], [nanmean(AOI_ridgecrest_MS_Dy_all(:)) nanmean(AOI_ridgecrest_MS_Dy_all(:))], 'r-', 'Linewidth', 2);
    else
        h2 = plot(x1_norm, AOI_ridgecrest_MS_Dy_all, 'k-', 'Linewidth', 2);
        h3 = plot([min(x1_norm) max(x1_norm)], [nanmean(AOI_ridgecrest_MS_Dy_all(:)) nanmean(AOI_ridgecrest_MS_Dy_all(:))], 'r-', 'Linewidth', 2);
    end
    legend([h1 h2 h3],{'\Delta ridgecrest profiles', ...
        'mean of all \Delta ridgecrest profiles', 'mean \Delta ridgecrest profiles'});
    grid on
    xlabel('normalized distance (by max. distance)', 'Fontsize', 12);
    ylabel('\Delta parabola - norm. elevation', 'Fontsize', 12);
    title_string = sprintf('%s: \Delta Ridgecrest profiles, parabola minus normalized elevation', ...
    DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14);
    export_fig(dem_plot_fname,quality_flag,'-pdf');
end
close all

figure
dem_plot_fname = sprintf('maps/%s_DEM_ridgecrest.pdf', DEM_basename_nodir);
if exist(dem_plot_fname, 'file') ~= 2
    fprintf(1,'\tgenerating ridgecrest maps\n');
    clf
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', PaperType_size);
    subplot(2,1,1,'align')
    imageschs(AOI_DEM, AOI_DEM, 'caxis', ...
        [floor(min(AOI_DEM(:))) ceil(max(AOI_DEM(:)))], ...
        'colormap',gray,'colorbar',false)
    hold on
    grid on
    for i = 1:length(AOI_ridgecrest_MS)
        scatter(AOI_ridgecrest_MS(i).x(AOI_ridgecrest_MS(i).Dy_mean_lg_1std), ...
            AOI_ridgecrest_MS(i).y(AOI_ridgecrest_MS(i).Dy_mean_lg_1std), ...
            abs(AOI_ridgecrest_MS(i).Dy_mean(AOI_ridgecrest_MS(i).Dyi_mean_lg_1std).*100), ...
            AOI_ridgecrest_MS(i).Dy_mean(AOI_ridgecrest_MS(i).Dyi_mean_lg_1std))
    end
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: 1-s ridgecrest anom. (mean) ', DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14);
    h1 = colorbar; h1.Label.String = '\Delta norm. elevation'; h1.Label.FontSize = 12;
    %caxis([prctile([AOI_STR_MS.ksn045], 5) prctile([AOI_STR_MS.ksn045], 95)])
    hold off
    
    subplot(2,1,2,'align')
    imageschs(AOI_DEM, AOI_DEM, 'caxis', [floor(min(AOI_DEM(:))) ...
        ceil(max(AOI_DEM(:)))], 'colormap',gray,'colorbar',false)
    hold on
    grid on
    for i = 1:length(AOI_ridgecrest_MS)
        scatter(AOI_ridgecrest_MS(i).x(AOI_ridgecrest_MS(i).Dy_mean_lg_2std), ...
            AOI_ridgecrest_MS(i).y(AOI_ridgecrest_MS(i).Dy_mean_lg_2std), ...
            abs(AOI_ridgecrest_MS(i).Dy_mean(AOI_ridgecrest_MS(i).Dyi_mean_lg_2std).*100), ...
            AOI_ridgecrest_MS(i).Dy_mean(AOI_ridgecrest_MS(i).Dyi_mean_lg_2std))
    end
    ylabel('UTM-Northing (m)', 'Fontsize', 12);
    xlabel('UTM-Easting (m)', 'Fontsize', 12);
    title_string = sprintf('%s: 2-s ridgecrest anom. (mean) ', DEM_basename_no_underscore);
    title(title_string, 'Fontsize', 14);
    h2 = colorbar; h2.Label.String = '\Delta norm. elevation'; h2.Label.FontSize = 12;
    %caxis([prctile([AOI_STR_MS.ks_adj], 5) prctile([AOI_STR_MS.ks_adj], 95)])
    hold off
    export_fig(dem_plot_fname,quality_flag,'-pdf');
end
close all
