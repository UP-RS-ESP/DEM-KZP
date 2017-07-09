%% Make figure for current tributary {j} Save figure in directory

basin_num = sprintf('%02d',i); % get current basin number
trib_num = sprintf('%02d',j); % get current tributary number
dem_plot_fname = sprintf('%s%s%s_LRP_b_%s_t_%s.pdf', KZP_parameters.KZP_plots_dirname, ...
    KZP_parameters.dir_sep, KZP_parameters.DEM_basename_nodir, basin_num, trib_num);
if exist(dem_plot_fname, 'file') ~= 2
    figure; set(gcf,'Visible', 'off');
    fprintf(1,'\tgenerating LRP tributary figure for basin %s and tributary %s\n', basin_num, trib_num);
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', KZP_parameters.PaperType_size);
    subplot(3,1,1,'align')
    
    plot(current_distance, current_elev,'k-', 'Linewidth', 2)
    hold on
    plot(kp_distance,kp_elev,'r.','MarkerSize',32)
    plot(kp_distance_bases,kp_elevation_bases,'bs','MarkerSize',32)
    grid on
    ylabel('Elevation (m)', 'Fontsize', 12)
    xlabel('Distance (m)', 'Fontsize', 12)
    title_s = sprintf('%s: LRP for basin %s and tributary %s\n', KZP_parameters.DEM_basename_no_underscore, basin_num, trib_num);
    title(title_s, 'Fontsize', 16)
    
    subplot(3,1,2)
    plot(current_chi, current_elev,'k-', 'Linewidth', 2)
    hold on
    plot(kp_chi,kp_elev,'r.','MarkerSize',32)
    plot(kp_chi_bases,kp_elevation_bases,'bs','MarkerSize',32)
    grid on
    ylabel('Elevation (m)', 'Fontsize', 12)
    xlabel('chi (distance scaled by m/n and A)', 'Fontsize', 12)
    title_s = sprintf('%s: Chi plot for basin %s and tributary %s\n', KZP_parameters.DEM_basename_no_underscore, basin_num, trib_num);
    title(title_s, 'Fontsize', 16)
    
    subplot(3,1,3)
    hold on
    plot(current_chi_resampled,smooth_detrended_elev_current,'color',[0 0 0.5], 'LineWidth',2)
    plot(current_chi_resampled, detrended_elev_current,'k-')
    plot(current_chi_resampled(sig_kps_cells_final),smooth_detrended_elev_current(sig_kps_cells_final),'r.','MarkerSize',32)
    plot(current_chi_resampled(sig_kps_bases_final),smooth_detrended_elev_current(sig_kps_bases_final),'bs','MarkerSize',32)
    grid on
    ylabel('Detrended elevation (m)', 'Fontsize', 12)
    xlabel('chi (distance scaled by m/n and A)', 'Fontsize', 12)
    title_s = sprintf('%s: Detrended Chi plot for basin %s and tributary %s\n', KZP_parameters.DEM_basename_no_underscore, basin_num, trib_num);
    title(title_s, 'Fontsize', 16)
    
    if exist('export_fig') == 2
        export_fig(dem_plot_fname,KZP_parameters.quality_flag,'-pdf');
    else
        saveas(gcf,dem_plot_fname, 'pdf');
    end
end