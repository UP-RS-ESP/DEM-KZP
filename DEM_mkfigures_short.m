function DEM_mkfigures_short(MATLABV, PaperType_size, quality_flag, DEM_basename_nodir, ...
    AOI_DEM, DEM_basename_no_underscore, AOI_rivers_STR_area2, AOI_DEM_gradient8, ...
    AOI_STR_MS, AOI_slopearea)

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
