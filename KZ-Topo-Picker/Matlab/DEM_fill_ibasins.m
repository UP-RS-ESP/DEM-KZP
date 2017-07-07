function [DEM2, DEM_SINKS_list] = DEM_fill_ibasins(DEM, basic_fill_value, area_threshold_a, area_threshold_b, fill_elevation_steps, dbasin_flag, figure_flag, min_drainage_area_to_process, PaperType_size, quality_flag, polygonize_cmd)
%function [DEM2, DEM_SINKS_list] = DEM_fill_ibasins(DEM, basic_fill_value, area_threshold_a, area_threshold_b, fill_elevation_steps, dbasin_flag, figure_flag, min_drainage_area_to_process, PaperType_size, quality_flag, polygonize_cmd)
%
% Fill internally drained basins in a DEM (e.g., Puna Plateau, Tibet)
%
% Bodo Bookhagen, April 2016
%
% DEM                 - input digital elevation model
% basic_fill_value    - first fill value, will always be applied to DEM
% area_threshold_a    - drainage area threshold for filled-basin removal,
%                       usually 0.1e6 (0.1 km^2)
% area_threshold_b    - drainage are threshold for basin processing (usually
%                       1e6 (1km^2)). All closed basins larger this area will
%                       be processed.
% fill_elevation_step - vector giving minimum and maximum elevation to be
%                       filled and the stepsize (e.g., [10:10:100])
% dbasin_flag         - drainage basin flag, output to geotiff (yes = true)
% figure_flag         - print figures (figure_flag = true)
% min_drainage_area_to_process - minimum drainage area for plots
%
%
if figure_flag == true && exist('DEM_nofill.pdf', 'file') ~= 2
    AOI_DEM = DEM;
    AOI_FD = FLOWobj(AOI_DEM,'preprocess','none', 'internaldrainage', true, 'mex', true);
    AOI_FAC = flowacc(AOI_FD);
    [AOI_dbasins, AOI_dbasins_outlet] = drainagebasins(AOI_FD);
    
    AOI_resolution = AOI_DEM.refmat(2,1);
    % get the DEM resolution, used to convert area threshold to meters squared
    minApix = area_threshold_b/(AOI_resolution*AOI_resolution);
    minApix = ceil(minApix);  % convert area threshold to meters specified above
    AOI_FAC_w = AOI_FAC > minApix;  % masks flow accum grid (above threshold)
    
    AOI_rivers_w = AOI_FAC.*(AOI_FAC.cellsize.^2) > ...
        min_drainage_area_to_process;
    AOI_STR_w = STREAMobj(AOI_FD,AOI_FAC_w);
    AOI_rivers_STR = STREAMobj(AOI_FD,AOI_rivers_w);
    
    figure;
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', PaperType_size);
    subplot(2,2,1,'align')
    imageschs(AOI_DEM, AOI_DEM,'caxis', [floor(min(AOI_DEM(:))) ...
        ceil(max(AOI_DEM(:)))])
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('DEM no fill', 'Fontsize', 14);
    
    subplot(2,2,2,'align')
    imageschs(AOI_DEM, shufflelabel(AOI_dbasins))
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('Basins no fill', 'Fontsize', 14);
    
    subplot(2,2,3,'align')
    imageschs(AOI_DEM, AOI_DEM-AOI_DEM, 'caxis', [0 10])
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('\Delta Elevation', 'Fontsize', 14);
    
    subplot(2,2,4,'align')
    imageschs(AOI_DEM, AOI_DEM-AOI_DEM)
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('\Delta Elevation', 'Fontsize', 14);
    export_fig('DEM_nofill.pdf',quality_flag,'-pdf');
    close
    clear AOI_DEM AOI_rivers_STR AOI_dbasins AOI_STR_w AOI_dbasins_outlet AOI_FAC AOI_FD AOI_FAC_w AOI_resolution AOI_rivers_w minApix
end

if figure_flag == true && exist('DEM_minfill.pdf', 'file') ~= 2
    AOI_DEM = fillsinks(DEM,basic_fill_value);
    AOI_FD = FLOWobj(AOI_DEM,'preprocess','none', 'internaldrainage', true, 'mex', true);
    AOI_FAC = flowacc(AOI_FD);
    [AOI_dbasins, AOI_dbasins_outlet] = drainagebasins(AOI_FD);
    
    AOI_resolution = AOI_DEM.refmat(2,1);
    % get the DEM resolution, used to convert area threshold to meters squared
    minApix = area_threshold_b/(AOI_resolution*AOI_resolution);
    minApix = ceil(minApix);  % convert area threshold to meters specified above
    AOI_FAC_w = AOI_FAC > minApix;  % masks flow accum grid (above threshold)
    
    AOI_rivers_w = AOI_FAC.*(AOI_FAC.cellsize.^2) > ...
        min_drainage_area_to_process;
    AOI_STR_w = STREAMobj(AOI_FD,AOI_FAC_w);
    AOI_rivers_STR = STREAMobj(AOI_FD,AOI_rivers_w);
    
    figure;
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', PaperType_size);
    subplot(2,2,1,'align')
    imageschs(AOI_DEM, AOI_DEM,'caxis', [floor(min(AOI_DEM(:))) ...
        ceil(max(AOI_DEM(:)))])
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('DEM min. fill', 'Fontsize', 14);
    
    subplot(2,2,2,'align')
    imageschs(AOI_DEM, shufflelabel(AOI_dbasins))
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('Basins min. fill', 'Fontsize', 14);
    
    subplot(2,2,3,'align')
    imageschs(AOI_DEM, AOI_DEM-AOI_DEM, 'caxis', [0 10])
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('\Delta Elevation', 'Fontsize', 14);
    
    subplot(2,2,4,'align')
    imageschs(AOI_DEM, AOI_DEM-AOI_DEM)
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('\Delta Elevation', 'Fontsize', 14);
    export_fig('DEM_minfill.pdf',quality_flag,'-pdf');
    close
    clear AOI_DEM AOI_rivers_STR AOI_dbasins AOI_STR_w AOI_dbasins_outlet AOI_FAC AOI_FD AOI_FAC_w AOI_resolution AOI_rivers_w minApix
end

% calculate filled basins
if figure_flag == true && exist('DEM_carve.pdf', 'file') ~= 2
    AOI_DEM = DEM;
    AOI_FD = FLOWobj(AOI_DEM,'preprocess','carve', 'internaldrainage', true, 'mex', true);
    AOI_FAC = flowacc(AOI_FD);
    [AOI_dbasins, AOI_dbasins_outlet] = drainagebasins(AOI_FD);
    
    AOI_resolution = AOI_DEM.refmat(2,1);
    % get the DEM resolution, used to convert area threshold to meters squared
    minApix = area_threshold_b/(AOI_resolution*AOI_resolution);
    minApix = ceil(minApix);  % convert area threshold to meters specified above
    AOI_FAC_w = AOI_FAC > minApix;  % masks flow accum grid (above threshold)
    
    AOI_rivers_w = AOI_FAC.*(AOI_FAC.cellsize.^2) > ...
        min_drainage_area_to_process;
    AOI_STR_w = STREAMobj(AOI_FD,AOI_FAC_w);
    AOI_rivers_STR = STREAMobj(AOI_FD,AOI_rivers_w);
    
    figure;
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', PaperType_size);
    subplot(2,2,1,'align')
    imageschs(AOI_DEM, AOI_DEM,'caxis', [floor(min(AOI_DEM(:))) ...
        ceil(max(AOI_DEM(:)))])
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('DEM carved', 'Fontsize', 14);
    
    subplot(2,2,2,'align')
    imageschs(AOI_DEM, shufflelabel(AOI_dbasins))
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('Basins carved', 'Fontsize', 14);
    
    subplot(2,2,3,'align')
    imageschs(AOI_DEM, AOI_DEM-DEM, 'caxis', [0 10])
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('\Delta Elevation', 'Fontsize', 14);
    
    subplot(2,2,4,'align')
    imageschs(AOI_DEM, AOI_DEM-DEM)
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('\Delta Elevation', 'Fontsize', 14);
    export_fig('DEM_carve.pdf',quality_flag,'-pdf');
    close
    clear AOI_DEM AOI_rivers_STR AOI_dbasins AOI_STR_w AOI_dbasins_outlet AOI_FAC AOI_FD AOI_DEM_FIL AOI_FAC_w AOI_resolution AOI_rivers_w minApix
end

if figure_flag == true && exist('DEM_fill.pdf', 'file') ~= 2
    AOI_DEM = DEM;
    AOI_FD = FLOWobj(AOI_DEM,'preprocess','fill', 'internaldrainage', true, 'mex', true);
    AOI_FAC = flowacc(AOI_FD);
    [AOI_dbasins, AOI_dbasins_outlet] = drainagebasins(AOI_FD);
    
    AOI_resolution = AOI_DEM.refmat(2,1);
    % get the DEM resolution, used to convert area threshold to meters squared
    minApix = area_threshold_b/(AOI_resolution*AOI_resolution);
    minApix = ceil(minApix);  % convert area threshold to meters specified above
    AOI_FAC_w = AOI_FAC > minApix;  % masks flow accum grid (above threshold)
    
    AOI_rivers_w = AOI_FAC.*(AOI_FAC.cellsize.^2) > ...
        min_drainage_area_to_process;
    AOI_STR_w = STREAMobj(AOI_FD,AOI_FAC_w);
    AOI_rivers_STR = STREAMobj(AOI_FD,AOI_rivers_w);
    
    figure;
    set(gcf,'units','normalized','position',[0 0 1 1]);
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperType', PaperType_size);
    subplot(2,2,1,'align')
    imageschs(AOI_DEM, AOI_DEM,'caxis', [floor(min(AOI_DEM(:))) ...
        ceil(max(AOI_DEM(:)))])
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('DEM filled', 'Fontsize', 14);
    
    subplot(2,2,2,'align')
    imageschs(AOI_DEM, shufflelabel(AOI_dbasins))
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('Basins filled', 'Fontsize', 14);
    
    subplot(2,2,3,'align')
    imageschs(AOI_DEM, AOI_DEM-DEM, 'caxis', [0 10])
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('\Delta Elevation', 'Fontsize', 14);
    
    subplot(2,2,4,'align')
    imageschs(AOI_DEM, AOI_DEM-DEM)
    hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);
    title('\Delta Elevation', 'Fontsize', 14);
    export_fig('DEM_fill.pdf',quality_flag,'-pdf');
    close
    clear AOI_DEM AOI_rivers_STR AOI_dbasins AOI_STR_w AOI_dbasins_outlet AOI_FAC AOI_FD AOI_DEM_FIL AOI_FAC_w AOI_resolution AOI_rivers_w minApix
end

fill_values = fill_elevation_steps;
DEM2 = fillsinks(DEM,basic_fill_value);
DEM_DIF10 = DEM.Z-DEM2.Z;
idx1 = find(DEM_DIF10 > 0);
idx0 = find(isnan(DEM_DIF10));
DEM_DIF_bw = DEM_DIF10;
DEM_DIF_bw(idx1) = 1;
DEM_DIF_bw(idx0) = 0;
DEM_DIF_bw = logical(DEM_DIF_bw);
DEM_DIF_bw_basins = logical(bwlabel(DEM_DIF_bw));
DEM_DIF_bw_basins_stats = regionprops(DEM_DIF_bw_basins, DEM_DIF_bw_basins, 'Area', 'PixelIdxList', 'MeanIntensity');
%DEM_DIF_bw_basins_topo_stats = regionprops(DEM_DIF_bw_basins, DEM2.Z, 'Area', 'MeanIntensity', 'MinIntensity', 'MaxIntensity');
    
%remove NaN, 1st element
idx2remove = find([DEM_DIF_bw_basins_stats.MeanIntensity] == 0);
if ~isempty(idx2remove)
    DEM_DIF_bw_basins_stats(idx2remove) = [];
end
    
%fill only cells with area < 0.1 km^2
area_threshold_pixels = ceil(area_threshold_a / (DEM2.cellsize^2));
idx = find([DEM_DIF_bw_basins_stats.Area] > area_threshold_pixels);
%this part could be optimized, for loop is not efficient to combine
%list of sinks
counter = 1;
for j = 1:length(idx)
    sink_list_length = length([DEM_DIF_bw_basins_stats(idx(j)).PixelIdxList]);
    sink_list(counter:counter+sink_list_length-1) = [DEM_DIF_bw_basins_stats(idx(j)).PixelIdxList]';
    counter = counter + sink_list_length;
end
DEM_DIF_bw_basins_lg_threshold = DEM_DIF_bw_basins;
DEM_SINKS = logical(zeros(size(DEM2.Z),'int8'));
if exist('sink_list', 'var') ~= 0
    DEM_DIF_bw_basins_lg_threshold(~sink_list) = 0;
    DEM_SINKS(sink_list) = 1;
end
idx = find([DEM_DIF_bw_basins_stats.Area] < area_threshold_pixels);
DEM_DIF_bw_basins_stats(idx) = []; 
DEM_SINKS_gr = DEM;
DEM_SINKS_gr.Z = DEM_SINKS;
clear DEM_SINKS
    
AOI_FD = FLOWobj(DEM2,'preprocess','carve', 'sinks', DEM_SINKS_gr, 'internaldrainage', true, 'mex', true);
AOI_FAC = flowacc(AOI_FD);
[AOI_dbasins_bfill, AOI_dbasins_bfill_outlet] = drainagebasins(AOI_FD);
save -v7.3 DEM_basicfill AOI_FD AOI_FAC AOI_dbasins_bfill AOI_dbasins_bfill_outlet DEM_SINKS_gr
dbasins_fname = sprintf('DEM_bfill%02d.tif', basic_fill_value);
dbasins_out_fname_shp = sprintf('DEM_bfill%02d_poly.shp', basic_fill_value);
GRIDobj2geotiff(AOI_dbasins_bfill, dbasins_fname)
eval([polygonize_cmd, ' -8 -mask ', dbasins_fname, ' ', dbasins_fname, ' -b 1 -f "ESRI Shapefile" ', dbasins_out_fname_shp, ' DBasin DBasin-ID &']);

warning off
DEM_SINKS_list = logical(zeros([size(DEM2.Z) length(fill_values)], 'int8'));

for i = 1:length(fill_values)
    fprintf(1,'at fill_value = %d (%d/%d)\n', fill_values(i), i, length(fill_values));
    tic;
    AOI_DEM_FIL = fillsinks(DEM2,fill_values(i));
    
    DEM_DIF10 = AOI_DEM_FIL.Z-DEM2.Z;
    idx1 = find(DEM_DIF10 > 0);
    idx0 = find(isnan(DEM_DIF10));
    DEM_DIF_bw = DEM_DIF10;
    DEM_DIF_bw(idx1) = 1;
    DEM_DIF_bw(idx0) = 0;
    DEM_DIF_bw = logical(DEM_DIF_bw);
    DEM_DIF_bw_basins = logical(bwlabel(DEM_DIF_bw));
    DEM_DIF_bw_basins_stats = regionprops(DEM_DIF_bw_basins, DEM_DIF_bw_basins, 'Area', 'PixelIdxList', 'MeanIntensity');
    DEM_DIF_bw_basins_topo_stats = regionprops(DEM_DIF_bw_basins, DEM2.Z, 'Area', 'MeanIntensity', 'MinIntensity', 'MaxIntensity');
    
    %remove NaN, 1st element
    idx2remove = find([DEM_DIF_bw_basins_stats.MeanIntensity] == 0);
    if ~isempty(idx2remove)
        DEM_DIF_bw_basins_stats(idx2remove) = [];
        DEM_DIF_bw_basins_topo_stats(idx2remove) = [];
    end
        
    %DEM_DIF_bw_basins_stats(1) = [];
    %fill only cells with area < 0.1 km^2
    area_threshold_pixels = ceil(area_threshold_a / (DEM2.cellsize^2));
    idx = find([DEM_DIF_bw_basins_stats.Area] > area_threshold_pixels);
    %this part could be optimized, for loop is not efficient to combine
    %list of sinks
    counter = 1;
    for j = 1:length(idx)
        sink_list_length = length([DEM_DIF_bw_basins_stats(idx(j)).PixelIdxList]);
        sink_list(counter:counter+sink_list_length-1) = [DEM_DIF_bw_basins_stats(idx(j)).PixelIdxList]';
        counter = counter + sink_list_length;
    end
    DEM_DIF_bw_basins_lg_threshold = DEM_DIF_bw_basins;
    DEM_SINKS = logical(zeros(size(DEM2.Z),'int8'));
    if exist('sink_list', 'var') ~= 0
        DEM_DIF_bw_basins_lg_threshold(~sink_list) = 0;
        DEM_SINKS(sink_list) = 1;
    end
    idx = find([DEM_DIF_bw_basins_stats.Area] < area_threshold_pixels);
    if length(idx) > 1
        DEM_DIF_bw_basins_stats(idx) = []; DEM_DIF_bw_basins_topo_stats(idx) = [];
    end
    DEM_SINKS_list(:,:,i) = DEM_SINKS;
    
    %merge original basins and SINK mask
    %DEM_dbasins_SINKS = uint32(double(AOI_dbasins_org.Z) .* DEM_SINKS);
    %plot
    %imageschs(DEM2, bwlabel(DEM_dbasins_SINKS), 'nancolor', [1 1 1]), colorbar
    DEM_SINKS_gr = AOI_DEM_FIL;
    DEM_SINKS_gr.Z = DEM_SINKS;
    clear DEM_SINKS

    AOI_FD = FLOWobj(AOI_DEM_FIL,'preprocess','carve', 'sinks', DEM_SINKS_gr, 'internaldrainage', true, 'mex', true);
    AOI_FAC = flowacc(AOI_FD);
    [AOI_dbasins_fill, AOI_dbasins_fill_outlet] = drainagebasins(AOI_FD);
    mat_fname = sprintf('DEM_fill%02d.mat', fill_values(i));
    save(mat_fname, 'AOI_DEM_FIL', 'AOI_FD', 'AOI_FAC', 'AOI_dbasins_fill', 'AOI_dbasins_fill_outlet', 'DEM_SINKS_gr', '-v7.3');
    
    dbasins_fname = sprintf('DEM_fill%02d.tif', fill_values(i));
    dbasins_out_fname_shp = sprintf('DEM_fill%02d_poly.shp', fill_values(i));
    dbasins_out_fname_kml = sprintf('DEM_fill%02d_poly.kml', fill_values(i));
    GRIDobj2geotiff(AOI_dbasins_fill, dbasins_fname)
    eval([polygonize_cmd, ' -8 -mask ', dbasins_fname, ' ', dbasins_fname, ' -b 1 -f "ESRI Shapefile" ', dbasins_out_fname_shp, ' DBasin DBasin-ID &']);
%    eval([ogr2ogr_cmd, ' ', dbasins_out_fname_kml, ' ', dbasins_out_fname_shp, ' -f "KML"']);
    toc
    fname = sprintf('DEM_it%02d_fill%02d.pdf', i, fill_values(i));
%    foo=[6.439397950293089e+05 7.737569070982749e+05 7.186685386932234e+06 7.349253188656373e+06];
    if figure_flag == true && exist(fname, 'file') ~= 2 || dbasin_flag == true
        AOI_FD = FLOWobj(DEM2,'preprocess','none', 'internaldrainage', true);
        AOI_FAC = flowacc(AOI_FD);
        [AOI_dbasins, AOI_dbasins_outlet] = drainagebasins(AOI_FD);
        dbasins_fname = sprintf('Basins_it%02d_fill%02d.tif', i, fill_values(i));
        
        if dbasin_flag == true
            GRIDobj2geotiff(AOI_dbasins, dbasins_fname);
        end
        AOI_resolution = DEM2.refmat(2,1);
        % get the DEM resolution, used to convert area threshold to meters squared
        minApix = area_threshold_b/(AOI_resolution*AOI_resolution);
        minApix = ceil(minApix);  % convert area threshold to meters specified above
        AOI_FAC_w = AOI_FAC > minApix;  % masks flow accum grid (above threshold)
        
        AOI_rivers_w = AOI_FAC.*(AOI_FAC.cellsize.^2) > ...
            min_drainage_area_to_process;
        AOI_STR_w = STREAMobj(AOI_FD,AOI_FAC_w);
        AOI_rivers_STR = STREAMobj(AOI_FD,AOI_rivers_w);
        
        if figure_flag == true && exist(fname, 'file') ~= 2
            figure;
            set(gcf,'units','normalized','position',[0 0 1 1]);
            set(gcf, 'PaperOrientation', 'landscape');
            set(gcf, 'PaperType', PaperType_size);
            subplot(2,2,1,'align')
            imageschs(DEM2, DEM2,'caxis', [floor(min(DEM2(:))) ...
                ceil(max(DEM2(:)))])
            hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1); axis(foo)
            s=sprintf('%d iteration, fillvalue = %d, DEM', i, fill_values(i));
            title(s, 'Fontsize', 14);
            
            subplot(2,2,2,'align')
            imageschs(DEM2, shufflelabel(AOI_dbasins))
            hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1); axis(foo)
            s=sprintf('%d iteration, fillvalue = %d, DEM Basins', i, fill_values(i));
            title(s, 'Fontsize', 14);
            
            if numel(DEM_DIF_bw_basins_stats) > 0
                subplot(2,2,3,'align')
                imageschs(DEM2, DEM_DIF_bw_basins)
                hold; plot(AOI_rivers_STR, 'k', 'Linewidth', 1);axis(foo)
                s=sprintf('%d iteration, nr. of basins with NaN = %d', i, numel(DEM_DIF_bw_basins_stats));
                title(s, 'Fontsize', 14);
            end
            subplot(2,2,4,'align')
            plot([DEM_DIF_bw_basins_topo_stats.Area]*(AOI_FAC.cellsize.^2), [DEM_DIF_bw_basins_topo_stats.MeanIntensity], '+')
            xlabel('Area (m^2)'), ylabel('elevation (m)'), grid
            s=sprintf('%d iteration, fillvalue = %d, D Elevation', i, fill_values(i));
            title(s, 'Fontsize', 14);
%            export_fig(fname,quality_flag,'-pdf');
%            close
        end
        clear AOI_DEM AOI_rivers_STR AOI_dbasins AOI_STR_w AOI_dbasins_outlet AOI_FAC AOI_FD AOI_FAC_w AOI_resolution AOI_rivers_w minApix
    end
end
