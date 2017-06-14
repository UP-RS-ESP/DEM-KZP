%% (4) Generate figures for each basin
% First we'll make the chiplots, longitudinal profiles, and the hillshades,
% then we'll add the knickpoints to them
% make chiplot for each individual basin

if KZP_parameters.show_figs == 1 || KZP_parameters.show_figs == 2
fprintf(1,'KZP identifying knickzones step 4 of 4: Making figures for all basins.\nAt basin # of %d: ', number_of_basins);
    for i = 1:number_of_basins2
        fprintf('%d, ', i);
        knickpoints_plot_fname = sprintf('plots/%s_basin_%d_knickpoints.pdf', KZP_parameters.DEM_basename_nodir, i);
        if exist(knickpoints_plot_fname, 'file') == 2
            continue
        end
        % turn the stream name into a text that can be pasted to the chiplot
        current_stream_id = sprintf('Stream # %d',i);
        str = sprintf('%s', current_stream_id);
        
        % this is a chi-plot of all of the basin streams to be analyzed. we will add
        % kps later
        
        elev_current_stream = elev_stored_trimmed_basin{i} ; % get the elevation data for current basin from (elev_stored_trimmed_basin) <- contains all basins
        chi_current_stream = chi_stored_trimmed_basin{i} ;% get the chi data for current basin from (elev_stored_trimmed_basin) <- contains all basins
        
        %prepare variables for knickpoint-plotting
        sig_kps_chi_all_tribs_all_kp_size_current_basin = kp_chi_stored_lips{i} ; % unpack knickpoints for each stream basin 1-13
        sig_kps_elevation_all_tribs_all_kp_size_current_basin = kp_elev_stored_lips{i};
        sig_kps_easting_all_tribs_all_kp_size_current_basin = kp_easting_stored_lips{i};
        sig_kps_northing_all_tribs_all_kp_size_current_basin = kp_northing_stored_lips{i};
        sig_kps_magnitude_all_tribs_all_kp_size_current_basin = kp_magnitude_matrix_stored{i};
        sig_kps_distance_all_tribs_all_kp_size_current_basin =kp_distance_stored_lips{i};
        troughs_chi_all_tribs_current_basin = kp_chi_stored_bases{i};
        troughs_elevation_all_tribs_current_basin = kp_elev_stored_bases{i};
        troughs_easting_all_tribs_current_basin = kp_easting_stored_bases{i};
        troughs_northing_all_tribs_current_basin = kp_northing_stored_bases{i};
        troughs_distance_all_tribs_current_basin = kp_distance_stored_bases{i};
        
        %trunk stream
        sig_kps_chi_all_tribs_all_kp_size_current_basin_t = kp_chi_stored_lips_t{i} ; % unpack knickpoints for each stream basin 1-13
        sig_kps_elevation_all_tribs_all_kp_size_current_basin_t = kp_elev_stored_lips_t{i};
        sig_kps_easting_all_tribs_all_kp_size_current_basin_t = kp_easting_stored_lips_t{i};
        sig_kps_northing_all_tribs_all_kp_size_current_basin_t = kp_northing_stored_lips_t{i};
        sig_kps_magnitude_all_tribs_all_kp_size_current_basin_t = kp_magnitude_matrix_stored_t{i};
        sig_kps_distance_all_tribs_all_kp_size_current_basin_t = kp_distance_stored_lips_t{i};
        troughs_chi_all_tribs_current_basin_t = kp_chi_stored_bases_t{i};
        troughs_elevation_all_tribs_current_basin_t = kp_elev_stored_bases_t{i};
        troughs_easting_all_tribs_current_basin_t = kp_easting_stored_bases_t{i};
        troughs_northing_all_tribs_current_basin_t = kp_northing_stored_bases_t{i};
        troughs_distance_all_tribs_current_basin_t = kp_distance_stored_bases_t{i};
        
        figure
        set(gcf,'units','normalized','position',[0 0 1 1]);
        set(gcf, 'PaperOrientation', 'landscape');
        set(gcf, 'PaperType', KZP_parameters.PaperType_size);
        %log area vs. log slope plot
        subplot(2,2,1,'align')
        loglog([min(AOI_STR_S_slopearea_dbasins{i}.a) max(AOI_STR_S_slopearea_dbasins{i}.a)], ...
            [max(AOI_STR_S_slopearea_dbasins{i}.g) min(AOI_STR_S_slopearea_dbasins{i}.g)], ...
            '-', 'LineWidth', 2, 'color', [0 0.4470 0.7410])
        hold on
        loglog(AOI_STR_S_slopearea_dbasins{i}.a, AOI_STR_S_slopearea_dbasins{i}.g, ...
            's', 'Markersize', 8, 'color', [0 0.4470 0.7410])
        loglog(AOI_STR_area_trunk_subset{i}, AOI_STR_slope_trunk_subset{i}, '+', 'Markersize', 2, 'color', [0.5 0.5 0.5])
        set(gca,'Children',flipud(get(gca,'Children')));
        xlabel('Drainage Area (m^2)', 'Fontsize', 12);
        ylabel('Slope (m/m)', 'Fontsize', 12);
        title_string = sprintf('DEM: %s, Basin #: %d, Trunk Stream (theta=0.45)', ...
            KZP_parameters.DEM_basename_no_underscore, i);
        title(title_string, 'Fontsize', 16), grid;
        %adding text: ksn, ks, theta, theta_ref, r2 and regression boundaries from SA plots
        s = sprintf('From slope-area analyses (binned-mean) (DA=%3.2f km^2): k_{sn} trunk = %3.1f with theta = %0.2f, k_{s} trunk = %3.1f with theta = %0.2f\nk_{s} trunk = %3.1f with theta from chi plot = %0.2f, k_{s} trunk = %3.1f +/- %3.1f with theta from robust regression = %0.2f +/- %0.2f\nr^2 = %0.2f, rmse=%2.2e, p(k_{s}) = %1.2e, p(theta) = %1.2e', ...
            max(AOI_STR_S_slopearea_dbasins{i}.a)/1e6, AOI_STR_S_slopearea_dbasins_trunk{i}.ks, AOI_STR_S_slopearea_dbasins_trunk{i}.theta, ...
            AOI_STR_S_slopearea_dbasins_trunk_adj{i}.ks, AOI_STR_S_slopearea_dbasins_trunk_adj{i}.theta, ...
            AOI_STR_S_slopearea_dbasins_trunk_chitheta{i}.ks, AOI_STR_S_slopearea_dbasins_trunk_chitheta{i}.theta, ...
            fitresult_logspace.a, abs(fitresult_logspace.a-ci_logspace(1)), fitresult_logspace.b, abs(fitresult_logspace.b-ci_logspace(3)), ...
            gof_logspace{i}.adjrsquare, gof_logspace{i}.rmse, p_logspace_a, p_logspace_b);
        foo=axis;
        text(0.01, 0.9, s, 'Fontsize', 8, 'Units', 'normalized')
        grid on
        legend('all data', 'mean-bins', 'regression', 'Location','southeast', 'FontSize',8,'FontWeight','bold')
        hold off
        
        
        %chi plot
        subplot(2,2,2,'align')
        plot(AOI_STR_S_chiplot{i}.chi, AOI_STR_S_chiplot{i}.elev, 'color',[.5 .5 .5], 'Linewidth', 1)
        hold on
        plot(AOI_STR_S_trunk_chiplot{i}.chi, AOI_STR_S_trunk_chiplot{i}.elev, 'color',[0 0 0], 'Linewidth', 2)
        xlabel('Chi', 'Fontsize', 12);
        ylabel('Elevation', 'Fontsize', 12);
        title_string = sprintf('DEM: %s, Basin #: %d, Trunk Stream (m/n = %0.2f)', ...
            KZP_parameters.DEM_basename_no_underscore, i, AOI_STR_S_trunk_chiplot{i}.mn);
        title(title_string, 'Fontsize', 16);
        s = sprintf('From Chi analysis (DA=%3.2f km^2):\nk_{s} trunk = %3.1f with m/n= %0.2f\nbeta trunk = %3.1f, betaSE = %1.2e, r^2 = %0.2f\nk_{s} (all streams)= %3.1f with m/n= %0.2f, r^2 = %0.2f', ...
            max(AOI_STR_S_slopearea_dbasins{i}.a)/1e6, AOI_STR_S_trunk_chiplot{i}.ks, AOI_STR_S_trunk_chiplot{i}.mn, ...
            AOI_STR_S_trunk_chiplot{i}.beta, AOI_STR_S_trunk_chiplot{i}.betase, AOI_STR_S_trunk_chiplot{i}.R2, ...
            AOI_STR_S_chiplot{i}.ks, AOI_STR_S_chiplot{i}.mn, AOI_STR_S_chiplot{i}.R2);
        foo=axis;
        text(0.55, 0.2, s, 'Fontsize', 10, 'Units', 'normalized');
        grid on
        hold off
        
        %Hillshade with adjusted K_sn and knickpoints
        symbolspec_ksn045 = makesymbolspec('line',...
            {'ksn045' [prctile([AOI_STR_MS.ksn045], 5) prctile([AOI_STR_MS.ksn045], 95)] 'color' jet(6)});
        symbolspec_ks_adj = makesymbolspec('line',...
            {'ks_adj' [prctile([AOI_STR_MS.ks_adj], 5) prctile([AOI_STR_MS.ks_adj], 95)] 'color' jet(6)});
        subplot3 = subplot(2,2,3, 'align');
        imageschs(AOI_DEM, AOI_DEM, 'caxis', [floor(min(AOI_DEM(:))) ceil(max(AOI_DEM(:)))], 'colormap',gray,'colorbar',false)
        hold on
        mapshow(AOI_STR_MS,'SymbolSpec',symbolspec_ksn045);
        ylabel('UTM-Northing (m)', 'Fontsize', 12);
        xlabel('UTM-Easting (m)', 'Fontsize', 12);
        title_string = sprintf('%s: K_{sn} -0.45 ', KZP_parameters.DEM_basename_no_underscore);
        title(title_string, 'Fontsize', 14), grid;
        colorbar
        caxis([prctile([AOI_STR_MS.ksn045], 5) prctile([AOI_STR_MS.ksn045], 95)])
        xrange = max(AOI_STR_dbasins_unique_subset{i}.x) - min(AOI_STR_dbasins_unique_subset{i}.x);
        yrange = max(AOI_STR_dbasins_unique_subset{i}.y) - min(AOI_STR_dbasins_unique_subset{i}.y);
        axis([min(AOI_STR_dbasins_unique_subset{i}.x)-(xrange/25) ...
            max(AOI_STR_dbasins_unique_subset{i}.x)+(xrange/25) ...
            min(AOI_STR_dbasins_unique_subset{i}.y)-(yrange/25) ...
            max(AOI_STR_dbasins_unique_subset{i}.y)+(yrange/25)])
        hold off
        
        % LRP
        subplot(2,2,4,'align')
        h = plotdz(AOI_STR_streams_dbasins_unique{i},AOI_DEM);
        set(h, 'color', [0.5 0.5 0.5], 'Linewidth', 1);
        hold on
        grid on
        h2 = plotdz(AOI_STR_all_streams_trunk{i},AOI_DEM);
        set(h2, 'color', 'k', 'Linewidth', 2);
        hold off
        
        %determine max. knickpointsize for plotting
        if ~isempty(sig_kps_magnitude_all_tribs_all_kp_size_current_basin)
%             foo = (cellfun(@(x) prctile(x(:),95), sig_kps_magnitude_all_tribs_all_kp_size_current_basin,'un',0));
            foo = (cellfun(@(x) max(x(:)), [sig_kps_magnitude_all_tribs_all_kp_size_current_basin],'un',0));
            sig_kps_magnitude_all_tribs_all_kp_size_current_basin_max = ceil(max(cell2mat([foo(:)]))); clear foo
            
            if sig_kps_magnitude_all_tribs_all_kp_size_current_basin_max > KZP_parameters.max_knickpointsize2plot
                %scale knickpoint for plotting
                kp_size2plot = sig_kps_magnitude_all_tribs_all_kp_size_current_basin_max/KZP_parameters.max_knickpointsize2plot;
            elseif sig_kps_magnitude_all_tribs_all_kp_size_current_basin_max < KZP_parameters.max_knickpointsize2plot
                kp_size2plot = sig_kps_magnitude_all_tribs_all_kp_size_current_basin_max/KZP_parameters.max_knickpointsize2plot;
            else
                kp_size2plot = 1;
            end
        else
            kp_size2plot = 1;
        end
        %add knickpoints - first trunk
        if ~isempty(sig_kps_chi_all_tribs_all_kp_size_current_basin_t)
            sig_kps_chi_tribs_current_trib_t = sig_kps_chi_all_tribs_all_kp_size_current_basin_t{1};  % chi
            sig_kps_elevation_tribs_current_trib_t = sig_kps_elevation_all_tribs_all_kp_size_current_basin_t{1};  % elevation
            sig_kps_easting_tribs_current_trib_t = sig_kps_easting_all_tribs_all_kp_size_current_basin_t{1};  % easting
            sig_kps_northing_tribs_current_trib_t = sig_kps_northing_all_tribs_all_kp_size_current_basin_t{1};  % northing
            sig_kps_magnitude_tribs_current_trib_t =  sig_kps_magnitude_all_tribs_all_kp_size_current_basin_t{1}; %kp magnitude
            sig_kps_distance_tribs_current_trib_t =  sig_kps_distance_all_tribs_all_kp_size_current_basin_t{1}; %kp magnitude
            troughs_chi_tribs_current_trib_t = troughs_chi_all_tribs_current_basin_t{1};
            troughs_elevation_tribs_current_trib_t = troughs_elevation_all_tribs_current_basin_t{1};
            troughs_easting_tribs_current_trib_t = troughs_easting_all_tribs_current_basin_t{1};  % easting
            troughs_northing_tribs_current_trib_t = troughs_northing_all_tribs_current_basin_t{1};  % easting
            troughs_distance_all_tribs_current_basin_t = troughs_distance_all_tribs_current_basin_t{1};
            
            subplot(2,2,2,'align') % elev/chi plot for the entire basin - trunk stream only
            hold on
            for k = 1:length(sig_kps_magnitude_tribs_current_trib_t);  % plot each knickpoint seperately ( # of kps in that trib = length of sig_kps_magnitude_tribs_current_trib )
                plot(sig_kps_chi_tribs_current_trib_t(k),sig_kps_elevation_tribs_current_trib_t(k),'r.','MarkerSize',(sig_kps_magnitude_tribs_current_trib_t(k)/sig_kps_magnitude_all_tribs_all_kp_size_current_basin_max * KZP_parameters.max_knickpointsize2plot)) % size of knickpoint marker is relative to kp_magnitude value
                plot(troughs_chi_tribs_current_trib_t(k),troughs_elevation_tribs_current_trib_t(k),'bs','MarkerSize',(sig_kps_magnitude_tribs_current_trib_t(k)/sig_kps_magnitude_all_tribs_all_kp_size_current_basin_max * KZP_parameters.max_knickpointsize2plot)) % size of knickpoint marker is relative to kp_magnitude value
            end
            legend('Tributaries', 'Trunk', 'knickpoint top', 'knickpoint bottom', 'Location','northwest', 'FontSize',8,'FontWeight','bold')
            hold off
        end
        
        subplot(2,2,3,subplot3) % spatial plot on regional hillshade (last figure after all the basin individual figures are done)
        hold on
        for k = 1:length(sig_kps_magnitude_tribs_current_trib_t);% plot each knickpoint seperately ( # of kps in that trib = length of sig_kps_magnitude_tribs_current_trib )
            mapshow(sig_kps_easting_tribs_current_trib_t(k),sig_kps_northing_tribs_current_trib_t(k),'DisplayType', 'Point', 'Color', 'red', 'Marker', 'o', 'Markersize', sig_kps_magnitude_tribs_current_trib_t(k)/sig_kps_magnitude_all_tribs_all_kp_size_current_basin_max * KZP_parameters.max_knickpointsize2plot)
            h=mapshow(troughs_easting_tribs_current_trib_t(k),troughs_northing_tribs_current_trib_t(k),'DisplayType', 'Point', 'Color', 'blue', 'Marker', 's', 'Markersize', sig_kps_magnitude_tribs_current_trib_t(k)/sig_kps_magnitude_all_tribs_all_kp_size_current_basin_max * KZP_parameters.max_knickpointsize2plot);
            h.MarkerEdgeColor='blue';
        end
        hold off
        
        subplot(2,2,4,'align')
        hold on
        for k = 1: length(sig_kps_magnitude_tribs_current_trib_t);  % plot each knickpoint seperately ( # of kps in that trib = length of sig_kps_magnitude_tribs_current_trib )
            plot(sig_kps_distance_tribs_current_trib_t(k), sig_kps_elevation_tribs_current_trib_t(k), 'r.','MarkerSize',(sig_kps_magnitude_tribs_current_trib_t(k)/sig_kps_magnitude_all_tribs_all_kp_size_current_basin_max * KZP_parameters.max_knickpointsize2plot)); % size of knickpoint marker is relative to kp_magnitude value
            plot(troughs_distance_all_tribs_current_basin_t(k),troughs_elevation_tribs_current_trib_t(k),'bs','MarkerSize',(sig_kps_magnitude_tribs_current_trib_t(k)/sig_kps_magnitude_all_tribs_all_kp_size_current_basin_max * KZP_parameters.max_knickpointsize2plot)) % size of knickpoint marker is relative to kp_magnitude value
        end
        
        
        %add knickpoints - all tributaries
        for j = 1:length(sig_kps_chi_all_tribs_all_kp_size_current_basin);
            % unpack kp cell arrays for each tributary
            % (varies depending on the stream)
            sig_kps_chi_tribs_current_trib = sig_kps_chi_all_tribs_all_kp_size_current_basin{j};  % chi
            sig_kps_elevation_tribs_current_trib = sig_kps_elevation_all_tribs_all_kp_size_current_basin{j};  % elevation
            sig_kps_easting_tribs_current_trib = sig_kps_easting_all_tribs_all_kp_size_current_basin{j};  % easting
            sig_kps_northing_tribs_current_trib = sig_kps_northing_all_tribs_all_kp_size_current_basin{j};  % northing
            sig_kps_magnitude_tribs_current_trib =  sig_kps_magnitude_all_tribs_all_kp_size_current_basin{j}; %kp magnitude
            sig_kps_distance_tribs_current_trib =  sig_kps_distance_all_tribs_all_kp_size_current_basin{j}; %kp magnitude
            troughs_chi_tribs_current_trib = troughs_chi_all_tribs_current_basin{j};
            troughs_elevation_tribs_current_trib = troughs_elevation_all_tribs_current_basin{j};
            troughs_easting_tribs_current_trib = troughs_easting_all_tribs_current_basin{j};
            troughs_northing_tribs_current_trib = troughs_northing_all_tribs_current_basin{j};
            % Add Knickpoints to previous plots
            
            %calculate size of bubbles
            %        bubsizes = [floor(min(sig_kps_magnitude_tribs_current_trib)); round(quantile(sig_kps_magnitude_tribs_current_trib,[0.25, 0.5, 0.75],1)); ceil(max(sig_kps_magnitude_tribs_current_trib))];
            
            %         legentry{1} = 'Trunk stream';
            %         for l = 1:length(bubsizes)
            %             legentry{l+1} = num2str(bubsizes(l));
            %         end
            %         scatter(sig_kps_chi_tribs_current_trib,sig_kps_elevation_tribs_current_trib, sig_kps_magnitude_tribs_current_trib, 'filled', 'r','MarkerFaceColor','red') % size of knickpoint marker is relative to kp_magnitude value
            %         legend(legentry)
            
            
            subplot(2,2,3, subplot3) % spatial plot on regional hillshade (last figure after all the basin individual figures are done)
            hold on
            for k = 1:length(sig_kps_magnitude_tribs_current_trib);% plot each knickpoint seperately ( # of kps in that trib = length of sig_kps_magnitude_tribs_current_trib )
                kp_chi = round(sig_kps_chi_tribs_current_trib(k)); % get the nearest chi value for the pt (used to color the dot)
                if (sig_kps_magnitude_tribs_current_trib(k)/5) < 100
                    h=mapshow(sig_kps_easting_tribs_current_trib(k),sig_kps_northing_tribs_current_trib(k),'DisplayType', 'Point', 'Color', 'k', 'Marker', 'o','MarkerSize', (sig_kps_magnitude_tribs_current_trib(k)/kp_size2plot)); % same marker for kps as above
                    h.MarkerEdgeColor='red';
                    %plot(troughs_easting_tribs_current_trib(k),troughs_northing_tribs_current_trib(k),'k.','MarkerSize',8) % same marker for kps as above
                end
            end
            
            % add knickpoints to longitudinal profile
            subplot(2,2,4,'align')
            hold on
            for k = 1: length(sig_kps_magnitude_tribs_current_trib);  % plot each knickpoint seperately ( # of kps in that trib = length of sig_kps_magnitude_tribs_current_trib )
                plot(sig_kps_distance_tribs_current_trib(k), sig_kps_elevation_tribs_current_trib(k), 'r.','MarkerSize',(sig_kps_magnitude_tribs_current_trib(k)/kp_size2plot)); % size of knickpoint marker is relative to kp_magnitude value
            end
        end
        
        % Save to PDF
        knickpoints_plot_fname = sprintf('plots/%s_basin_%d_knickpoints.pdf', KZP_parameters.DEM_basename_nodir, i);
        if exist(knickpoints_plot_fname, 'file') ~= 2
            export_fig(knickpoints_plot_fname,KZP_parameters.quality_flag,'-pdf');
        end
        close all
    end
end
fprintf('\nfinished with KZP_knickzone_processing\n');
