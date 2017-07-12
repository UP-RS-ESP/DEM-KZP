
if KZP_parameters.show_figs == 1 || KZP_parameters.show_figs == 2
    fprintf('\nMaking figures for all drainage basins.\nAt basin # of %d: ', number_of_basins);
    for i = 1:number_of_basins
        fprintf('%d, ', i);
        knickpoints_plot_fname = sprintf('%s%s%s_knickpoints_basin_%02d.pdf', KZP_parameters.KZP_plots_dirname, KZP_parameters.dir_sep, KZP_parameters.DEM_basename_nodir, i);
        if exist(knickpoints_plot_fname, 'file') == 2
            continue
        end
        if ~isempty(AOI_STR_S_slopearea_dbasins{i})
            % turn the stream name into a text that can be pasted to the chiplot
            current_stream_id = sprintf('Stream # %d',i);
            str = sprintf('%s', current_stream_id);

            idx_basin = AOI_dbasins_unique(i);
            AOI_dbasins_current = AOI_dbasins;
            idx_ncurrent_basin = find(AOI_dbasins_current.Z ~= idx_basin);
            AOI_dbasins_current.Z(idx_ncurrent_basin) = 0;
            AOI_DEM_c = crop(AOI_DEM, AOI_dbasins_current, 0);
            idx=find(AOI_DEM_c.Z == 0);
            AOI_DEM_c.Z(idx) = NaN;
            
            % this is a chi-plot of all of the basin streams to be analyzed. we will add
            % kps later
            elev_current_stream = elev_stored_trimmed_basin{i} ; % get the elevation data for current basin from (elev_stored_trimmed_basin) <- contains all basins
            chi_current_stream = chi_stored_trimmed_basin{i} ;% get the chi data for current basin from (elev_stored_trimmed_basin) <- contains all basins

            %prepare variables for knickpoint-plotting
            current_kz_lips_elevation = kz_elev_all_tribs_lips{i};
            current_kz_lips_chi = kz_chi_all_tribs_lips{i};
            current_kz_lips_distance = kz_distance_all_tribs{i};
            current_kz_lips_easting = kz_easting_all_tribs_lips{i};
            current_kz_lips_northing = kz_northing_all_tribs_lips{i};

            %bases
            current_kz_bases_elevation = kz_elev_all_tribs_bases{i};
            current_kz_bases_chi = kz_chi_all_tribs_bases{i};
            current_kz_bases_distance = kz_distance_all_tribs_bases{i};
            current_kz_bases_easting = kz_easting_all_tribs_bases{i};
            current_kz_bases_northing = kz_northing_all_tribs_bases{i};

            % geometry
            current_kz_magnitude = kz_magnitude_all_tribs{i};
            current_kz_relief = kz_relief_all_tribs{i};

            figure
            set(gcf,'units','normalized','position',[0 0 1 1]);
            set(gcf, 'PaperOrientation', 'landscape');
            set(gcf, 'PaperType', KZP_parameters.PaperType_size);
            set(gcf,'Visible', 'off');
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
            fitresult_logspace.a = AOI_dbasins_stats_min_area(i,25);
            fitresult_logspace.b = AOI_dbasins_stats_min_area(i,24);
            p_logspace_a = AOI_dbasins_stats_min_area(i,28);
            p_logspace_b = AOI_dbasins_stats_min_area(i,29);
            s = sprintf('From slope-area analyses (binned) (DA=%3.2f km^2): k_{sn} trunk = %3.1f with theta = %0.2f, k_{s} trunk = %3.1f with theta = %0.2f\nk_{s} trunk = %3.1f with theta from chi plot = %0.2f, k_{s} trunk = %3.1f +/- %3.1f with theta from robust regression = %0.2f +/- %0.2f\nr^2 = %0.2f, rmse=%2.2e, p(k_{s}) = %1.2e, p(theta) = %1.2e', ...
                max(AOI_STR_S_slopearea_dbasins{i}.a)/1e6, AOI_STR_S_slopearea_dbasins_trunk{i}.ks, AOI_STR_S_slopearea_dbasins_trunk{i}.theta, ...
                AOI_STR_S_slopearea_dbasins_trunk_adj{i}.ks, AOI_STR_S_slopearea_dbasins_trunk_adj{i}.theta, ...
                AOI_STR_S_slopearea_dbasins_trunk_chitheta{i}.ks, AOI_STR_S_slopearea_dbasins_trunk_chitheta{i}.theta, ...
                fitresult_logspace.a, abs(fitresult_logspace.a-ci_logspace(1)), fitresult_logspace.b, abs(fitresult_logspace.b-ci_logspace(3)), ...
                gof_logspace{i}.adjrsquare, gof_logspace{i}.rmse, p_logspace_a, p_logspace_b);
            foo=axis;
            text(0.01, 0.9, s, 'Fontsize', 8, 'Units', 'normalized')
            grid on
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


            % LRP
            subplot(2,2,4,'align')
            h = plotdz(AOI_STR_streams_dbasins_unique{i},AOI_DEM);
            set(h, 'color', [0.5 0.5 0.5], 'Linewidth', 1);
            hold on
            grid on
            h2 = plotdz(AOI_STR_all_streams_trunk{i},AOI_DEM);
            set(h2, 'color', 'k', 'Linewidth', 2);
            hold off
            
            if length(current_kz_lips_elevation) == 0
                %no knickpoint found
                fprintf('i = %d, no knickpoint found\n', i);
            end

        % THE FOLLOWING scales knickpoint magnitude in plot by all knickpoint sizes for the entire DEM.        
        %         %determine max. knickpointsize for plotting
        %         if ~isempty(kz_magnitude_all_tribs)
        %             foo = (cellfun(@(x) prctile(x(:),95), kz_magnitude_all_tribs,'un',0));
        %             kz_magnitude_c_basin_max = ceil(max(cell2mat([foo(:)]))); clear foo
        %             if kz_magnitude_c_basin_max > KZP_parameters.max_knickpointsize2plot
        %                 %scale knickpoint for plotting
        %                 kp_size2plot = kz_magnitude_c_basin_max/KZP_parameters.max_knickpointsize2plot;
        %             elseif kz_magnitude_c_basin_max < KZP_parameters.max_knickpointsize2plot
        %                 kp_size2plot = kz_magnitude_c_basin_max/KZP_parameters.max_knickpointsize2plot;
        %             else
        %                 kp_size2plot = 1;
        %             end
        %         else
        %             kp_size2plot = 1;
        %         end
        %         %increase size by factor 10
        %         kp_size2plot = kp_size2plot * 0.1;

        %This scales knickpoint plot size for the current catchment only
            %determine max. knickpointsize for plotting
            if ~isempty(current_kz_magnitude)
                kz_magnitude_c_basin_max = ceil(max(current_kz_magnitude));
                if kz_magnitude_c_basin_max > KZP_parameters.max_knickpointsize2plot
                    %scale knickpoint for plotting
                    kp_size2plot = kz_magnitude_c_basin_max/KZP_parameters.max_knickpointsize2plot;
                elseif kz_magnitude_c_basin_max < KZP_parameters.max_knickpointsize2plot
                    kp_size2plot = kz_magnitude_c_basin_max/KZP_parameters.max_knickpointsize2plot;
                else
                    kp_size2plot = 1;
                end
            else
                kp_size2plot = 1;
            end
            %increase size by factor 10
            %kp_size2plot = kp_size2plot * 0.1;

            %add knickpoints - first trunk
            subplot(2,2,2,'align') % elev/chi plot for the entire basin - trunk stream only
            hold on
            for k = 1:length(current_kz_magnitude)
                % plot each knickpoint seperately ( # of kps in that trib = length of sig_kps_magnitude_tribs_current_trib )
                plot(current_kz_lips_chi(k),current_kz_lips_elevation(k),'r.','MarkerSize',(current_kz_magnitude(k)/kp_size2plot)) % size of knickpoint marker is relative to kp_magnitude value
                plot(current_kz_bases_chi(k),current_kz_bases_elevation(k),'bs','MarkerSize',(current_kz_magnitude(k)/kp_size2plot)) % size of knickpoint marker is relative to kp_magnitude value
            end
            hold off

            subplot(2,2,3,'align') % spatial plot on regional hillshade (last figure after all the basin individual figures are done)
            hold on
            %Hillshade with adjusted K_sn and knickpoints
            symbolspec_ksn045 = makesymbolspec('line',...
                {'ksn045' [prctile([AOI_STR_MS.ksn045], 5) prctile([AOI_STR_MS.ksn045], 95)] 'color' jet(6)});
            symbolspec_ks_adj = makesymbolspec('line',...
                {'ks_adj' [prctile([AOI_STR_MS.ks_adj], 5) prctile([AOI_STR_MS.ks_adj], 95)] 'color' jet(6)});
            subplot(2,2,3, 'align')
            imageschs(AOI_DEM, AOI_DEM, 'caxis', ...
                [nanmin((AOI_DEM_c.Z(:))) ceil(nanmax(AOI_DEM_c.Z(:)))],...
                'colormap',gray,'colorbar',false);
            %imageschs(AOI_DEM, AOI_DEM, 'caxis', [floor(min(AOI_DEM(:))) ceil(max(AOI_DEM(:)))], 'colormap',gray,'colorbar',false)
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
            for k = 1:length(current_kz_lips_easting) % plot each knickpoint seperately ( # of kps in that trib = length of sig_kps_magnitude_tribs_current_trib )
                plot(current_kz_lips_easting(k),current_kz_lips_northing(k),'r.','MarkerSize',(current_kz_magnitude(k)/kp_size2plot)) % size of knickpoint marker is relative to kp_magnitude value
                plot(current_kz_bases_easting(k),current_kz_bases_northing(k),'bs','MarkerSize',(current_kz_magnitude(k)/kp_size2plot)) % size of knickpoint marker is relative to kp_magnitude value
            end
            hold off

            subplot(2,2,4,'align')
            hold on
            for k = 1: length(current_kz_magnitude)  % plot each knickpoint seperately ( # of kps in that trib = length of sig_kps_magnitude_tribs_current_trib )
                plot(current_kz_lips_distance(k),current_kz_lips_elevation(k),'r.','MarkerSize',(current_kz_magnitude(k)/kp_size2plot)) % size of knickpoint marker is relative to kp_magnitude value
                plot(current_kz_bases_distance(k),current_kz_bases_elevation(k),'bs','MarkerSize',(current_kz_magnitude(k)/kp_size2plot)) % size of knickpoint marker is relative to kp_magnitude value
            end


            % Save to PDF
            if exist(knickpoints_plot_fname, 'file') ~= 2
                if exist('export_fig') == 2
                    export_fig(knickpoints_plot_fname,KZP_parameters.quality_flag,'-pdf');
                else
                    saveas(gcf,knickpoints_plot_fname, 'pdf')
                end
            end
            close all
        end
    end
end

