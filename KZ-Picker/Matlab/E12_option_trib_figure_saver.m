%% Make figure for current tributary {j} Save figure in directory

mkdir trib_long_profile_figs;

handle_trib_fig = figure;set(handle_trib_fig,'units','inches','Position', [0 0 8 10], 'Visible', 'off')
% don't display figure

subplot(3,1,1)
plot(current_distance, current_elev,'k-')
hold on
plot(kp_distance,kp_elev,'r.','MarkerSize',16)
plot(kp_distance_bases,kp_elevation_bases,'b.','MarkerSize',14)
grid on
ylabel('elevation (m)')
xlabel('distance (m)')

subplot(3,1,2)
plot(current_chi, current_elev,'k-')
hold on
plot(kp_chi,kp_elev,'r.','MarkerSize',16)
plot(kp_chi_bases,kp_elevation_bases,'b.','MarkerSize',14)
grid on
ylabel('elevation (m)')
xlabel('chi (distance scaled by m/n and A)')

subplot(3,1,3)

hold on
plot(current_chi_resampled,smooth_detrended_elev_current,'color',[0 0 0.5], 'LineWidth',1.5)
plot(current_chi_resampled, detrended_elev_current,'k-')
plot(current_chi_resampled(sig_kps_cells_final),smooth_detrended_elev_current(sig_kps_cells_final),'r.','MarkerSize',16)
plot(current_chi_resampled(sig_kps_bases_final),smooth_detrended_elev_current(sig_kps_bases_final),'b.','MarkerSize',14)
grid on
ylabel('Detrended elevation (m)')
xlabel('chi (distance scaled by m/n and A)')

basin_num = num2str(i); % get current basin number
trib_num = num2str(j); % get current tributary number
long_prof_name = strcat('long_profile','_b_',basin_num,'t_',trib_num,'.png');
% name the figure according to basin and tributary number
saveas(handle_trib_fig,[pwd '/trib_long_profile_figs/',long_prof_name]);
