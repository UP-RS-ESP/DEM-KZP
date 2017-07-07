function [sgolay_window_c_basin, sgolayfilt_order_c_basin, lumping_window_c_basin,...
    min_kp_size1_c_basin,min_kp_size2_c_basin,min_kp_slope_c_basin,Min_trib_size_c_basin,...
    Min_DA_threshold_c_basin,theta_ref_c_basin,ks_chiplot_c_basin,ksn_SA_c_basin,...
    theta_bf_c_basin,stream_id_c_basin,ks_chiplot_all,ksn_SA_all,theta_bf_all,...
    stream_id_all,sgolay_window_all,sgolayfilt_order_all,lumping_window_all,...
    min_kp_size1_all,min_kp_size2_all,min_kp_slope_all,Min_trib_size_all,...
    Min_DA_threshold_all,theta_ref_all] = ...
    E13_parameter_prep_csv(smoothing_window,sgolayfilt_order,...
    lumping_distance_upstream,min_kp_size1,min_kp_size2_magnitude,...
    min_kp_size2_relief,kp_relief_filter_option,kp_magnitude_filter_option,min_kp_slope,Min_trib_size,...
    Min_DA_threshold,theta_ref,ks_chiplot,ksn_SA_store,theta_bf,kz_up_da_c_basin_lips,ks_chiplot_all,...
    ksn_SA_all,theta_bf_all,stream_id_all,sgolay_window_all,sgolayfilt_order_all,...
    lumping_window_all,min_kp_size1_all,min_kp_size2_all,min_kp_slope_all,...
    Min_trib_size_all,Min_DA_threshold_all,theta_ref_all,stream_id_c_basin,...
    ksn_SA_c_basin,ks_chiplot_c_basin,theta_bf_c_basin,sgolay_window_c_basin,...
    sgolayfilt_order_c_basin,lumping_window_c_basin,min_kp_size1_c_basin,...
    min_kp_size2_c_basin,min_kp_slope_c_basin,Min_trib_size_c_basin,...
    Min_DA_threshold_c_basin)

% inputs: cell array containing knickzone attributes (to get number of
% knickzones in the basin)
% (current_num_kz, parameter values, basinwide metrics, vectors to store 
% parameter values for all knickzones)

% outputs: cell arrays containing parameter values for each basin {i}. list
% _all contains parameter values and basinwide metrics which will be stored
% for all knickzones

for i = 1:length(kz_up_da_c_basin_lips) % for each basin
    current_num_kz = length(kz_up_da_c_basin_lips{i}); % num kz in basin
% savitzky golay filter window size
    sgolay_window_basin_current = ones(current_num_kz,1)*smoothing_window;
    sgolay_window_c_basin{i} = sgolay_window_basin_current;
    
    % savitzky golay filter polynomial size
    sgolayfilt_order_basin_current = ones(current_num_kz,1)*sgolayfilt_order;
    sgolayfilt_order_c_basin{i} = sgolayfilt_order_basin_current;

    % chi lumping window
    lumping_window_basin_current = ones(current_num_kz,1)*lumping_distance_upstream;
    lumping_window_c_basin{i} = lumping_window_basin_current;
    
    % minimum knickzone magnitude 1
    min_kp_size1_basin_current = ones(current_num_kz,1)*min_kp_size1;
    min_kp_size1_c_basin{i} = min_kp_size1_basin_current;
    
    % minimum knickzone height 2 (depends on whether they chose to use
    % magntiude or releif to filter small knickzones
    if kp_magnitude_filter_option == 1
        min_kp_size2_basin_current = ones(current_num_kz,1)*min_kp_size2_magnitude;
        min_kp_size2_c_basin{i} = min_kp_size2_basin_current;
    end 
    
    if kp_relief_filter_option == 1
        min_kp_size2_basin_current = ones(current_num_kz,1)*min_kp_size2_relief;
        min_kp_size2_c_basin{i} = min_kp_size2_basin_current;
    end 
    
    % min kz slope
    min_kp_slope_basin_current = ones(current_num_kz,1)*min_kp_slope;
    min_kp_slope_c_basin{i} = min_kp_slope_basin_current;
    
    % minmum tributary size (cells)
    Min_trib_size_basin_current = ones(current_num_kz,1)*Min_trib_size;
    Min_trib_size_c_basin{i} = Min_trib_size_basin_current;
    
    % minmum drainage area (square meters)
    Min_DA_threshold_basin_current = ones(current_num_kz,1)*Min_DA_threshold;
    Min_DA_threshold_c_basin{i} = Min_DA_threshold_basin_current;
    
    % theta ref
    theta_ref_basin_current = ones(current_num_kz,1)*theta_ref;
    theta_ref_c_basin{i} = theta_ref_basin_current;
    
    % ks chiplot
    ks_chiplot_basin_current = ones(current_num_kz,1)*ks_chiplot(i);
    ks_chiplot_c_basin{i} = ks_chiplot_basin_current;
    
    % ksn chiplot
    ksn_SA_basin_current = ones(current_num_kz,1)*ksn_SA_store(i);
    ksn_SA_c_basin{i} = ksn_SA_basin_current;
    
    % theta chiplot
    theta_bf_basin_current = ones(current_num_kz,1)*theta_bf(i);
    theta_bf_c_basin{i} = theta_bf_basin_current;
    
    % stream id
    stream_id_basin_current = ones(current_num_kz,1)*(i);
    stream_id_c_basin{i} = stream_id_basin_current;
    
    % save basin specific statistics for all knickzones (total list)
    ks_chiplot_all = vertcat(ks_chiplot_all,ks_chiplot_c_basin{i});
    ksn_SA_all = vertcat(ksn_SA_all,ksn_SA_c_basin{i});
    theta_bf_all = vertcat(theta_bf_all,theta_bf_c_basin{i});
    stream_id_all = vertcat(stream_id_all,stream_id_c_basin{i});
    
    % save prameters used for all knickzones (total list) 
    sgolay_window_all = vertcat(sgolay_window_all,sgolay_window_c_basin{i});
    sgolayfilt_order_all = vertcat(sgolayfilt_order_all,sgolayfilt_order_c_basin{i});
    lumping_window_all = vertcat(lumping_window_all,lumping_window_c_basin{i});
    min_kp_size1_all = vertcat(min_kp_size1_all,min_kp_size1_c_basin{i});
    min_kp_size2_all = vertcat(min_kp_size2_all,min_kp_size2_c_basin{i});
    min_kp_slope_all = vertcat(min_kp_slope_all,min_kp_slope_c_basin{i});
    Min_trib_size_all = vertcat(Min_trib_size_all,Min_trib_size_c_basin{i});
    Min_DA_threshold_all = vertcat(Min_DA_threshold_all,Min_DA_threshold_c_basin{i});
    theta_ref_all = vertcat(theta_ref_all,theta_ref_c_basin{i});
end