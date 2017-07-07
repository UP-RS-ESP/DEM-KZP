function [kp_id_num,sgolay_window_size, sgolay_polynomial, max_kp_search_dist,...
    min_kp_size_pre_combination,min_kp_size_post_combination,...
    min_kp_slope_d,min_trib_size_cells,min_DA] = ...
    KZP_calib_save_parameters(all_kp_elevations, SG_smoothing_calib,l,...
    lumping_search_distance_calib,c,min_kp_size1_calib,a,...
    min_kp_size2_calib,b,min_kp_slope,Min_trib_size,sgolayfilt_order,...
    Min_DA_threshold)

% ^^ inputs (all_kp_elevations is used to total the total number of
% knickzones found.  The rest of the inputs are parameter values.

% outputs: vector that assigns each knickzone (row) the parameter value
% used to make the knickzone selection (column). Since these are constants
% for the current calibration iteration, the parameter value will be the
% same throughout each row.


%
total_number_of_kps = length(all_kp_elevations);
% Total  number of knickzones

kp_id_number = 1:1:total_number_of_kps;
% create vector to label each knickzone boundary by # (1 -> total #
% knickzones)
kp_id_num = transpose(kp_id_number); 

% parameters used  
sgolay_window_size = zeros(total_number_of_kps,1);
sgolay_polynomial = zeros(total_number_of_kps,1);
max_kp_search_dist = zeros(total_number_of_kps,1);
min_kp_size_pre_combination = zeros(total_number_of_kps,1);
min_kp_size_post_combination = zeros(total_number_of_kps,1);
min_kp_slope_d = zeros(total_number_of_kps,1);
min_trib_size_cells = zeros(total_number_of_kps,1);
min_DA = zeros(total_number_of_kps,1);

for k = 1 : total_number_of_kps 
% make a vector as long as the # of knickzones found recording the parameter values used for the current knickzone selections
    sgolay_window_size(k) = SG_smoothing_calib(l);     
    max_kp_search_dist(k) = lumping_search_distance_calib(c) ;  
    min_kp_size_pre_combination(k)= min_kp_size1_calib(a) ;   
    min_kp_size_post_combination(k)= min_kp_size2_calib(b);     
    % ^^ change during calibration script
    
    min_kp_slope_d (k) = min_kp_slope ;  
    min_trib_size_cells(k) = Min_trib_size; 
    sgolay_polynomial(k)= sgolayfilt_order;
    min_DA(k) = Min_DA_threshold;
    % ^^ do not change during calibration script (could recode to optimize
    % these parameters if this is desired)
end

