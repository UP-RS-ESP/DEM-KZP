function [algorithm_kp_easting, algorithm_kp_northing,...
    algorithm_kp_relief,algorithm_kp_sgolay_smooth_current,...
    algorithm_kp_chi_lump_current,algorithm_kp_minkp_1_pre_current,...
    algorithm_kp_minkp_2_post_current,num_alg_kps_current] = ...
    F3_calib_load_algorithm_output_table(Knickpoint_Database_Lips_all_params,i)

% ^^ inputs: File structure array (contains filename), iteration #
% outputs: easting, northing, relief, and parameter values of knickzone
% lips selected by the algorithm for the current iteration (each iteration
% is a different parameter result)


% load that file
current_algorithm_result_data=Knickpoint_Database_Lips_all_params{i};

% record important metrics from tables

algorithm_kp_easting = current_algorithm_result_data(1:end,12);
algorithm_kp_northing = current_algorithm_result_data(1:end,13);
% use these  to find TP, FP, and FN

algorithm_kp_relief = current_algorithm_result_data(1:end,11);
% use this to Calculate 'G'

% parameters used in prior simulation
algorithm_kp_sgolay_smooth_current = current_algorithm_result_data(2,18);
algorithm_kp_chi_lump_current = current_algorithm_result_data(2,20);
algorithm_kp_minkp_1_pre_current = current_algorithm_result_data(2,21);
algorithm_kp_minkp_2_post_current = current_algorithm_result_data(2,22);
% record the parameters used for that simulation

num_alg_kps_current = length(algorithm_kp_easting);
% record the # of algorithm knickpoints for current param combination    