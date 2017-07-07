function [TP, FP, FN,Confirmed_kp_pair] = ...
    F4_calib_alg_tp_or_fp_and_fn(error_radius, num_alg_kps_current,...
    numb_calib_kz_ref,calib_northing_str_loc,...
    algorithm_kp_northing,calib_easting_str_loc,algorithm_kp_easting)

% inputs: allowable spatial error, number of algorithm knickzones (lips or
% bases), number of calibration knickzones (lips or bases), northing and
% easting locations of calibration knickzone lips/bases and algorithm
% lips/bases

% outputs: number of true-positives, false positives, false negatives, and
% the cells of knickzones that produced true positives and the residual
% distance between these knickzone boundaries.

% create 2 vectors to store the where knickpoints are TP and FP and FN
    Rel_alg = zeros(num_alg_kps_current,1);
    Rel_calib = zeros(numb_calib_kz_ref,1);
    
    % this matrix stores exactly which knickpoint pairs created a TP and
    % the distance between these pairs
    Confirmed_kp_pair = [];
    
    n =1; % counter for number of True Positive knickpoints 
    
    for  j = 1:num_alg_kps_current % loop through calibration points
        for k = 1:numb_calib_kz_ref  % loop through algorithm points
            D_lips =(((calib_northing_str_loc(k) - algorithm_kp_northing(j))^2)+...
            (calib_easting_str_loc(k) - algorithm_kp_easting(j))^2)^0.5;
             if error_radius > D_lips  
            % if the error radius is larger than the distance between the
            % current algorithm (j) and current calibration (k) knickpoint,
            % record a TP
            Rel_alg(j) = 1;
            Rel_calib(k) = 1;
            
            Confirmed_kp_pair(n,:) = [j k D_lips]; 
            % save the alg. kp #, calib kp # and the distance btwn the two
            n = n+1;
            % increase the TP counter
             end
             % if the error radius is smaller than the distance between the
            % current algorithm (j) and current calibration (k) knickpoint,
            % a FP will be recorded becase the value of Rel_alg will still
            % be 0
        end
    end
    
    idx_Rel_alg_TP = find(Rel_alg ==1);
    % record # of TP (true positives), algorithm kpt, with calib kpt nearby
    TP = length(idx_Rel_alg_TP);
    
    % record # of FP (false positives), algorithm kpt with no calib kpt
    % (total # knickzone lips - # tp
    FP = num_alg_kps_current - TP;
    
    % find all the calibration knickzones that were missed (still have
    % value of 0
    idx_calib_FN =  find(Rel_calib ==0);
    FN = length(idx_calib_FN);
    % caluclate number of false negatives (missed calibration knickpoint)
    