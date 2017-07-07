function [G] = F5_calib_geometric_accuracy(Confirmed_kp_pair, num_alg_kps_current,...
algorithm_kp_relief,Relief_calib)

% inputs: index number of algorithm knickzone and calibration knickzones
% that retrieved true positives, and the distance between all of those
% knickzones (within 'confirmed_kp_pair').  Number of algorithm knickzones
% in the current iteration (with the parameter values used)

% outputs: geometric accuracy 1-0

Alg_kp_num = Confirmed_kp_pair(:,1);
% alg kpt numb that was identified as a TP
Calib_kp_num = Confirmed_kp_pair(:,2);
% calib kpt numb that was identified as a TP
Dist = Confirmed_kp_pair(:,3); % distance recorded btwn knickpoints

knickpoint_match = []; % clear this matrix

n = 1;
for j = 1:num_alg_kps_current
    % loop through each algorithm knickpoint

    conf_alg_kpt_idx = find(Alg_kp_num == j); 
    % find if alg kpt confirmed multiple calib kpts

    if ~isempty(conf_alg_kpt_idx)
        % if there was a true positive (confirmed alg. kpt)

        current_kp = Confirmed_kp_pair(conf_alg_kpt_idx,:);
        % find all the instances of TP (maybe more than 1 if alg 
        %knickpoint sits next to multiple calibration knickpoints)

        [min_value,best_match_idx] = min(current_kp(:,3)); 
        % get minimum distance out of these (because there may be
        % more than 1)

        idx_best_match_kp = find(Dist ==min_value);
        % find index of calib kpt that created this minimum distance

        knickpoint_match(n)=idx_best_match_kp(1);
        % store this index, will record the algorithm kp # and
        % corresponding calibration kp # (these are a knickpoint
        % pair)

        n = n+1; 
        % increase the counter to store the next index
    end

end

TP_Alg_kp_for_comp = Alg_kp_num(knickpoint_match);
TP_calib_kp_for_comp = Calib_kp_num(knickpoint_match);
% Use the minimum distances to find the best TP algorithm and
% calibration knickpoint pairs (the ones that are closest matches)

relief_TP_alg = algorithm_kp_relief(TP_Alg_kp_for_comp);
% relief of TP algorithm knickpoints

relief_TP_calib = Relief_calib(TP_calib_kp_for_comp);
% relief of corresponding TP calibration knickpoints

% calculate the geometric error for all knickpoints (error between
% measured knickpoint relief values)
Geometric_error = 1-abs(relief_TP_alg - relief_TP_calib)./(relief_TP_calib+relief_TP_alg);

% take the mean of this (score of 1 is perfect)
G = mean(Geometric_error);
