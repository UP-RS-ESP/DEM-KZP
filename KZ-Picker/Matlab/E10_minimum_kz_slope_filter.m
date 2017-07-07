function [sig_kps_cells_final, sig_kps_bases_final, kp_magnitude_matrix_final...
    ,kp_dist_final,kp_face_slope_final] = ...
    E10_minimum_kz_slope_filter(current_chi_resampled, sig_kps_cells_lumped_filtered,...
    current_chi,sig_kps_bases_lumped_filtered,current_distance,current_elev,...
    min_kp_slope,kp_magnitude_matrix_lumped_filtered)

% create empty cell to store the upstream distance values for the potential
% knickzone locations (need to rereference the pre-interpolated positions
% to calculate the slope of the knickzone in d(elevation)/d(distance)
kp_cells_dist = [];
kp_bases_dist = [];

kp_cells_chi_equispaced_slope = current_chi_resampled(sig_kps_cells_lumped_filtered); 
% chi coordinates for knickpoint lips
% reference the original chi dataset, find where the values are closest to
% the knickpoint chi coordinates above

for k = 1:length(kp_cells_chi_equispaced_slope)
    [~,idx] = min(abs(current_chi-kp_cells_chi_equispaced_slope(k))) ;  % find cell with the minimum difference (usually pretty small)
    kp_cells_dist(k) = idx; % these are the cells in the original data which contain the knickpoints
end

kp_bases_chi_equispaced_slope = current_chi_resampled(sig_kps_bases_lumped_filtered); 
% chi coordinates for knickpoint bases

for k = 1:length(kp_bases_chi_equispaced_slope)
    [~,idx] = min(abs(current_chi-kp_bases_chi_equispaced_slope(k))) ;  
    % find cell with the minimum difference (usually pretty small)
    kp_bases_dist(k) = idx; 
    % these are the cells in the original data which contain the
    % knickzone boundaries located (remember we had to interpolate
    % to use a constant chi-stepsie before) so the knickzone
    % boundary containing cells found are different from the cells
    % containing the original chi-elev position of the knickzone
    % boundaries
end

% index the knickzone distance upstream coordinates from the
% original vector that stored the distance upstream for each
% coordinate in the tributary

kp_length_dist = current_distance(kp_cells_dist)- current_distance(kp_bases_dist);  
% calculate length of knickpoint
kp_height = current_elev(kp_cells_dist) - current_elev(kp_bases_dist);

kp_face_slope_d = atand(kp_height./kp_length_dist) ; % divide drop in elevation by change in distance downstream

% This is an optional filter that would filter out
% knickpoints with face slopes that are too gradual

for k = 1:length(kp_length_dist) % for the number of potential knickpoitns we still have
    if kp_face_slope_d(k) < min_kp_slope;% Parmeter line 105

        sig_kps_cells_lumped_filtered(k) = NaN; % mark cells that belong to too gradual knickpoint
        sig_kps_bases_lumped_filtered(k)=NaN;  % mark cells that belong to too gradual knickpoints
        kp_magnitude_matrix_lumped_filtered (k) = NaN;  % mark cells that belong to too gradual knickpoint
        kp_length_dist(k) = NaN;  % mark cells that belong to too gradual knickpoint
        kp_face_slope_d(k) = NaN;  % mark cells that belong to too gradual knickpoint

    end
end

% remove knickpoint cells that were too gradual
idx_kps_tribs_4= ~isnan(sig_kps_cells_lumped_filtered);
sig_kps_cells_final = sig_kps_cells_lumped_filtered(idx_kps_tribs_4);  % index out knickpoint lips to knickpoints that were to gradual

idx_bases_tribs_4= ~isnan(sig_kps_bases_lumped_filtered);
sig_kps_bases_final = sig_kps_bases_lumped_filtered(idx_bases_tribs_4); % respective bases

idx_magnitude_matrix_3 = ~isnan(kp_magnitude_matrix_lumped_filtered);
kp_magnitude_matrix_final  = kp_magnitude_matrix_lumped_filtered(idx_magnitude_matrix_3);  % respective knickpoint magnitudes

idx_kp_dist_length = ~isnan(kp_length_dist);
kp_dist_final = kp_length_dist(idx_kp_dist_length );                       % respective knickpoint lengths (in units of chi)

idx_kp_face_slope_dist_space_detrended  = ~isnan(kp_face_slope_d);
kp_face_slope_final= kp_face_slope_d(idx_kp_face_slope_dist_space_detrended);   % respective knickpoint slopes  detrended elev/chi
