function [sig_kps_cells_lumped_filtered, sig_kps_bases_lumped_filtered,...
    kp_magnitude_matrix_lumped_filtered] = ...
    E9_calib_min_kz_height2_post_lumping(smooth_detrended_elev_current, sig_kps_cells_lumped,...
    sig_kps_bases_lumped,min_kp_size2_calib,b,kp_magnitude_filter_option,...
    kp_relief_filter_option,current_elev_resampled)
% calculate the magnitude of the knickpoints after lumping (knickpoints
% that lumped together will have a larger magnitude now)
% eliminate knickzones that are still too small

% can filter based on minimum knickzone magnitude or minimum knickzone
% relief

kp_magnitude_matrix_lumped_equispaced = smooth_detrended_elev_current(sig_kps_cells_lumped) - smooth_detrended_elev_current(sig_kps_bases_lumped);

if kp_magnitude_filter_option == 1
 for k = 1:length(sig_kps_cells_lumped); 
     % for how many potential knickpoints we still have

    if kp_magnitude_matrix_lumped_equispaced(k) < min_kp_size2_calib(b)
% specify a 2nd minimum knickpoint size if you'd like (after summing closely spaced knickpoints together)

        kp_magnitude_matrix_lumped_equispaced(k) = NaN;  
% remove the stored magnitude for knickpoints that are too small  
        sig_kps_bases_lumped(k) =NaN; 
% remove the stored base position for knickpoints that are too small 
        sig_kps_cells_lumped(k) = NaN; 
% remove the stored lips for knickpoints that are too small 
    end
 end

% remove knickpoints that aren't big enough even after lumping close
% proximity knickpoints
idx_kps_tribs_3= ~isnan(sig_kps_cells_lumped);
sig_kps_cells_lumped_filtered = sig_kps_cells_lumped(idx_kps_tribs_3);   
% index out knickpoints that after lumping still didn't get big enough

idx_bases_tribs_3= ~isnan(sig_kps_bases_lumped);
sig_kps_bases_lumped_filtered = sig_kps_bases_lumped(idx_bases_tribs_3);  
% index out respective bases to those knickpoints

idx_magnitude_matrix_2 = ~isnan(kp_magnitude_matrix_lumped_equispaced);
kp_magnitude_matrix_lumped_filtered  = kp_magnitude_matrix_lumped_equispaced(idx_magnitude_matrix_2);  % magnitudes of those knickponts that didn't get big enough
end

% if you want to filter by knickzone relief rather than magnitude
if kp_relief_filter_option == 1
% remember the knickzone relief increases upstream
% systematically in a concave up stream profile (stream
% gradient increases upstream)
current_rel = current_elev_resampled(sig_kps_cells_lumped) - current_elev_resampled(sig_kps_bases_lumped);

for k = 1:length(sig_kps_cells_lumped); 
    % for how many potential knickpoints we still have
    if current_rel(k) < min_kp_size2_calib(b) 
% specify a 2nd minimum knickpoint size if you'd like (after summing closely spaced knickpoints together)

        kp_magnitude_matrix_lumped_equispaced(k) = NaN;  
% remove the stored magnitude for knickpoints that are too small  
        sig_kps_bases_lumped(k) =NaN; 
% remove the stored base position for knickpoints that are too small 
        sig_kps_cells_lumped(k) = NaN; 
% remove the stored lips for knickpoints that are too small 
    end
 end

% remove knickpoints that aren't big enough even after lumping close
% proximity knickpoints
idx_kps_tribs_3= ~isnan(sig_kps_cells_lumped);
sig_kps_cells_lumped_filtered = sig_kps_cells_lumped(idx_kps_tribs_3);   % index out knickpoints that after lumping still didn't get big enough

idx_bases_tribs_3= ~isnan(sig_kps_bases_lumped);
sig_kps_bases_lumped_filtered = sig_kps_bases_lumped(idx_bases_tribs_3);  % index out respective bases to those knickpoints

idx_magnitude_matrix_2 = ~isnan(kp_magnitude_matrix_lumped_equispaced);
kp_magnitude_matrix_lumped_filtered  = kp_magnitude_matrix_lumped_equispaced(idx_magnitude_matrix_2);  
% relief of those knickponts that didn't get big enough
end