function [sig_kps_cells, sig_kps_bases] = ...
    E7_min_kz_height_1_prelumping(kps_lips_length, smooth_detrended_elev_current,...
    kps_cells_current_trib,base_cells_current_trib,min_kp_size1)

n2 = kps_lips_length;  % for all the potential knickpoints found
for k =1:1:n2;
    if smooth_detrended_elev_current(kps_cells_current_trib (k))-...
            smooth_detrended_elev_current (base_cells_current_trib(k)) < min_kp_size1 ;
        % see uf the elevation difference between the upstream knickpoint 
        %lip and downstream base exceeds a minimum knickpoint size

        kps_cells_current_trib(k) = NaN;  % sets the knickpoints that are too small to NaN
        base_cells_current_trib(k) = NaN;  % does the same thing with the bases corresponding to those knickpoints
    end
end

idx_kps_tribs= ~isnan(kps_cells_current_trib); % index out those inflections that were too small
sig_kps_cells= kps_cells_current_trib(idx_kps_tribs);  %cells of potential knickpoints

idx_base= ~isnan(base_cells_current_trib);  % index out those inflections that were too small
sig_kps_bases = base_cells_current_trib(idx_base);  %cells of potential respective bases