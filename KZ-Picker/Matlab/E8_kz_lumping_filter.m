function [sig_kps_cells_lumped, sig_kps_bases_lumped] = ...
    E8_kz_lumping_filter(smooth_detrended_elev_current, sig_kps_bases,...
    sig_kps_cells,current_chi_resampled,current_chi,current_distance,...
    lumping_distance_upstream)

% If knickpoints are within a certain
% distance from one another
% their magnitudes get summed and the downstream knickpoint
% lip is erased. also the upstream knickpoint base is
% erased.  this is to compare the elevation drop from the
% upstream most knickpoint lip, and the downstream most
% base.

% this filter is a bit tough to follow, but it works

% get the detrended elevation of the potential lips and bases
kp_det_elev_bases = smooth_detrended_elev_current(sig_kps_bases); % OPTIONAL
kp_det_elev_lips = smooth_detrended_elev_current(sig_kps_cells); 

% need to find the distance upstream for the potential knickzone boundaries
% because we resampled the chi data, this is a bit complicated
kp_cells_d = []; % store cells containing potential kz lips
kp_cells_db = []; % store cells containing potential kz bases

kp_cells_chi_equispaced_filtering = current_chi_resampled(sig_kps_cells); 
% chi coordinates for knickpoint lips

% reference the original chi dataset, find where the values are closest to
% the knickpoint chi coordinates above (this is to find
% cells in the un-resampled data that contain the potential
% knickzone boundaries)
for k = 1:length(kp_cells_chi_equispaced_filtering)
    [~,idx] = min(abs(current_chi-kp_cells_chi_equispaced_filtering(k))) ;  
    % find cell with the minimum difference (usually pretty small)
    kp_cells_d(k) = idx; 
    % these are the cells in the original data which contain the knickpoints
end

% do same for bases
kp_cells_chi_equispaced_filtering_b = current_chi_resampled(sig_kps_bases); % chi coordinates for knickpoint lips
% reference the original chi dataset, find where the values are closest to
% the knickpoint chi coordinates above
for k = 1:length(kp_cells_chi_equispaced_filtering_b)
    [~,idx] = min(abs(current_chi-kp_cells_chi_equispaced_filtering_b(k))) ;  
    % find cell with the minimum difference (usually pretty small)
    kp_cells_db(k) = idx; 
    % these are the cells in the original data which contain the knickpoints
end

sig_kps_distance = current_distance(kp_cells_d); 
% get the distance upstream coordinate for each potential knickzone lip
sig_kps_distance_b = current_distance(kp_cells_db); 
% get the distance upstream coordinate for each potential knickzone base


% create empty matrix to store lumped knickzones found
% (these are kz lips or bases that will be deleted as the
% lumped kz lips or bases are combined together)
lumpcell_base = zeros(length(sig_kps_distance),1);
lumpcell_lips = zeros(length(sig_kps_distance),1);

for k = 1:length(sig_kps_distance)-1 

    if  sig_kps_distance_b(k) - sig_kps_distance(k+1) < lumping_distance_upstream  
        % if the knickzone lip of the downstream knickzone and base of the upstream knickzone
        % are close together than the lumpingwindow (distance upstream)
        % combine both knickzones and sum their magnitudes

        lumpcell_base(k) = NaN;  % set the downstream kp to NaN 
        lumpcell_lips(k+1) = NaN  ;   % set the base for the upstream kp to NaN 
    end 
    %( we want to compare the upstream knickpoint elev to
      %the downstream base elev to get the total elev drop)
end


% remove the knickpoints that got lumped (nan's)
idx_kps_tribs_2 = ~isnan(lumpcell_lips);
sig_kps_cells_lumped = sig_kps_cells(idx_kps_tribs_2);   % remove the knickpoint lip locations of knickpoints that got lumped
idx_bases_tribs_2= ~isnan(lumpcell_base);
sig_kps_bases_lumped  = sig_kps_bases(idx_bases_tribs_2);  % remove respective bases of knickpoints that got lumped
