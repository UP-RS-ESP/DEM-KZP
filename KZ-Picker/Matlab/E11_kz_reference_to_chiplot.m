function [knickpoint_relief, kp_up_da, kp_distance,kp_easting,kp_northing,...
    kp_chi,kp_elev,kp_up_da_bases,kp_distance_bases, kp_easting_bases,...
    kp_northing_bases,kp_chi_bases,kp_elevation_bases] = ...
    E11_kz_reference_to_chiplot(sig_kps_cells_final, sig_kps_bases_final,current_chi,...
    current_chi_resampled,current_elev,current_upDA,current_distance,...
    current_easting,current_northing)

% create empty vector to store cells that contain the knickzone lips and
% bases
kp_cells = [];
kp_cells_bases=[];

kp_cells_chi_equispaced = current_chi_resampled(sig_kps_cells_final); 
% chi coordinates for knickpoint lips (resampled space)
kp_cells_bases_chi_equispaced = current_chi_resampled(sig_kps_bases_final); 
% chi coordinates for knickpoint bases (resampled space)

% reference the original chi dataset, find where the values are closest to
% the knickpoint chi coordinates above
for k = 1:length(kp_cells_chi_equispaced)
    [~,idx] = min(abs(current_chi-kp_cells_chi_equispaced(k))) ;  % find cell with the minimum difference (usually pretty small)
    kp_cells(k) = idx; % these are the cells in the original data which contain the knickpoints
end


% Do the same thing for the bases
for k = 1:length(kp_cells_bases_chi_equispaced)
    [~,idx_bases] = min(abs(current_chi-kp_cells_bases_chi_equispaced(k)));
    kp_cells_bases(k) = idx_bases; % these are the cells in the original data which contain the bases
end

% "kp_cells" and "kp_cells_bases" are the list of cells that contain the
% lips and bases of all the knickpoints.  with these cell numbers, we can
% get the northing/easting position of the knickpoints and the upstream
% drainage area, distance upstream, ect...

% Get the other metrics we want (upstream DA, distance, ect..
% Index out the attributes for the knickpoint cells and knickpoint base cells


% knickpoint relief
current_kp_lip_elev = current_elev(kp_cells);
current_kp_bases_elev = current_elev(kp_cells_bases);
knickpoint_relief = current_kp_lip_elev - current_kp_bases_elev;


% knickpoint lips
kp_up_da = current_upDA(kp_cells);  % drainage area upstream of knickpoint
kp_distance = current_distance(kp_cells); % knickpoint distance upstream
kp_easting = current_easting(kp_cells);
kp_northing = current_northing(kp_cells);
kp_chi = current_chi(kp_cells); % chi coordinate
kp_elev = current_elev(kp_cells); % elevation of knickpoint


% knickpoint bases
kp_up_da_bases = current_upDA(kp_cells_bases);  % drainage area upstream of knickpoint
kp_distance_bases = current_distance(kp_cells_bases); % knickpoint distance upstream
kp_easting_bases = current_easting(kp_cells_bases);
kp_northing_bases = current_northing(kp_cells_bases);
kp_chi_bases = current_chi(kp_cells_bases); % chi coordinate
kp_elevation_bases = current_elev(kp_cells_bases); % elevation of knickpoint
