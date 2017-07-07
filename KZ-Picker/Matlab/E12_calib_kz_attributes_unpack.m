function [all_kp_elevations, all_kp_chi,all_kp_easting,all_kp_northing,...
    all_kp_upstream_DA,all_kp_distance,all_kp_magnitude,all_kp_slope,...
    all_kp_length_distance,all_kp_relief,all_kp_elevations_bases,...
    all_kp_chi_bases,all_kp_easting_bases,all_kp_northing_bases,...
    all_kp_upstream_DA_bases,all_kp_distance_bases] = ...
    E12_calib_kz_attributes_unpack(kp_elev_stored_lips, kp_chi_stored_lips,...
    kp_easting_stored_lips, kp_northing_stored_lips,...
    kp_up_da_stored_lips, kp_distance_stored_lips, kp_magnitude_matrix_stored,...
    kp_face_slope_stored,kp_length_stored,kp_relief_stored,kp_elev_stored_bases,...
    kp_chi_stored_bases,kp_easting_stored_bases,kp_northing_stored_bases,...
    kp_up_da_stored_bases,kp_distance_stored_bases)

% messier version of reading stored outputs from knickzone selection. This script
% is only used in the calibration code because an older version of storing
% outputs is used in the calibration code

% ^^ inputs (cell arrays containing the knickzones found in each tributary
% {j} within each basin {i}

% code will organize these into a compiled list of knickzones with
% attributes (each row contains information for a knickzone lip or base and
% each column contains an attribute for that knickzone lip or base)

% need to first load in the current basin {i}, then vertically concatinate the all the
% information from each tributary {j} into one list, then vertically concatinate all of
% these from each basin {i} into one list



number_of_basins2 = length(kp_elev_stored_lips);
% get the number of basins that contained knickzones 

for i = 1:number_of_basins2 ;% % for all stream basins
    
        %for all the knickpoint lips
        current_basin_elevation = kp_elev_stored_lips{i};
        current_basin_chi = kp_chi_stored_lips{i};
        current_basin_easting = kp_easting_stored_lips{i};
        current_basin_northing = kp_northing_stored_lips{i};
        current_basin_upstream_DA = kp_up_da_stored_lips{i};
        current_basin_distance = kp_distance_stored_lips{i};
        % get the cell containing all the knickzones

        % measurements between lips and bases
        current_basin_magnitude = kp_magnitude_matrix_stored{i};
        current_basin_kp_slope = kp_face_slope_stored{i};
        current_basin_kp_length_distance = kp_length_stored{i};
        current_basin_kp_relief = kp_relief_stored{i}; 

        % for knickpoint bases
        current_basin_elevation_bases = kp_elev_stored_bases{i};
        current_basin_chi_bases = kp_chi_stored_bases{i};
        current_basin_easting_bases = kp_easting_stored_bases{i};
        current_basin_northing_bases = kp_northing_stored_bases{i};
        current_basin_upstream_DA_bases = kp_up_da_stored_bases{i};
        current_basin_distance_bases = kp_distance_stored_bases{i};


        % get the cell array containing all the knickpoints in all the
        % tributaries for the current basin (i)

        % create empty matricies where we will store the vertically
        % concatinated lists of knickpoint information for each basin

        %lips
        current_basin_all_kps_elevation = [];
        current_basin_all_kps_chi = [];
        current_basin_all_kps_easting = [];
        current_basin_all_kps_northing = [];
        current_basin_all_kps_upstream_DA = [];
        current_basin_all_kps_distance = [];

        %measurements
        current_basin_all_kps_magnitude = [];
        current_basin_all_kps_slope = [];
        current_basin_all_kps_length_distance = [];
        current_basin_all_kps_relief = [];

        %bases
        current_basin_all_elevation_bases = [];
        current_basin_all_chi_bases = [];
        current_basin_all_easting_bases = [];
        current_basin_all_northing_bases = [];
        current_basin_all_upstream_DA_bases = [];
        current_basin_all_distance_bases = [];


        for j = 1:length(current_basin_elevation);  
            % for each tributary, vertically concatinate all the knickpoints from each tributary within the basin
            
            %lips
            current_basin_all_kps_elevation = vertcat(current_basin_all_kps_elevation,current_basin_elevation{j});
            current_basin_all_kps_chi = vertcat(current_basin_all_kps_chi,current_basin_chi{j});
            current_basin_all_kps_easting = vertcat(current_basin_all_kps_easting,current_basin_easting{j});
            current_basin_all_kps_northing = vertcat(current_basin_all_kps_northing,current_basin_northing{j});
            current_basin_all_kps_upstream_DA = vertcat(current_basin_all_kps_upstream_DA,current_basin_upstream_DA{j});
            current_basin_all_kps_distance = vertcat(current_basin_all_kps_distance,current_basin_distance{j});
            %measurements
            current_basin_all_kps_magnitude = vertcat(current_basin_all_kps_magnitude,current_basin_magnitude{j});
            current_basin_all_kps_slope = vertcat(current_basin_all_kps_slope,current_basin_kp_slope{j});
            current_basin_all_kps_length_distance = vertcat(current_basin_all_kps_length_distance,current_basin_kp_length_distance{j});
            current_basin_all_kps_relief = vertcat(current_basin_all_kps_relief,current_basin_kp_relief{j});

            %bases
            current_basin_all_elevation_bases = vertcat(current_basin_all_elevation_bases,current_basin_elevation_bases{j});
            current_basin_all_chi_bases = vertcat(current_basin_all_chi_bases,current_basin_chi_bases{j});
            current_basin_all_easting_bases = vertcat(current_basin_all_easting_bases,current_basin_easting_bases{j});
            current_basin_all_northing_bases= vertcat(current_basin_all_northing_bases,current_basin_northing_bases{j});
            current_basin_all_upstream_DA_bases = vertcat(current_basin_all_upstream_DA_bases,current_basin_upstream_DA_bases{j});
            current_basin_all_distance_bases = vertcat(current_basin_all_distance_bases,current_basin_distance_bases{j});

        end
        % now all the knickpoints in the different tributaries have been
        % organized into 1 list.  but there are still seperate lists for
        % each seperate basin

        % lips
        all_basins_all_kps_elevation{i} = current_basin_all_kps_elevation;
        all_basins_all_kps_chi{i} = current_basin_all_kps_chi;
        all_basins_all_kps_easting{i} = current_basin_all_kps_easting;
        all_basins_all_kps_northing{i} = current_basin_all_kps_northing;
        all_basins_all_kps_upstream_DA{i} = current_basin_all_kps_upstream_DA;
        all_basins_all_kps_distance{i} = current_basin_all_kps_distance;

        %measurements
        all_basins_all_kps_magnitude{i} = current_basin_all_kps_magnitude;
        all_basins_all_kps_slope{i} = current_basin_all_kps_slope;
        all_basins_all_kps_length_distance{i} = current_basin_all_kps_length_distance;
        all_basins_all_kps_relief{i} = current_basin_all_kps_relief;

        %bases
        all_basins_all_kps_elevation_bases{i} = current_basin_all_elevation_bases;
        all_basins_all_kps_chi_bases{i} = current_basin_all_chi_bases;
        all_basins_all_kps_easting_bases{i} = current_basin_all_easting_bases;
        all_basins_all_kps_northing_bases{i} = current_basin_all_northing_bases;
        all_basins_all_kps_upstream_DA_bases{i} = current_basin_all_upstream_DA_bases;
        all_basins_all_kps_distance_bases{i} = current_basin_all_distance_bases;

end

% now lets combine the lists for each different basin, so we have a list of
% knickpoints for all basins
% create empty matrix to store total list

%lips
all_kp_elevations = [];
all_kp_chi = [];
all_kp_easting = [];
all_kp_northing = [];
all_kp_upstream_DA = [];
all_kp_distance = [];
% measurements
all_kp_slope = [];
all_kp_length_distance = [];
all_kp_magnitude = [];
all_kp_relief = [];

% bases
all_kp_elevations_bases = [];
all_kp_chi_bases = [];
all_kp_easting_bases = [];
all_kp_northing_bases = [];
all_kp_upstream_DA_bases = [];
all_kp_distance_bases = [];

% concatinate all the kp lists from each basin together into 1 list
for i = 1:number_of_basins2 % % for all stream basin
    % lips
    all_kp_elevations = vertcat(all_kp_elevations,all_basins_all_kps_elevation{i}); % 6
    all_kp_chi = vertcat(all_kp_chi,all_basins_all_kps_chi{i}); %5
    all_kp_easting = vertcat(all_kp_easting,all_basins_all_kps_easting{i}); %8
    all_kp_northing = vertcat(all_kp_northing,all_basins_all_kps_northing{i}); %9
    all_kp_upstream_DA = vertcat(all_kp_upstream_DA,all_basins_all_kps_upstream_DA{i}); %10
    all_kp_distance = vertcat(all_kp_distance,all_basins_all_kps_distance{i}); %11
    % measurements
    all_kp_magnitude = vertcat(all_kp_magnitude,all_basins_all_kps_magnitude{i}); %7
    all_kp_slope = vertcat(all_kp_slope,all_basins_all_kps_slope{i}); % 12
    all_kp_length_distance = vertcat(all_kp_length_distance,all_basins_all_kps_length_distance{i}); %13
    all_kp_relief = vertcat(all_kp_relief,all_basins_all_kps_relief{i}); 

    % bases
    all_kp_elevations_bases = vertcat(all_kp_elevations_bases,all_basins_all_kps_elevation_bases{i}); %6
    all_kp_chi_bases = vertcat(all_kp_chi_bases,all_basins_all_kps_chi_bases{i}); %5
    all_kp_easting_bases = vertcat(all_kp_easting_bases,all_basins_all_kps_easting_bases{i}); %8
    all_kp_northing_bases = vertcat(all_kp_northing_bases,all_basins_all_kps_northing_bases{i}); %9
    all_kp_upstream_DA_bases = vertcat(all_kp_upstream_DA_bases,all_basins_all_kps_upstream_DA_bases{i}); %10
    all_kp_distance_bases = vertcat(all_kp_distance_bases,all_basins_all_kps_distance_bases{i}); %11

end

% the resulting vectors contain a list of all the knickzone lips/bases elevations
% found.. all the chi coordinates for the knickzone lips/bases found.. ect

% each one of these attributes can be horizontally concatinated to make a
% table that contains all of the information for each knickzone. So if you
% wanted to plot the position of each knickzone in arcmap, you would use
% the 'all_kp_easting/northing' for the lips or bases to get the position
% information, and you could use the 'all_kp_magnitude' to plot the size of
% the knickzones for example.  

