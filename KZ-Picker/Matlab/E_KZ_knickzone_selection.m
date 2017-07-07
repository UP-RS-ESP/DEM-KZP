%% Select Knickzones Using the Chi and Elevation Information from each Basin

mkdir maps_and_chi_plots;
mkdir csv_KZ_databases;

fprintf(1,'step 4 of 4: Selecting Knickzones' )
try
    number_of_basins = length(basin_index);
catch ME
    number_of_basins = 1;
end

for i = 1:number_of_basins
    % record the maximum elevation of the entire streamobj (needed for
    % creating histogram bins and easily plotting text on figures)
    current_AOI_chiplot = AOI_chiplot{i};
    max_elev_current(i) = max(current_AOI_chiplot.elev);
    min_elev_current(i) = min(current_AOI_chiplot.elev);
    min_east_current(i) = min(current_AOI_chiplot.x);
    max_north_current(i)= max(current_AOI_chiplot.y);
    max_dist_current(i)= max(current_AOI_chiplot.distance);
    
    % record the maximum chi of the entire streamobj (needed for
    % creating histogram bins and easily plotting text on figures)
    max_chi_current(i) = max(current_AOI_chiplot.chi);
    min_chi_current(i) = min(current_AOI_chiplot.chi);
end

% counters used to keep track of how many kncikzones are selected in each
% basin
kz_num_counter_prior = 0; % used to index the prior number of knickzones 
%(figure out how many were added after another iteration 'after analyzing
%another basin')
kz_num = 1;  % count knickzones (used as index to store information)

%% Create Vectors to store output knickzone attributes for each basin

% 1: knickpoint #
% 2: stream id #
% 3: stream ksn_SA
% 4: stream ks_chiplot
% 5: stream m/n (from chiplot)
% 6: tributary id #
% 7: tributary ks (from chiplot)
% 8: chi coordinate of knickzone lip/base
% 9: elevation coordinate of knickzone lip/base (m)
% 10: knickzone magnitude (detrended elevation drop (m))
% 11: knickzone relief (elevation drop across knickpoint (m))
% 12: easting coordinate of knickzone lip/base meters utm
% 13: northing coordinate of knickzone lip/base meters utm
% 14: upstream drainage area coordinate of knickzone lip/base m2
% 15: knickzone distance upstream coordinate of knickzone lip/base (m)
% 16: knickzone slope (elev drop/distance upstream)
% 17: knickzone length (distance usptream)
% 18: sgolay smoothing window size
% 19: sgolay polynomial order
% 20: knickpoint lumping search window size
% 21: minimum knickzone size pre-lumping
% 22: minimum knickzone size post-lumping (final minimum knickpoitn size)
% 23: minimum slope
% 24: minimum stream size for analysis (cells)
% 25: minimum drainage area for analysis (m2)
% 26: long profile for current tributary pathname ??? i don't know
% how to save the pathname to the longitudinal profile figure in
% the .csv file.  This is what we need to do so that users can
% click on icons in arcmap and look at a long-profile

% will go in database table as well as parameter values (which are
% constant) Parameters are compiled after the script selects knickzones



% We will store knickzone attributes in database 
% create vectors to store attributes (basinwide and ID numbers) % column #
kz_num_stored = []; % knickzone id number (for current basin)                       1
stream_id_all = []; % stream id number                                              2
ksn_SA_all = []; % ksn from slope area plot (uses ref concavity)                    3
ks_chiplot_all = []; % ks from chiplot (either fixed ref concavity or best-fit)     4
theta_bf_all = []; % theta used for chiplot construction (fixed or best fit)        5
trib_id_num_stored= []; % tributary id number (for current basin)                   6
ks_current_trib_stored = []; % steepness of tributary (for current basin)           7

% knickzone attributes
kz_chi_all_tribs_lips = []; % chi coordinate                                        8
kz_elev_all_tribs_lips = []; % elevation                                            9
kz_magnitude_all_tribs = []; % knickzone magnitude (detrended elevation change)     10
kz_relief_all_tribs = []; % relief of knickzone (m)                                 11
kz_easting_all_tribs_lips = []; % easting utm                                       12
kz_northing_all_tribs_lips = []; % northing utm                                     13
kz_up_da_all_tribs_lips = []; % drainage area upstream                              14
kz_distance_all_tribs_lips = []; % distance from river mouth                        15
kz_face_slope_all_tribs = []; % slope of knickzone (degrees)                        16
kz_length_all_tribs = []; % length of knickzone (m)                                 17


% bases
kz_chi_all_tribs_bases = []; %8
kz_elev_all_tribs_bases = []; %9
kz_easting_all_tribs_bases = [];  %12 % aslo records magnitude, relief, lenght and slope
kz_northing_all_tribs_bases = []; %13
kz_up_da_all_tribs_bases = [];  %14
kz_distance_all_tribs_bases = []; %15

% create vectors to stor parameters used
sgolay_window_all = []; % smoothing window size                                     18
sgolayfilt_order_all = []; % polynomial order                                       19
lumping_window_all = []; % lumping window size                                      20
min_kp_size1_all = []; % min kp size 1                                              21
min_kp_size2_all = []; % min kp size 2                                              22
min_kp_slope_all = []; % min kp slope                                               23
Min_trib_size_all = []; % min tributary size                                        24
Min_DA_threshold_all = []; % min drainage area for stream analysis                  25
theta_ref_all = []; % ref concavity for SA analysis (not saved)



% THESE ATRIBUTES ^^^ will also be saved for each basin
kz_up_da_c_basin_lips = {};
kz_distance_c_basin_lips = {};
kz_easting_c_basin_lips = {};
kz_northing_c_basin_lips = {};
kz_elev_c_basin_lips = {};
kz_chi_c_basin_lips = {};

% bases
kz_up_da_c_basin_bases = {};
kz_distance_c_basin_bases = {};
kz_easting_c_basin_bases = {};
kz_northing_c_basin_bases = {};
kz_elev_c_basin_bases = {};
kz_chi_c_basin_bases = {};

% geometry
kz_magnitude_c_basin = {};
kz_length_c_basin = {};
kz_face_slope_c_basin = {};
kz_relief_c_basin = {};

% tributary info and knickzone ids
kz_num_stored_basin = {};
ks_current_trib_stored_basin = {};
trib_id_num_stored_basin = {};

% paramter values and basinwide statistics
stream_id_c_basin = {};
ksn_SA_c_basin = {};
ks_chiplot_c_basin = {};
theta_bf_c_basin = {};
sgolay_window_c_basin = {};
sgolayfilt_order_c_basin = {};
lumping_window_c_basin = {};
min_kp_size1_c_basin = {};
min_kp_size2_c_basin = {};
min_kp_slope_c_basin = {};
Min_trib_size_c_basin = {};
Min_DA_threshold_c_basin = {};


% now begin knickzone selection
for i = 1:number_of_basins
        
    % First Step, Reorganize chiplot and streamobj array so each tributary can be analyzed seperately
            current_basin_chiplot = AOI_chiplot{i}; % chiplot structure array for calibration basin
            current_basin_streamobj = AOI_STR_largest{i};
            
            
            % use function to reorganize chiplot
            [elev_stored_trimmed, chi_stored_trimmed, eastingUTM11_stored_trimmed,...
            northingUTM11_stored_trimmed,upstream_DA_stored_trimmed,distance_stored_trimmed] = ...
            E1_sort_chiplot_to_tributaries(current_basin_chiplot,Min_trib_size);
            % within AOI_chiplot contains chiplots for all basins analyzed
            % within current_basin_chiplot, all of the chi values
            % of the stream network are listed.  Each stream segment (or tributary) is separated by a
            % NaN value which is the confluence between that tributary and the larger tributary/trunkstream.
            %
            % Cell 1 is the headwaters of the trunk
            % stream.  The first NaN is the mouth of the trunkstream.  The cell after
            % that is the headwaters of the first tributary.  The next NaN is the
            % confluence between that tributary and the trunk stream.  ECT for all tributaries....

            % we want different tributaries to each have their own
            % individual vector within a cell array so tributaries can be
            % analyzed separately
    
    % save the tributary vectors within the '_trimmed' cell arrays inside
    % of larger cell arrays for each basin {i}
    elev_stored_trimmed_basin{i} = elev_stored_trimmed ;
    chi_stored_trimmed_basin{i} = chi_stored_trimmed ;
    eastingUTM11_stored_trimmed_basin{i} = eastingUTM11_stored_trimmed ;
    northingUTM11_stored_trimmed_basin{i} = northingUTM11_stored_trimmed ;
    upstream_DA_stored_trimmed_basin{i} = upstream_DA_stored_trimmed ;
    distance_stored_trimmed_basin{i} = distance_stored_trimmed ;
    
    % now each basin (i), has an array of tributaries (j) that we can
    % pick knickzones in seperately.
end    
    
% Selecting Knickzones
fprintf('Selecting knickzones for all tributaries. At basin # of %d: ', ...
    number_of_basins);

%% For each basin, cycle through tributaries and identify knickzones

for i = 1:number_of_basins
    fprintf('%d, ', i);
    if mod(i,10) == 0
        fprintf('\n Selecting knickzones for all tributaries. At basin # of %d: ', ...
            number_of_basins);
    end
    
    % load cell array containing tributary vectors for the current basin
    chi_stored_current = chi_stored_trimmed_basin{i}; % chi info for each tributary within the stream {i}
    elev_stored_current  = elev_stored_trimmed_basin{i}; % elev info for each tributary within the stream {i}
    eastingUTM11_stored_current  = eastingUTM11_stored_trimmed_basin{i}; % ect..
    northingUTM11_stored_current  = northingUTM11_stored_trimmed_basin{i};
    upstream_DA_stored_current  = upstream_DA_stored_trimmed_basin{i};
    distance_stored_current  = distance_stored_trimmed_basin{i};
    
    % empty the cell arrays after each stream iteration or else data from
    % large catchments with lots of tributaries (j) are not fully overwritten by small streams
    % with less tributaries ({j} isn't as long, so less cells)

    % step below clears cell arrays that contain knickzone
    % information. needs to be done because one tributary may have
    % less knickzones than the next, and then each cell won't get
    % totally overwritten
    kp_magnitude_matrix_final={};
    kp_dist_final = {};
    kp_face_slope_final = {};
    kp_up_da= {};
    kp_distance= {};
    kp_easting = {};
    kp_northing = {};
    kp_chi = {};
    kp_elev = {};

    % knickpoint bases
    kp_up_da_bases = {};
    kp_distance_bases = {};
    kp_easting_bases = {};
    kp_northing_bases = {};
    kp_chi_bases = {};
    kp_elevation_bases = {};
    
    % Now the stream catchment is organized so each cell 
    % within array '_strored_trimmed' is a tributary we want to look at
    
    % Preform knickpoint picking on each tributary, cycle through all of the
    % tributaries and save the results.  This will be a loop that iterates for
    % the number of tributaries in the basin.  We will interpolate each
    % tributary, de-trend it, differentiate it, and select knickpoints
    
    % record number of tributaries of current catchment (number of cells
    % within the cell array containing all of the elevation.. or chi.. 
    % or distance.. or ect values for each tributary 
    number_of_tributaries = length(elev_stored_current);
    
    % loop through each tributary of the current catchment
    for j = 1:number_of_tributaries
        % get vector of chi values for the tributary
        current_chi = chi_stored_current{j};
        
        %verify if this chi vector is correct -- if it contains unusual chi
        %values, continue to next tributary
        if nanmax(current_chi) > 1e10
            continue
        end
        
        % get the vector of other important stored values for the current
        % tributary
        current_elev = elev_stored_current{j};
        current_easting = eastingUTM11_stored_current{j};
        current_northing = northingUTM11_stored_current{j};
        current_upDA = upstream_DA_stored_current{j};
        current_distance = distance_stored_current{j};
        
        % remove NaNs
        idx_chi = ~isnan(current_chi);
        current_chi_trib = current_chi(idx_chi) ; % get chi list  ( last value is NaN so don't include that)
        idx_elev = ~isnan(current_elev);
        current_elev_trib = current_elev(idx_elev); % get elev list ( last value is NaN so don't include that)
        idx_easting = ~isnan(current_easting);
        current_easting_trib = current_easting(idx_easting); % get easting list ( last value is NaN so don't include that)
        idx_northing = ~isnan(current_northing);
        current_northing_trib = current_northing(idx_northing); % get northing list ( last value is NaN so don't include that)
        idx_upDA = ~isnan(current_upDA);
        current_upDA_trib = current_upDA(idx_upDA); % get upstream drainage area list ( last value is NaN so don't include that)
        idx_distance = ~isnan(current_distance);
        current_distance_trib = current_distance(idx_distance); % get distance upstream list ( last value is NaN so don't include that)
                
        %% New function interpolates chi-elev data
            % interpolates the points in the chiplot because the chi 
            % values do not increase with constant intervals between 
            % successive values.  You shouldn�t use a Savitzky Golay 
            % filter on data that has a different stepsize between data
            % values, so the data must be interpolated over a constant 
            % measurement interval

            [current_chi_resampled, current_elev_resampled,current_stepsize_mean] = ...
            E2_chi_interpolate_stepsize(current_chi_trib, current_elev_trib);

        
        %%   Work with resampled data (use new function)
            % second step: calculate reference profile through headwater point
            % and outlet point (in chi elev space) and detrend the data
            [beta_current_trib, detrended_elev_current] = ...
            E3_reference_profile_construct_detrend(current_elev_resampled,...
            current_chi_resampled,current_stepsize_mean);  
        
        %% Smooth Data and Differentiate the detrended profile
            % use smoothing function            
            [smooth_detrended_elev_current,diff_smooth_detrended_elev_current] = ...
            E4_detrend_smooth_differentiate(current_chi_resampled, detrended_elev_current,...
            smoothing_window,sgolayfilt_order,smoothing_option);
            % smooths profile with sgolay filter and smoothing window
            % size

            % performs differentiation to find inflections in detrended
            % chi plot
        
       %% Now start picking out some knickpoints
            % new function searches the differentiated smoothed
            % detrended chi plot and finds inflections (from
            % oversteepened to understeepened and vice versa)
            [kps_cells_current_trib, base_cells_current_trib, kps_lips_length, kps_base_length,stream_mouth_cell] = ...
             E5_find_inflections_in_smooth_diff_detrended_chiplot(diff_smooth_detrended_elev_current);
        
        %% ORGANIZE THESE MATRICIES to make sure lips are compared to bases
            % reorganize matrices of inflection points so you know you are comparing
            % the upstream knickpoint lip to the downstream-adjacent knickpoint base

            % there's four possible conditions:

            % if there are more lips than bases, we need to compare the last lip to the
            % end of the stream.

            % if there are more bases than lips, we need to eliminate the upstream-most base.

            % or there could be the same number of lips and bases, but the first
            % inflection is a base, so we have to eliminate the upstream most base,
            % and add a base at the bottom of the stream.

            % Or everything could be fine, starting with a lip, ending with a
            % base.

            % we need to scan the stream from UPSTREAM to DOWNSTREAM and make
            % adjustments so we always have the last condition listed above

            % USE NEW FUNCTION 'potential_kps_inflections_sort'
            [base_cells_current_trib] = ...
             E6_potential_kps_inflections_sort(kps_cells_current_trib,kps_lips_length,...
             kps_base_length, base_cells_current_trib,stream_mouth_cell);
            %  All of the knickpoint cells are organized so the upstream 
            % knickzone lips are compared to a respective downstream knickpoint base
        
       %% Now We are ready to use the Knickpoint magnitude Filter (prelumping)
            % each knickpoint lip elevation
            % can be compared to the base directly downstream of it.
            % if the magnitude of detreneded elevation change doesn't exceed a
            % threshhold, we disregard it

            [sig_kps_cells, sig_kps_bases] = ...
             E7_min_kz_height_1_prelumping(kps_lips_length, smooth_detrended_elev_current,...
             kps_cells_current_trib,base_cells_current_trib,min_kp_size1);

            % done with the first knickpoint magnitude filter
        
%        Second Filter
        % LUMP KNICKPOINTS: looks for closely spaced
        % knickpoints and lumps them together into a bigger
        % knickpoint.

        % If knickpoints are within a certain
        % distance from one another
        % their magnitudes get summed and the downstream knickpoint
        % lip is erased. also the upstream knickpoint base is
        % erased.  this is to compare the elevation drop from the
        % upstream most knickpoint lip, and the downstream most
        % base.

        % this filter is a bit tough to follow, but it works

        [sig_kps_cells_lumped, sig_kps_bases_lumped] = ...
        E8_kz_lumping_filter(smooth_detrended_elev_current, sig_kps_bases,...
        sig_kps_cells,current_chi_resampled,current_chi,current_distance,...
        lumping_distance_upstream);
        
        
        %% 3rd filter KZ magnitude filter post lumping
                % calculate the magnitude of the knickpoints after lumping (knickpoints
                % that lumped together will have a larger magnitude now)
                % eliminate knickzones that are still too small
                
                % option to filter based on minimum knickzone magnitude or 
                % minimum knickzone relief
                
                [sig_kps_cells_lumped_filtered, sig_kps_bases_lumped_filtered,...
                 kp_magnitude_matrix_lumped_filtered] = ...
                E9_min_kz_height2_post_lumping(smooth_detrended_elev_current, sig_kps_cells_lumped,...
                sig_kps_bases_lumped,min_kp_size2_magnitude,min_kp_size2_relief,kp_magnitude_filter_option,...
                kp_relief_filter_option,current_elev_resampled);
            
        %% Final Filter (minimum steepness filter) OPTIONAL!

        %% Final Filter (minimum slope filter) OPTIONAL!

            % If the slope across the knickpoint doesn't exceed some value you can
            % get rid of that knickpoint here.  

            % tan-1(elev/distance) is the threshold slope in degrees here

            [sig_kps_cells_final, sig_kps_bases_final, kp_magnitude_matrix_final...
            ,kp_dist_final,kp_face_slope_final] = ...
            E10_minimum_kz_slope_filter(current_chi_resampled, sig_kps_cells_lumped_filtered,...
            current_chi,sig_kps_bases_lumped_filtered,current_distance,current_elev,...
            min_kp_slope,kp_magnitude_matrix_lumped_filtered);

            % All the filtering is done! We have our knickpoints identified!
            % Time to re-reference knickpoints back to orignal data
        
        %% Need to find cells of non-resampled data containing the knickzone positions and elevations
            % find the chi values for each knickpoint in the equispaced data, use that
            % to find the corresponding cell with the same chi value in the original
            % data.  then once you have all those cells that contain the knickpoints,
            % you can extract all the rest of the attributes for each knickponit
            % (upstream da, distance, northing easting ect..)

            % new function: kz_rereference_to_chiplot

            [knickpoint_relief, kp_up_da, kp_distance,kp_easting,kp_northing,...
            kp_chi,kp_elev,kp_up_da_bases,kp_distance_bases, kp_easting_bases,...
            kp_northing_bases,kp_chi_bases,kp_elevation_bases] = ...
            E11_kz_reference_to_chiplot(sig_kps_cells_final, sig_kps_bases_final,current_chi,...
            current_chi_resampled,current_elev,current_upDA,current_distance,...
            current_easting,current_northing);
        
    %% add option to plot and save a longitudinal profile for the current tributary
    
    if create_long_prof_4_tributaries == 1
        % function saves a long profile for each tributary in a new
        % directory. Labels these by the i and j # of the tributary
        % variable 'long_prof_name' contains the name of this label and can
        % be used to reference the long profile from the 'shapefile that
        % can be written for each tributary'.. I don't know how to write a
        % shapefile though
        E12_option_trib_figure_saver;
        
        current_trib_fig_path = ([pwd,'\calib_lips\',long_prof_name]);
        % this is the path to the figure that is written for each tributary
        % if we add this to a table i don't know if arcmap can open this
        % with the lightning bolt tool.  Also saving such a long filename
        % to a table?? does this work??
        
        store_trib_fig_path{j} = current_trib_fig_path;
        % store this for current tributary
        store_basin_trib_fig_path{i} = store_trib_fig_path;
        % store this for current basin
    end
    % would be used to display the longitudinal profile in arcmap if linked
    % to a shapefile.  
    
    %% ADD function below to write a shapefile for the
    % current tributary using the northing and easting vectors and the DEM
    % projection!  ^^^ use the PATH 'current_trib_fig_path' from the steps 
    % above to hyperlink the longitudinal profile figure to the shapefile
    % of the current tributary! 
    
    % I DON"T KNOW HOW TO DO THIS, but you'd want to use vectors: 'current_easting'
    % and 'current_northing' (these still contain NaN value at the outlet
    % of the tributary). If we don't to use the NaN values, use: 
    % 'current_easting_trib' and 'current_northing_trib'    
        
        
    
    
    %% store a list of ID numbers for each knickzone (also tributary channel steepness)
    
    % this is a bit complicated
    
        if ~isempty(kp_elevation_bases) 
            % if there are knickzones in the current tributary
            
            % create counter vectors to store number of knickzones in
            % current tributary
            kz_num_counter = []; 
            ks_current_trib = []; 
            trib_id_num = [];  % these are cleared after each tributary
            
            for k =1:length(kp_elevation_bases); 
                % for the number of knickzones found in the current tributary
                
                kz_num_counter(k,1) = kz_num ;
                % count the number of kncikzones in current tributary
                kz_num = kz_num +1;
                % increase the counter by 1
                ks_current_trib(k,1) = beta_current_trib;
                % record the steepness of the current tributary for each
                % knickzone occurance
                trib_id_num(k,1) = j;
                % record the tributary number for each knickzone occurance
            end
            
            kz_num_stored = vertcat(kz_num_stored,kz_num_counter);
            % add the number of knickzone occurances to the list of
            % knickzone occurances
            ks_current_trib_stored = vertcat(ks_current_trib_stored,ks_current_trib);
            % add the number of knickzone occurances to the list of
            % knickzone occurances
            trib_id_num_stored = vertcat(trib_id_num_stored,trib_id_num);
            % add the number of knickzone occurances to the list of
            % knickzone occurances
        end
        
        %% compile list of knickzone attributes (will be complete list of
        % all knickzones found)
        kz_up_da_all_tribs_lips = vertcat(kz_up_da_all_tribs_lips,kp_up_da);
        kz_distance_all_tribs_lips = vertcat(kz_distance_all_tribs_lips,kp_distance);
        kz_easting_all_tribs_lips = vertcat(kz_easting_all_tribs_lips,kp_easting);
        kz_northing_all_tribs_lips = vertcat(kz_northing_all_tribs_lips,kp_northing);
        kz_chi_all_tribs_lips = vertcat(kz_chi_all_tribs_lips,kp_chi);
        kz_elev_all_tribs_lips = vertcat(kz_elev_all_tribs_lips,kp_elev);
        
        % bases
        kz_up_da_all_tribs_bases = vertcat(kz_up_da_all_tribs_bases,kp_up_da_bases);
        kz_distance_all_tribs_bases = vertcat(kz_distance_all_tribs_bases,kp_distance_bases);
        kz_easting_all_tribs_bases = vertcat(kz_easting_all_tribs_bases,kp_easting_bases);
        kz_northing_all_tribs_bases = vertcat(kz_northing_all_tribs_bases,kp_northing_bases);
        kz_chi_all_tribs_bases = vertcat(kz_chi_all_tribs_bases,kp_chi_bases);
        kz_elev_all_tribs_bases = vertcat(kz_elev_all_tribs_bases,kp_elevation_bases);
        
        % geometry
        kz_magnitude_all_tribs = vertcat(kz_magnitude_all_tribs,kp_magnitude_matrix_final);
        kz_length_all_tribs = vertcat(kz_length_all_tribs,kp_dist_final);
        kz_face_slope_all_tribs = vertcat(kz_face_slope_all_tribs,kp_face_slope_final);
        kz_relief_all_tribs = vertcat(kz_relief_all_tribs,knickpoint_relief);
         
    end % end tributaries loop
    
    % save knickzone attributes for current basin only (use only the
    % knickzones added in the last iteration)!!! 
    current_kz_num = max(kz_num_stored);
    if kz_num_counter_prior < current_kz_num 
        % if more knickzones have been added from the current basin
        
        % add those knickzones to a new cell in a cell array so we can 
        % create knickzone stats for each basin individually
        
        % lips
        kz_up_da_c_basin_lips{i} = kz_up_da_all_tribs_lips((kz_num_counter_prior+1):current_kz_num);
        kz_distance_c_basin_lips{i} = kz_distance_all_tribs_lips((kz_num_counter_prior+1):current_kz_num);
        kz_easting_c_basin_lips{i} = kz_easting_all_tribs_lips((kz_num_counter_prior+1):current_kz_num);
        kz_northing_c_basin_lips{i} = kz_northing_all_tribs_lips((kz_num_counter_prior+1):current_kz_num);
        kz_elev_c_basin_lips{i} = kz_elev_all_tribs_lips((kz_num_counter_prior+1):current_kz_num);
        kz_chi_c_basin_lips{i} = kz_chi_all_tribs_lips((kz_num_counter_prior+1):current_kz_num);

        % bases
        kz_up_da_c_basin_bases{i} = kz_up_da_all_tribs_bases((kz_num_counter_prior+1):current_kz_num);
        kz_distance_c_basin_bases{i} = kz_distance_all_tribs_bases((kz_num_counter_prior+1):current_kz_num);
        kz_easting_c_basin_bases{i} = kz_easting_all_tribs_bases((kz_num_counter_prior+1):current_kz_num);
        kz_northing_c_basin_bases{i} = kz_northing_all_tribs_bases((kz_num_counter_prior+1):current_kz_num);
        kz_elev_c_basin_bases{i} = kz_elev_all_tribs_bases((kz_num_counter_prior+1):current_kz_num);
        kz_chi_c_basin_bases{i} = kz_chi_all_tribs_bases((kz_num_counter_prior+1):current_kz_num);

        % geometry
        kz_magnitude_c_basin{i} = kz_magnitude_all_tribs((kz_num_counter_prior+1):current_kz_num);
        kz_length_c_basin{i} = kz_length_all_tribs((kz_num_counter_prior+1):current_kz_num);
        kz_face_slope_c_basin{i} = kz_face_slope_all_tribs((kz_num_counter_prior+1):current_kz_num);
        kz_relief_c_basin{i} = kz_relief_all_tribs((kz_num_counter_prior+1):current_kz_num);

        % tributary info and knickzone ids
        kz_num_stored_basin{i} = kz_num_stored((kz_num_counter_prior+1):current_kz_num);
        ks_current_trib_stored_basin{i} = ks_current_trib_stored((kz_num_counter_prior+1):current_kz_num);
        trib_id_num_stored_basin{i} = trib_id_num_stored((kz_num_counter_prior+1):current_kz_num);
       
    else % add no new knickzones to the cell array for the current basin
        % lips
        kz_up_da_c_basin_lips{i} = [];
        kz_distance_c_basin_lips{i} = [];
        kz_easting_c_basin_lips{i} = [];
        kz_northing_c_basin_lips{i} = [];
        kz_elev_c_basin_lips{i} = [];
        kz_chi_c_basin_lips{i} = [];

        % bases
        kz_up_da_c_basin_bases{i} = [];
        kz_distance_c_basin_bases{i} = [];
        kz_easting_c_basin_bases{i} = [];
        kz_northing_c_basin_bases{i} = [];
        kz_elev_c_basin_bases{i} = [];
        kz_chi_c_basin_bases{i} = [];

        % geometry
        kz_magnitude_c_basin{i} = [];
        kz_length_c_basin{i} = [];
        kz_face_slope_c_basin{i} = [];
        kz_relief_c_basin{i} = [];
        
        % tributary info and knickzone ids
        kz_num_stored_basin{i} = [];
        ks_current_trib_stored_basin{i} = [];
        trib_id_num_stored_basin{i} = [];
        
    end
        
        % increase the 'prior' counter by the amount of knickzones found in the
        % past iteration
        kz_num_counter_prior = kz_num_counter_prior + (max(kz_num_counter)-kz_num_counter_prior);
    

end % end basin loop

%% ^^ ALL knickzones have been selected and stored


%% make vectors with parameter values used (for storing in database)

% for each basin and for all basins (these could be pre-allocated to make
% the script faster)


% function organizes parameters and basinwide metrics into vectors with
% same lenght as # kz in each basin and total # kz.  These are used to
% write the .csv file databases for knickzones selected
[sgolay_window_c_basin, sgolayfilt_order_c_basin, lumping_window_c_basin,...
min_kp_size1_c_basin,min_kp_size2_c_basin,min_kp_slope_c_basin,Min_trib_size_c_basin,...
Min_DA_threshold_c_basin,theta_ref_c_basin,ks_chiplot_c_basin,ksn_SA_c_basin,...
theta_bf_c_basin,stream_id_c_basin,ks_chiplot_all,ksn_SA_all,theta_bf_all,...
stream_id_all,sgolay_window_all,sgolayfilt_order_all,lumping_window_all,...
min_kp_size1_all,min_kp_size2_all,min_kp_slope_all,Min_trib_size_all,...
Min_DA_threshold_all,theta_ref_all] = ...
E13_parameter_prep_csv(smoothing_window,sgolayfilt_order,...
lumping_distance_upstream,min_kp_size1,min_kp_size2_magnitude,...
min_kp_size2_relief,kp_relief_filter_option,kp_magnitude_filter_option,min_kp_slope,Min_trib_size,...
Min_DA_threshold,theta_ref,ks_chiplot,ksn_SA_store,theta_bf,kz_up_da_c_basin_lips,ks_chiplot_all,...
ksn_SA_all,theta_bf_all,stream_id_all,sgolay_window_all,sgolayfilt_order_all,...
lumping_window_all,min_kp_size1_all,min_kp_size2_all,min_kp_slope_all,...
Min_trib_size_all,Min_DA_threshold_all,theta_ref_all,stream_id_c_basin,...
ksn_SA_c_basin,ks_chiplot_c_basin,theta_bf_c_basin,sgolay_window_c_basin,...
sgolayfilt_order_c_basin,lumping_window_c_basin,min_kp_size1_c_basin,...
min_kp_size2_c_basin,min_kp_slope_c_basin,Min_trib_size_c_basin,...
Min_DA_threshold_c_basin);

% ^ a bit bulky, but there are a lot of parameters we want to save

%% Now all knickzones and attributes are stored and ready for plotting/saving

% long profile and chiplot for each basin, also save database table for
% each basin

if ~isempty(kz_up_da_all_tribs_lips) 
% if knickzones were located
    for i = 1:number_of_basins
        % for each basin
        
        % get current basin ID number
        current_basin_id = basin_index(i);
        
        str_basin_id = num2str(current_basin_id); 
        % convert to text (from index from plot&select script (C_)
        
        % construct chi plot (plot each tributary seperately)
        
        current_basin_chi_values = chi_stored_trimmed_basin{i};
        current_basin_elev_values = elev_stored_trimmed_basin{i};
        
        % plot knickzone lips (weighted by size of knickzone, so need to
        % know which criteria user used to measure minimum knickzone size)
        
        % get knickzone attributes for current iteration needed for plots
        % lips
        current_kz_lips_elevation = kz_elev_c_basin_lips{i};
        current_kz_lips_chi = kz_chi_c_basin_lips{i};
        current_kz_lips_distance = kz_distance_c_basin_lips{i};
        %bases
        current_kz_bases_elevation = kz_elev_c_basin_bases{i};
        current_kz_bases_chi = kz_chi_c_basin_bases{i};
        current_kz_bases_distance = kz_distance_c_basin_bases{i};
        % geometry        
        current_kz_magnitude = kz_magnitude_c_basin{i};
        current_kz_relief = kz_relief_c_basin{i};
         
        H3 = figure;set(H3, 'Visible', 'off'); 
         
        for j = 1:length(current_basin_chi_values)
            plot(current_basin_chi_values{j},current_basin_elev_values{j},'k-')
            hold on
        end
        
         
        for k = 1:length(current_kz_lips_elevation)
            if kp_magnitude_filter_option ==1 % used knickzone magnitude
            plot(current_kz_lips_chi(k),current_kz_lips_elevation(k),...
                'k.','MarkerSize',(((current_kz_magnitude(k))^kp_plot_size)+10))
            end

            if kp_relief_filter_option == 1 % used knickzone relief
            plot(current_kz_lips_chi(k),current_kz_lips_elevation(k),...
                'k.','MarkerSize',(((current_kz_relief(k))^kp_plot_size)+10))
            end
        end

        % plot knickzone bases
        for k = 1:length(current_kz_bases_elevation)
            plot(current_kz_bases_chi(k),current_kz_bases_elevation(k),...
                'b.','MarkerSize',8)
        end
        
        theta_string = [' Theta Chiplot = ', num2str(theta_bf(i))];
        text(0,max(current_basin_elev_values{1}),theta_string);
        ksn_SA_string = [' SA ksn = ', num2str(ksn_SA_store(i))];
        text(0,max(current_basin_elev_values{1})-(max(current_basin_elev_values{1})/20),ksn_SA_string);
        ks_chiplot_current = [' Chiplot ks = ', num2str(ks_chiplot(i))];
        text(0,max(current_basin_elev_values{1})-(max(current_basin_elev_values{1})/10),ks_chiplot_current);
        theta_ref_string = [' Ref Concavity = ', num2str(theta_ref)];
        text(0,max(current_basin_elev_values{1})-(max(current_basin_elev_values{1})/7),theta_ref_string);


        xlabel('Chi')
        ylabel('Elev (m)')
        title('Chi Plot with Knickzones Scaled to Knickzone Height')
        grid on
        
        % change basin_id so saves as different filename for each basin
        chiplot_fig_name = ['Chi_plot_w_KZ_basin_', str_basin_id];
        % need to add path to directory where we store figures
        saveas(H3,[pwd '/maps_and_chi_plots/',chiplot_fig_name,'.png']);
        % save figure in directory
        clf
        
        %% construct longitudinal profile (plot each tributary seperately)
        current_basin_dist_values = distance_stored_trimmed_basin{i};
        current_basin_elev_values = elev_stored_trimmed_basin{i};

        H4 = figure;set(H4, 'Visible', 'off') ;
             
        for j = 1:length(current_basin_dist_values)
            plot(current_basin_dist_values{j},current_basin_elev_values{j},'k-')
            hold on
        end
        % plot knickzone lips
        for k = 1:length(current_kz_lips_distance)
            if kp_magnitude_filter_option ==1
            plot(current_kz_lips_distance(k),current_kz_lips_elevation(k),...
                'k.','MarkerSize',(((current_kz_magnitude(k))^kp_plot_size)+10))
            end

            if kp_relief_filter_option ==1
            plot(current_kz_lips_distance(k),current_kz_lips_elevation(k),...
                'k.','MarkerSize',(((current_kz_relief(k))^kp_plot_size)+10))
            end
        end

        % plot knickzone bases
        for k = 1:length(current_kz_lips_distance)
            plot(current_kz_bases_distance(k),current_kz_bases_elevation(k),...
                'b.','MarkerSize',8)
        end

        smoothing_window_str = [' Smoothing Window = ', num2str(smoothing_window),' cells'];
        text(0,max(current_basin_elev_values{1}),smoothing_window_str);
        lumping_window_str = [' lumping Window = ', num2str(lumping_distance_upstream),' m'];
        text(0,max(current_basin_elev_values{1})-(max(current_basin_elev_values{1})/20),lumping_window_str);
        min_kz1_str = ['Min Kz 1 = ', num2str(min_kp_size1),' m'];
        text(0,max(current_basin_elev_values{1})-(max(current_basin_elev_values{1})/10),min_kz1_str );
        
        if kp_magnitude_filter_option == 1
            min_kz2_str  = [' Min Kz 2 = ', num2str(min_kp_size2_magnitude)];
            text(0,max(current_basin_elev_values{1})-(max(current_basin_elev_values{1})/7),min_kz2_str );
        end
        
        if kp_relief_filter_option == 1;
            min_kz2_str  = [' Min Kz 2 = ', num2str(min_kp_size2_relief),' m'];
            text(0,max(current_basin_elev_values{1})-(max(current_basin_elev_values{1})/7),min_kz2_str );
        end
            
        xlabel('Distance (m)')
        ylabel('Elev (m)')
        title('Long. Profile with Knickzones Scaled to Knickzone Height')
        grid on
        % Export figure (should save in a different directory)
        %Generate name for figure
        long_prof_filename = ['Long_profile_KZ_basin_', str_basin_id];
        % need to add path to directory where we store figures
        saveas(H4,[pwd '/maps_and_chi_plots/',long_prof_filename,'.png']);
        % save figure in directory
        clf
        
        
        %% save database files
        % into a .csv file.  
        
        % list of database measurments and parameters (separate table for lips and
        % bases)
        % 1: knickpoint #
        % 2: stream id #
        % 3: stream ksn_SA
        % 4: stream ks_chiplot
        % 5: stream m/n (from chiplot)
        % 6: tributary id #
        % 7: tributary ks (from chiplot)
        % 8: chi coordinate of knickzone lip/base
        % 9: elevation coordinate of knickzone lip/base (m)
        % 10: knickzone magnitude (detrended elevation drop (m))
        % 11: knickzone relief (elevation drop across knickpoint (m))
        % 12: easting coordinate of knickzone lip/base meters utm
        % 13: northing coordinate of knickzone lip/base meters utm
        % 14: upstream drainage area coordinate of knickzone lip/base m2
        % 15: knickzone distance upstream coordinate of knickzone lip/base (m)
        % 16: knickzone slope (elev drop/distance upstream)
        % 17: knickzone length (distance usptream)
        % 18: sgolay smoothing window size
        % 19: sgolay polynomial order
        % 20: knickpoint lumping search window size
        % 21: minimum knickzone size pre-lumping
        % 22: minimum knickzone size post-lumping (final minimum knickpoitn size)
        % 23: minimum slope
        % 24: minimum stream size for analysis (cells)
        % 25: minimum drainage area for analysis (m2)
        % 26: long profile for current tributary pathname ??? i don't know
        % how to save the pathname to the longitudinal profile figure in
        % the .csv file.  This is what we need to do so that users can
        % click on icons in arcmap and look at a long-profile
        
        Knickpoint_Database_Lips_CB = horzcat (kz_num_stored_basin{i},...
            stream_id_c_basin{i},ksn_SA_c_basin{i},ks_chiplot_c_basin{i},...
            theta_bf_c_basin{i},trib_id_num_stored_basin{i},...
            ks_current_trib_stored_basin{i},kz_chi_c_basin_lips{i},...
            kz_elev_c_basin_lips{i},kz_magnitude_c_basin{i},...
            kz_relief_c_basin{i},kz_easting_c_basin_lips{i},kz_northing_c_basin_lips{i},...
            kz_up_da_c_basin_lips{i},kz_distance_c_basin_lips{i},kz_face_slope_c_basin{i},...
            kz_length_c_basin{i},sgolay_window_c_basin{i},sgolayfilt_order_c_basin{i},...
            lumping_window_c_basin{i},min_kp_size1_c_basin{i},min_kp_size2_c_basin{i},...
            min_kp_slope_c_basin{i},Min_trib_size_c_basin{i},Min_DA_threshold_c_basin{i});
        
        
        Knickpoint_Database_Bases_CB = horzcat (kz_num_stored_basin{i},...
            stream_id_c_basin{i},ksn_SA_c_basin{i},ks_chiplot_c_basin{i},...
            theta_bf_c_basin{i},trib_id_num_stored_basin{i},...
            ks_current_trib_stored_basin{i},kz_chi_c_basin_bases{i},...
            kz_elev_c_basin_bases{i},kz_magnitude_c_basin{i},...
            kz_relief_c_basin{i},kz_easting_c_basin_bases{i},kz_northing_c_basin_bases{i},...
            kz_up_da_c_basin_bases{i},kz_distance_c_basin_bases{i},kz_face_slope_c_basin{i},...
            kz_length_c_basin{i},sgolay_window_c_basin{i},sgolayfilt_order_c_basin{i},...
            lumping_window_c_basin{i},min_kp_size1_c_basin{i},min_kp_size2_c_basin{i},...
            min_kp_slope_c_basin{i},Min_trib_size_c_basin{i},Min_DA_threshold_c_basin{i});
         
        % write .csv file for each basin
        
        oldfolder = pwd; % save the working directory
        
        % create a header for the .csv files

        hdr = {'1kp_id', '2stream_id', '3ksn_SA', '4ks_chi', '5m_n', '6trib_id', ...
                        '7trib_ks', '8chi', '9Elev', '10kz_mag', '11kz_rel', '12kz_E', ...
                        '13kz_N', '14kz_DA', '15kz_dFm', '16kz_slp', '17kz_leng', ...
                        '18sg_wz', '19sg_po', '20lwsz', '21mkz_1', ...
                        '22mkz_2','23mkz_slp','24m_trib','25min_DA'};
                    txt = sprintf('%s, ',hdr{:}); % change from cell array to text
                    txt(end)='';
                    
        % convert basin id # to string so you can save in filename
        str_basin_id = num2str(basin_index(i));          
        % basin_index is the vector inputted by the user to specify which 
        % basin numbers they wanted to analyze. From script:
        % 'C_KZ_drainage_basin_plot_and_select'
              
        % specify filename for knickzone lips
        kp_lips_fn = strcat('kp_lips_basin', str_basin_id,'.csv');             
        %give filename
        cd([oldfolder,'/csv_KZ_databases/'])
        fid = fopen(kp_lips_fn, 'w');
        fprintf(fid, txt);
        fclose(fid);

        % write .csv file containing knickzone lip information for
        % current parameters (need one row offset for the header)
        dlmwrite([oldfolder,'/csv_KZ_databases/',kp_lips_fn],Knickpoint_Database_Lips_CB,'-append','delimiter',',', 'precision', '%.5f','roffset',1);
        
        % specify filename for knickzone bases
        kp_bases_fn = strcat('kp_bases_basin', str_basin_id,'.csv');             
        %give filename
        fid = fopen(kp_bases_fn, 'w');
        fprintf(fid, txt);
        fclose(fid);
        % ^^ add legend
        
        % write .csv file containing knickzone lip information for
        % current parameters
        dlmwrite([oldfolder,'/csv_KZ_databases/',kp_bases_fn],Knickpoint_Database_Bases_CB,'-append','delimiter',',', 'precision', '%.5f','roffset',1);
        
        cd(oldfolder); % return to working directory
    end
end


% Write .csv files for all of the knickzones (so these can be plotted in
% arcmap all at once)

Knickpoint_Database_Lips = horzcat (kz_num_stored,stream_id_all,...
        ksn_SA_all, ks_chiplot_all, theta_bf_all,trib_id_num_stored,...
        ks_current_trib_stored,kz_chi_all_tribs_lips,kz_elev_all_tribs_lips,...
        kz_magnitude_all_tribs,kz_relief_all_tribs,kz_easting_all_tribs_lips,...
        kz_northing_all_tribs_lips,kz_up_da_all_tribs_lips,kz_distance_all_tribs_lips,...
        kz_face_slope_all_tribs,kz_length_all_tribs,sgolay_window_all,...
        sgolayfilt_order_all,lumping_window_all,min_kp_size1_all,...
        min_kp_size2_all,min_kp_slope_all,Min_trib_size_all, Min_DA_threshold_all);
     
Knickpoint_Database_Bases = horzcat (kz_num_stored,stream_id_all,...
        ksn_SA_all, ks_chiplot_all, theta_bf_all,trib_id_num_stored,...
        ks_current_trib_stored,kz_chi_all_tribs_bases,kz_elev_all_tribs_bases,...
        kz_magnitude_all_tribs,kz_relief_all_tribs,kz_easting_all_tribs_bases,...
        kz_northing_all_tribs_bases,kz_up_da_all_tribs_bases,kz_distance_all_tribs_bases,...
        kz_face_slope_all_tribs,kz_length_all_tribs,sgolay_window_all,...
        sgolayfilt_order_all,lumping_window_all,min_kp_size1_all,...
        min_kp_size2_all,min_kp_slope_all,Min_trib_size_all, Min_DA_threshold_all);
 
% write the table

% specify filename for knickzone lips
kp_lips_fn = strcat('kz_lips_all','.csv');             
%give filename
cd([oldfolder,'/csv_KZ_databases/'])
fid = fopen(kp_lips_fn, 'w');
fprintf(fid, txt);
fclose(fid);

% write .csv file containing knickzone lip information for
% current parameters (need one row offset for the header)
dlmwrite([oldfolder,'/csv_KZ_databases/',kp_lips_fn],Knickpoint_Database_Lips,'-append','delimiter',',', 'precision', '%.5f','roffset',1);

% specify filename for knickzone bases
kp_bases_fn = strcat('kp_bases_all','.csv');             
%give filename
fid = fopen(kp_bases_fn, 'w');
fprintf(fid, txt);
fclose(fid);
% ^^ add legend

% write .csv file containing knickzone lip information for
% current parameters
dlmwrite([oldfolder,'/csv_KZ_databases/',kp_bases_fn],Knickpoint_Database_Bases,'-append','delimiter',',', 'precision', '%.5f','roffset',1);

cd(oldfolder); 
    
