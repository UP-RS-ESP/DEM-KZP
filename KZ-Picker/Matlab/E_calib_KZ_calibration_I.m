%% Knickpoint calibration preparation (1 basin, Basinwide)

% code by Al Neely 10/4/2015 and Bodo Bookhagen

% this code takes ~1hr to run and could be optimized.... a lot. 

% GOAL:

% gather knickpoint positions and attributes while looping through various
% parameter combinations.  

% This code analyzes a single drainage basin and selects knickpoints through each iteration
% using a different parameter combination.  

% The knickpoints selected with the particular parameter combination are
% saved and archived after each iteration in a .csv file.

% later the code will reopen this archive and analyze which
% parameter combination produces the lowest misfit with the calibration
% knickpoints spread across the entire catchment (calibration_II)

% keep i = 1 because we're only working with 1 drainage basin

% loop through different parameter combinations and store output CSV files
% that contain the dimensions and positions for each knickpoint for that
% parameter combination

% these .csv files will be read by the next script to identify their
% accuracy with the calibration knickpoints

mkdir calib_bases;
mkdir calib_lips;
mkdir calib_figure_outputs;
mkdir calib_database_outputs;

num = 1; %counter

for b = 1:length(min_kp_size2_calib) % minimum knickpoint size2
    for a = 1:length(min_kp_size1_calib) % minimum knickpoint size1
        for c = 1:length(lumping_search_distance_calib)  % lumping window
            for l = 1:length(SG_smoothing_calib)   % smoothing window
            
            iteration_number = num % tells you progress (this takes a while)
            
            %% Name Parameters
            
            % I change only the smoothing windows and the minimum
            % knickpoint sizes
            
            % other parameters are fairly constant and have pretty clear
            % affects on results (e.i. if you increast the min_trib_size
            % you only parts of stream reaches at larger drainage areas)
            
            
            %% Start picking knickpoints
            % First Step, Reorganize chiplot and streamobj array so each tributary can be analyzed seperately
            current_basin_chiplot = AOI_chiplot{1}; % chiplot structure array for calibration basin
            current_basin_streamobj = AOI_STR_largest{1};
            
            
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
            
            
            % next step is just to keep naming consistent between
            % calibration script and the processing script that uses the
            % same function (SLIGHTLY DIFFERENT IN CALIB. SCRIPT BC ONLY 
            % ONE BASIN ANALYZED AT A TIME)    

            
            % save these in cell arrays for each basin (Would contain
            % nested cell array if more than 1 basin is analyzed)
            elev_stored_current = elev_stored_trimmed ;
            chi_stored_current = chi_stored_trimmed ;
            eastingUTM11_stored_current = eastingUTM11_stored_trimmed ;
            northingUTM11_stored_current = northingUTM11_stored_trimmed ;
            upstream_DA_stored_current = upstream_DA_stored_trimmed ;
            distance_stored_current = distance_stored_trimmed ;
            % each of these cell arrays contains a cell for each tributary
            % in the basin: ex: elev_stored_current contains cells {j} of
            % all of the elevation values for tributary {j}
            
            % empty the cell arrays after each stream iteration or else data from
            % large catchments with lots of tributaries (j) are not fully overwritten by small streams
            % with less tributaries ({j} isn't as long, so less cells)
            
            % step below clears cell arrays that contain knickzone
            % information. needs to be done because one tributary may have
            % less knickzones than the next, and then each cell won't get
            % totally overwritten
            kp_magnitude_matrix_final={};
            kp_chi_dist_final = {};
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

            % Preform knickpoint picking on each tributary, cycle through all of the
            % tributaries and save the results.  This will be a loop that iterates for
            % the number of tributaries in the basin.  We will interpolate each
            % tributary, de-trend it, differentiate it, and select knickpoints

            number_of_tributaries = length(elev_stored_current);
            % record number of tributaries

            for j = 1:number_of_tributaries
                % get the information for the tributary
                current_chi = chi_stored_current{j};
                %verify if this chi vector is correct -- if it contains unusual chi
                %values, continue to next tributary
                if nanmax(current_chi) > 1e10
                    continue
                end
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
                E4_detrend_smooth_differentiate_calib(current_chi_resampled, detrended_elev_current,...
                l,SG_smoothing_calib,sgolayfilt_order,smoothing_option);
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
                 E7_calib_min_kz_height_1_prelumping(kps_lips_length, smooth_detrended_elev_current,...
                 kps_cells_current_trib,base_cells_current_trib,a,min_kp_size1_calib);

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
                E8_calib_kz_lumping_filter(smooth_detrended_elev_current, sig_kps_bases,...
                sig_kps_cells,current_chi_resampled,current_chi,current_distance,...
                lumping_search_distance_calib,c);
            
            
                %% 3rd filter KZ magnitude filter post lumping
                % calculate the magnitude of the knickpoints after lumping (knickpoints
                % that lumped together will have a larger magnitude now)
                % eliminate knickzones that are still too small
                
                % option to filter based on minimum knickzone magnitude or 
                % minimum knickzone relief
                
                [sig_kps_cells_lumped_filtered, sig_kps_bases_lumped_filtered,...
                kp_magnitude_matrix_lumped_filtered] = ...
                E9_calib_min_kz_height2_post_lumping(smooth_detrended_elev_current, sig_kps_cells_lumped,...
                sig_kps_bases_lumped,min_kp_size2_calib,b,kp_magnitude_filter_option,...
                kp_relief_filter_option,current_elev_resampled);
            
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

                % REMEMBER TO SAVE EVERYTHING FOR THE CURRENT TRIBUTARY {j}
                
                % save attributes for knickzone lips found
                kp_up_da_tribs_lips{j} = kp_up_da;
                kp_distance_tribs_lips{j} = kp_distance;
                kp_easting_tribs_lips{j} = kp_easting;
                kp_northing_tribs_lips{j} = kp_northing;
                kp_chi_tribs_lips{j} = kp_chi;
                kp_elev_tribs_lips{j} = kp_elev;

                % save attributes for knickzone bases found
                kp_up_da_tribs_bases{j} = kp_up_da_bases;
                kp_distance_tribs_bases{j} = kp_distance_bases;
                kp_easting_tribs_bases{j} = kp_easting_bases;
                kp_northing_tribs_bases{j} = kp_northing_bases;
                kp_chi_tribs_bases{j} = kp_chi_bases;
                kp_elev_tribs_bases{j} = kp_elevation_bases;

                % save attributes for knickzone geometry
                kp_magnitude_matrix_tribs{j}=kp_magnitude_matrix_final;
                kp_length_tribs{j}=kp_dist_final;
                kp_face_slope_tribs{j} =kp_face_slope_final;
                kp_relief_tribs{j} = knickpoint_relief; %(store knickpoint relief measure)
                
                % save slope of current tributary
                beta_all_tribs{j} = beta_current_trib;  % for each tributary
                
                % Save all the important attributes in cell arrays, each cell representing the current tributary analyzed
                % Save attributes for current basin being analyzed 
                % (only one basin analyzed for calibration)
                kp_up_da_stored_lips{1} = kp_up_da_tribs_lips;
                kp_distance_stored_lips{1} = kp_distance_tribs_lips;
                kp_easting_stored_lips{1} = kp_easting_tribs_lips;
                kp_northing_stored_lips{1} = kp_northing_tribs_lips;
                kp_chi_stored_lips{1} = kp_chi_tribs_lips;
                kp_elev_stored_lips{1} = kp_elev_tribs_lips;

                % Knickpoint Bases
                kp_up_da_stored_bases{1} = kp_up_da_tribs_bases;
                kp_distance_stored_bases{1} = kp_distance_tribs_bases;
                kp_easting_stored_bases{1} = kp_easting_tribs_bases;
                kp_northing_stored_bases{1} = kp_northing_tribs_bases;
                kp_chi_stored_bases{1} = kp_chi_tribs_bases;
                kp_elev_stored_bases{1} = kp_elev_tribs_bases;

                % Knickpoint measurements
                kp_magnitude_matrix_stored{1}=kp_magnitude_matrix_tribs;
                kp_length_stored{1}=kp_length_tribs;
                kp_face_slope_stored{1} =kp_face_slope_tribs;
                kp_relief_stored{1} = kp_relief_tribs; %(store knickpoint relief measure)
                
                % save steepness of tributaries within basin
                beta_all_tribs_all_basins{1} = beta_all_tribs;  % for each basin
            end


            fprintf('\n');

            % All the knickpoints are saved in the above cell arrays
            % Organize all knickpoints into a lists (unpack the cell arrays and
            % organize these attributes into columns in a matrix so they 
            % can be coverted into a table.csv file)

            % So take all the knickpoints that are stored in the different cells of
            % the cell arrays and move them all into 1 list
            
            % use function 'kz_attributes_unpack'
    [all_kp_elevations, all_kp_chi,all_kp_easting,all_kp_northing,...
    all_kp_upstream_DA,all_kp_distance,all_kp_magnitude,all_kp_slope,...
    all_kp_length_distance,all_kp_relief,all_kp_elevations_bases,...
    all_kp_chi_bases,all_kp_easting_bases,all_kp_northing_bases,...
    all_kp_upstream_DA_bases,all_kp_distance_bases] = ...
    E12_calib_kz_attributes_unpack(kp_elev_stored_lips, kp_chi_stored_lips,...
    kp_easting_stored_lips, kp_northing_stored_lips,...
    kp_up_da_stored_lips, kp_distance_stored_lips, kp_magnitude_matrix_stored,...
    kp_face_slope_stored,kp_length_stored,kp_relief_stored,kp_elev_stored_bases,...
    kp_chi_stored_bases,kp_easting_stored_bases,kp_northing_stored_bases,...
    kp_up_da_stored_bases,kp_distance_stored_bases);

            % NEXT
% find the tributary that each knickpoint came from
% this is tricky because some tributaries have multiple knickpoints, some
% have none. Also store the steepness of the reference profile for each
% tributary. Store the basin id # for each knickpoint (which basin did the
% knickpoint come from), and store the ksn from SA plot, ks from chiplot,
% and m/n from chiplot for each basin.

% use function 'trib_id_num_for_kp'
[all_trib_ids, stream_id, ksn_SA_list, ks_chiplot_list,...
    m_n_basinwide_list, all_ks_trib_list] = ...
    E13_calib_trib_id_num_for_kp(beta_all_tribs_all_basins, kp_elev_stored_lips,theta_bf,...
    ks_chiplot,ksn_SA_store);


% record a running list to give each knickzone found an id # (kp_id_num)
% record parameters used to make the knickzone selections

% use function, 'save_parameters'
[kp_id_num,sgolay_window_size, sgolay_polynomial, max_kp_search_dist,...
    min_kp_size_pre_combination,min_kp_size_post_combination,...
    min_kp_slope_d,min_trib_size_cells,min_DA] = ...
    E14_calib_save_parameters(all_kp_elevations, SG_smoothing_calib,l,...
    lumping_search_distance_calib,c,min_kp_size1_calib,a,...
    min_kp_size2_calib,b,min_kp_slope,Min_trib_size,sgolayfilt_order,...
    Min_DA_threshold);
% NOTE ^^ in calibration script, parameters will change their values after
% each iteration in the a, b, c, and l loops!

%% write .csv file where each row represents a knickzone found, and each 
% column is the following attributes (saves information for each knickzone):

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

% DON'T do this for calibration script, because this would sloooooow down
% the runtime! but we SHOULD add this to the normal knickzone selection
% script (Read Below:)

% use these attributes, N-E-Projection to write shapefile? ^^ would add to
% function below! (just need to add input which is the path to the long
% profile figure generated from script: 'trib_figure_saver.m')

% would want to add 1 more column containing the PATH to the longitudinal
% profile figure of the tributary: so in ArcMAP the user could click on the
% plotted point and use the hyperlink lightning bolt tool to pop up the
% longitudinal profile!  SUPER USEFUL!

%% make the database that stores knickzone lip/base attributes listed above

            % knickpoint lips
            Knickpoint_Database_Lips = horzcat (kp_id_num,stream_id,...
                ksn_SA_list, ks_chiplot_list, m_n_basinwide_list,all_trib_ids,...          
                all_ks_trib_list,all_kp_chi,all_kp_elevations,all_kp_magnitude,...
                all_kp_relief,all_kp_easting,all_kp_northing,all_kp_upstream_DA,...
                all_kp_distance,all_kp_slope,all_kp_length_distance,...
                sgolay_window_size, sgolay_polynomial, max_kp_search_dist,...
                min_kp_size_pre_combination,min_kp_size_post_combination,...
                min_kp_slope_d,min_trib_size_cells,min_DA);
            

            %knickpoint bases (database containing all information about the bottom of
            %the knickpoints)
            Knickpoint_Database_Bases = horzcat (kp_id_num,stream_id,...
                ksn_SA_list, ks_chiplot_list, m_n_basinwide_list,all_trib_ids,...          
                all_ks_trib_list,all_kp_chi_bases,all_kp_elevations_bases,...
                all_kp_magnitude,all_kp_relief,all_kp_easting_bases,...
                all_kp_northing_bases,all_kp_upstream_DA_bases,...
                all_kp_distance_bases,all_kp_slope,all_kp_length_distance,...
                sgolay_window_size, sgolay_polynomial, max_kp_search_dist,...
                min_kp_size_pre_combination,min_kp_size_post_combination,...
                min_kp_slope_d,min_trib_size_cells,min_DA);
            
            
            % save parameters for current parameter combination
            Knickpoint_Database_Lips_all_params{num} = Knickpoint_Database_Lips;
            
            % save parameters for current parameter combination
            Knickpoint_Database_Bases_all_params{num} = Knickpoint_Database_Bases;
                        
            num = num+1; % increase counter for iteration!
            
            % clear everything except for the loop number, the data you
            % had to load, parameter choices from parameters script and 
            % data useful for generating plots and figures in the next script
            % also the loop counter 'num'
            clearvars -except a b c l AOI_chiplot AOI_STR_largest current_basin_streamobj...
                current_basin_chiplot SG_smoothing_calib lumping_search_distance_calib...
                min_kp_size1_calib min_kp_size2_calib AOI_DEM Min_trib_size...
                sgolayfilt_order error_radius Calibration_option...
                kp_magnitude_filter_option kp_relief_filter_option... 
                KZ_lips_calib_fname KZ_bases_calib_fname min_kp_slope calibration_snapping_tolerance...
                kp_plot_size KZ_calib_northing_column_num KZ_calib_easting_column_num...
                KZ_calib_relief_column_num do_you_have_calibration_KZ_bases smoothing_option...
                theta_bf_option theta_ref num theta_bf ks_chiplot ksn_SA_store Min_DA_threshold...
                Knickpoint_Database_Lips_all_params Knickpoint_Database_Bases_all_params...
            
            end
        end
    end
end