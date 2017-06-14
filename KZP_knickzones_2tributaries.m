%% (2) Iterate through all basin tributaries and extract knickpoints
%
fprintf(1,'KZP identifying knickzones step 2 of 4: Extracting knickpoints for all streams.\nAt basin # of %d: ', number_of_basins);
for i = 1:number_of_basins
    fprintf('%d, ', i);
    if mod(i,10) == 0
        fprintf('\nAt basin # of %d: ', number_of_basins);
    end
    % First Step, Reorganize chiplot and streamobj array so each tributary can be analyzed seperately
    
    % within SCI_example_LiDAR1m_STR_5e3_S_chiplot.chi, all of the chi values
    % of the stream network are listed.  Each stream segment (or tributary) is separated by a
    % NaN value which is the confluence between that tributary and the larger tributary/trunkstream.
    %
    % Cell 1 is the headwaters of the trunk
    % stream.  The first NaN is the mouth of the trunkstream.  The cell after
    % that is the headwaters of the first tributary.  The next NaN is the
    % confluence between that tributary and the trunk stream.  ECT for all tributaries....
    
    % you can use .elev, .chi, .dist to reorganize the structure array, here we just
    % chose to use .chi,
    
    current_basin_chiplot = AOI_STR_S_chiplot{i};
    current_basin_streamobj = AOI_STR_streams_dbasins_unique{i};
    
    idx_confluences = ~isnan(current_basin_chiplot.chi);  % sets cells that are NaN's to 0   (the confluences are where cells within .chi = NaN)
    
    confluences = find(idx_confluences == 0) ;  % FINDS NAN VALUES (CONFLUENCES)
    
    confluences = [1;confluences];    % GETS CELL LOCATIONS OF WHERE CONFLUENCES ARE   % Need to add the first cell in the structure array to the NaN's (start of trunkstream)
    
    % We want to arrange this into a cell array, where each cell within the
    % cell array is the chi information for only one tributary.  We will do
    % this for the .elev informaion, the .x information and the .y
    % information... ect
    % Each one of those variables is organized with the same format within the
    % structure array, so what works for .chi works for .elev, .x, and .y likewise.
    
    tributaries= {};  % create empty cell array to store tributaries in
    % ( also erases tributary information from previous loops) streams with
    % less tributaries will not overwrite all the cells created from
    % streams with more tributaries as this loop goes from basin to basin.
    % we want to make sure information from one basin doesn't get carried
    % over into the wrong cell array
    
    chi_stored = {}; % erase previous data
    elev_stored = {}; % erase previous data
    eastingUTM11_stored = {}; % erase previous data
    northingUTM11_stored={} ;% erase previous data
    upstream_DA_stored={} ;% erase previous data
    distance_stored = {};
    
    %  This starts from the headwaters of the trunk stream (cell 1), measures to the next NaN value (confluence)
    for j = 1:(length(confluences)-1);
        tributaries{j} = [confluences(j):1:confluences(j+1)] ;
        % This stores the cells between the each confluence
        %(from the headwaters to the mouth of the tributary). The trunk stream
        % is included. The values are stored in tributaries{j}. Each cell in
        % tributaries{j} contains information for 1 tributary (including ts, which is just the first cell in "tributaries")
        
        % Use the cells for each tributary found in the loop above to index and
        % sort all the .chi information into a cell array.  Each cell in chi_stored
        % contains the .chi information for first the trunk stream (chi_stored{1}) and then the
        % tributaries
        chi_stored{j} = current_basin_chiplot.chi(tributaries{j})  ;
        % Do the same thing for the .elev data
        elev_stored{j} = current_basin_chiplot.elev(tributaries{j});
        % Do the same thing for the .x data (easting UTM 11)
        eastingUTM11_stored{j} = current_basin_chiplot.x(tributaries{j})   ;
        % Do the same thing for the .y data (northing UTM 11)
        northingUTM11_stored{j} = current_basin_chiplot.y(tributaries{j}) ;
        % drainge area info
        upstream_DA_stored{j} = current_basin_chiplot.area(tributaries{j}) ;
        % distance infor
        distance_stored{j} = current_basin_chiplot.distance(tributaries{j}) ;
        % done
    end
    
    % get rid of tributaries that are too small
    % find the tributaries that are too small and mark them so we can remove them
    % this step is not necessary, but can help you generate clearer
    % figures, and increase run speed for the rest of the code.
    
    % we usually aren't interested in tributaries that are less than a
    % certian size.  Here we remove tributaries shorter than a certain
    % length of cells.
    %
    % This isn't the best way to do this because # of
    % cells isn't exactly the length of the stream or the drainage area.
    % But this is simple and works suffistream_namely.  We can upgrade this if
    % necessary.
    
    for j = 1:length(tributaries);
        % Mark cells for short tributaries in the stored elevation cell array (same
        % cells as chi, .x, .y, ect..
        if length(chi_stored{j}) < KZP_parameters.min_trib_size ; % measure the length of the chi vector for each tributary. if that vector is too small we ignore that stream segment
            chi_stored{j} = NaN ;
            elev_stored{j} = NaN;
            eastingUTM11_stored{j} = NaN;
            northingUTM11_stored{j} = NaN;
            upstream_DA_stored{j} = NaN;
            distance_stored{j} = NaN;
        end
    end
    % this is a good method though to avoid having to increase the min drainage area for
    % streams we care about.  This way we can keep a low minimum drainage
    % area, 1e4m^2, but get rid of tributaries that are very short.
    % Headwaters are not far from the confluence with another tributary
    
    % now index out tributaries that are too small (tributaries within
    % "chi_stored" (for example) that we marked NaN
    for j = 1:length(tributaries);
        if isnan(chi_stored{j})
            chi_stored{j} = [];
            elev_stored{j} = [];
            eastingUTM11_stored{j} = [];
            northingUTM11_stored{j} = [];
            upstream_DA_stored{j} = [];
            distance_stored{j} = [];
        end
    end
    % Remove the Empty cells from each cell array (we turned all the
    % tributaries that were too small into empty cells)
    
    elev_stored_trimmed = elev_stored(~cellfun(@isempty, elev_stored))   ;
    chi_stored_trimmed = chi_stored(~cellfun(@isempty, chi_stored))    ;
    eastingUTM11_stored_trimmed = eastingUTM11_stored(~cellfun(@isempty, eastingUTM11_stored))    ;
    northingUTM11_stored_trimmed = northingUTM11_stored(~cellfun(@isempty, northingUTM11_stored)) ;
    upstream_DA_stored_trimmed = upstream_DA_stored(~cellfun(@isempty, upstream_DA_stored)) ;
    distance_stored_trimmed = distance_stored(~cellfun(@isempty, distance_stored)) ;
    % new cell arrays with just tributaries of interest  ^^^ (for information;
    % of interest:  .chi, .elev, .x, .y)
    
    % save these in cell arrays for each basin
    elev_stored_trimmed_basin{i} = elev_stored_trimmed ;
    chi_stored_trimmed_basin{i} = chi_stored_trimmed ;
    eastingUTM11_stored_trimmed_basin{i} = eastingUTM11_stored_trimmed ;
    northingUTM11_stored_trimmed_basin{i} = northingUTM11_stored_trimmed ;
    upstream_DA_stored_trimmed_basin{i} = upstream_DA_stored_trimmed ;
    distance_stored_trimmed_basin{i} = distance_stored_trimmed ;
    
    % now each basin (i), has an array of tributaries (j) that we can
    % study seperately.
end
fprintf('\n');

% Postprocessing knickpoints
fprintf(1, 'KZP identifying knickzones step 2 of 4: Postprocessing knickpoints for all tributaries.\nAt basin # of %d: ', ...
    number_of_basins);
for i = 1:number_of_basins
    fprintf('%d, ', i);
    if mod(i,10) == 0
        fprintf('\nAt basin # of %d: ', ...
            number_of_basins);
    end
    
    chi_stored_current = chi_stored_trimmed_basin{i}; % chi info for each tributary within the stream {i}
    elev_stored_current  = elev_stored_trimmed_basin{i}; % elev info for each tributary within the stream {i}
    eastingUTM11_stored_current  = eastingUTM11_stored_trimmed_basin{i}; % ect..
    northingUTM11_stored_current  = northingUTM11_stored_trimmed_basin{i};
    upstream_DA_stored_current  = upstream_DA_stored_trimmed_basin{i};
    distance_stored_current  = distance_stored_trimmed_basin{i};
    
    % empty the cell arrays after each stream iteration or else data from
    % large catchments with lots of tributaries (j) are not fully overwritten by small streams
    % with less tributaries (j isn't as long, so less cells)
    kp_magnitude_matrix_final={};
    kp_chi_length_final = {};
    kp_face_slope_chi_final = {};
    kp_up_da= {};
    kp_distance= {};
    kp_easting = {};
    kp_northing = {};
    kp_chi = {};
    kp_elev = {};
    knickpoint_relief = {}; %(erase knickpoint relief from previous iterations)
    
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
    
    number_of_tributaries = length(elev_stored_current);
    
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
        
        % Now We are ready to Interpolate
        % get chi stepsize
        for d = 2:length(current_chi_trib)
            current_stepsize(d) = current_chi_trib(d) - current_chi_trib(d-1);
        end
        % takes the difference between each successive point and stores them in
        % current_stepsize
        current_stepsize(1) = 0; % add cell to the begining because we lost a cell in that loop
        
        %the next line calculates averaged stepsize for chi:
        current_stepsize_mean = mean(abs(current_stepsize));
        
        % make a new vector with average stepsize for resampling/interpolation
        current_chi_resampled = min(current_chi_trib):current_stepsize_mean:max(current_chi_trib);
        
        %perform 1D interpolation
        current_elev_resampled = interp1(current_chi_trib,current_elev_trib,current_chi_resampled,'linear');
        
        % need to flip these and transpose them because they come out upside down
        current_chi_resampled = flipud(current_chi_resampled');
        current_elev_resampled = flipud(current_elev_resampled');
        
        %   Work with resampled data
        % second step: calculate linear regression and detrend the data
        mx_b = (max(current_elev_resampled)-min(current_elev_resampled))/(max(current_chi_resampled) - min(current_chi_resampled)); % calculate linear regression
        blf_y_current_tribs = (mx_b*current_chi_resampled + (min(current_elev_resampled))); % get the y coordinates for the trendline
        
        % save slope of chiplot
        beta_current_trib = mx_b;
        detrended_elev_current = current_elev_resampled - blf_y_current_tribs;  % perform detrending on resampled data
        
        % Smooth Data (THIS IS WHERE WE COULD PUT A DIFFERENT FILTER/SMOOTHER)
        %perform smoothing
        if length(current_chi_resampled) < 10
            continue
        end
        if license('test', 'curve_fitting_toolbox') == 1
            smooth_detrended_elev_current = smooth(current_chi_resampled,detrended_elev_current,KZP_parameters.smoothing_window,'sgolay');
        else
            smooth_detrended_elev_current = sgolayfilt(detrended_elev_current, KZP_parameters.sgolayfilt_order, KZP_parameters.smoothing_window);
        end
        % Differentiate
        diff_smooth_detrended_elev_current = diff(smooth_detrended_elev_current); % differentiate
        
        local_max=[];
        local_min=[];
        
        % Now start picking out some knickpoints
        % first check if there are potential knickpoints
        %         if max(diff_smooth_detrended_elev_current) < 1
        %             continue;
        %         end
        
        % find inflections (zeros in the differentiated, smoothed detrended curve)
        n = length(diff_smooth_detrended_elev_current);
        for k = 1:n-1
            if diff_smooth_detrended_elev_current(k)>=0 && diff_smooth_detrended_elev_current(k+1) <0
                local_max(k) = 1;      % < -- stores potential knickpoint lips
                % scans "diff_smooth_detrend_elev_current_trib" from headwaters to
                % mouth.  if the cell closer to the headwaters is + and the
                % adjastream_name next downstream cell is -.
                % We mark the cell as 1, a potential knickpoint.  This is an
                % inflection in detrended gradient from positive to
                % negative. (convexity, concave down)
            end
            
            if diff_smooth_detrended_elev_current(k)<=0 && diff_smooth_detrended_elev_current(k+1) >0
                local_min(k) = 1;     % < -- stores potential knickpoint bases
                % scans "diff_smooth_detrend_elev_current_trib2" from headwaters to
                % mouth.  if the cell closer to the headwaters is - and the
                % adjastream_name next downstream cell is +.
                % We mark the cell as 1, a base (knickpoint base)
                
                % this is the other type of inflection point
                % (convexity, concave up)
            end
        end
        
        % Organize inflections so that upstream knickpoint lips are compared to downstream knickpoint
        % bases. Then index them within the chi/elev matricies to get their
        % positions/elevations
        % find the cells within the current tributary that are potential knickpoints
        kps_cells_current_trib = find(local_max== 1)'; % cells containing knickpoint lips
        base_cells_current_trib = find(local_min== 1)'; % cells containing knickpoint bases
        
        kps_lips_length = length(kps_cells_current_trib);  % find how many concave down inflections are found
        kps_base_length = length(base_cells_current_trib); % find how many concave up inflections are found
        
        % NOW ORGANIZE THESE MATRICIES
        % reorganize matrices of inflection points so you know you are comparing
        % the upstream knickpoint lip to the downstream knickpoint base
        
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
        
        %
        % CASE 1
        % downstream        <------              Upstream
        %  knickpoint - base - knickpoint - base
        % not what we want.. make adjustments
        if ~isempty(kps_cells_current_trib)
            if kps_lips_length == kps_base_length && min(kps_cells_current_trib)> min(base_cells_current_trib) ;
                % if the stream starts out with the first "0" detrend gradient as a
                % base (near headwaters)^^ and then ends with the
                % last "0" detrend gradient as a knickpoint (near mouth
                % of stream)
                
                base_cells_current_trib = base_cells_current_trib(2:kps_base_length,1);
                % ^^ then we need to remove the first base and add the mouth of the stream so all the
                % knickpoint lips will be compared to their respective
                % knickpoint bases, and the last knickpoint near the
                % mouth of the stream will be compared to the mouth of the stream
                
                base_cells_current_trib = [base_cells_current_trib;n];
                % Lastly, we must add the mouth of the stream, so the last
                % knickpoint will be compared to the mouth of the stream "n"
            end
        end
        
        % End with:
        % downstream        <------              Upstream
        %  mouth stream - knickpoint - base - knickpoint
        % ^^ that will work
        
        %
        % CASE 2
        % downstream        <------              Upstream
        %  knickpoint - base - knickpoint - base - knickpoint
        % this needs to be fixed too, the knickpoint at the end of the stream needs
        % a reference base
        
        if kps_lips_length> kps_base_length ;
            base_cells_current_trib = [base_cells_current_trib;n];
        end
        % if there are more knickpoints than bases, this loop
        % adds the last cell
        % point (mouth of the stream 'n') to the end of the bases vector to make them
        % the same size. So we compare the elevation of the last knickpoint
        % (furthest downstream) to the elevation of the mouth of the stream.
        
        % End with:
        % downstream        <------              Upstream
        %  mouth stream - knickpoint - base - knickpoint
        %^^ that will work
        
        %
        % CASE 3
        % downstream        <------              Upstream
        %  base - knickpoint - base - knickpoint - base
        % ^ this won't work because the first inflection is a base
        
        if kps_lips_length < kps_base_length;
            base_cells_current_trib = base_cells_current_trib(2:kps_base_length(1),1);
        end
        
        % if there are more base than knickpoints, this loop removes the first
        % base (the extra base that occurs near the headwaters of the stream)
        
        % End with:
        % downstream        <------              Upstream
        %  base - knickpoint - base - knickpoint
        % ^ that will work
        
        % Case 4
        % downstream        <------              Upstream
        %                     base
        % ^ here there is only 1 base and no knickpoints. Just remove the
        % base because we aren't interested in that
        
        if kps_base_length == 1 && kps_lips_length == 0
            base_cells_current_trib = [];
        end
        %
        % ALWAYS WANT EITHER:
        % downstream        <------              Upstream
        %  base - knickpoint - base - knickpoint
        % or
        % downstream        <------              Upstream
        %  mouth stream - knickpoint - base - knickpoint
        
        %  All of the knickpoint cells are organized so the upstream knickpoint lips are compared to a respective downstream knickpoint base
        % First Filter (knickpoint magnitude)
        
        % Now We are ready to use the Knickpoint magnitude Filter,
        % each knickpoint lip elevation
        % can be compared to the base directly downstream of it.
        % if the magnitude of detreneded elevation change doesn't exceed a
        % threshhold, we disregard it
        
        n2 = kps_lips_length;  % for all the potential knickpoints found
        for k =1:1:n2;
            if smooth_detrended_elev_current(kps_cells_current_trib (k))-...
                    smooth_detrended_elev_current (base_cells_current_trib (k)) < KZP_parameters.min_kp_size1 ;
                % see uf the elevation difference between the upstream knickpoint lip and downstream base exceeds a minimu knickpoint size
                
                kps_cells_current_trib(k) = NaN;  % sets the knickpoints that are too small to NaN
                base_cells_current_trib(k) = NaN;  % does the same thing with the bases corresponding to those knickpoints
            end
        end
        
        idx_kps_tribs= ~isnan(kps_cells_current_trib); % index out those inflections that were too small
        sig_kps_cells= kps_cells_current_trib(idx_kps_tribs);  %cells of potential knickpoints
        
        idx_base= ~isnan(base_cells_current_trib);  % index out those inflections that were too small
        sig_kps_bases = base_cells_current_trib(idx_base);  %cells of potential respective bases
        
        % done with the first knickpoint magnitude filter
        
        %%        Second Filter
        % LUMP KNICKPOINTS: looks for closely spaced
        % knickpoints and lumps them together into a bigger
        % knickpoint.
        
        % If knickpoints are within a certain chi
        % distance from one another
        % their magnitudes get summed and the downstream knickpoint
        % lip is erased. also the upstream knickpoint base is
        % erased.  this is to compare the elevation drop from the
        % upstream most knickpoint lip, and the downstream most
        % base.
        
        % this filter is a bit tough to follow, but it works
        if KZP_parameters.chi_lump_option == 1
            kp_elev_bases = smooth_detrended_elev_current(sig_kps_bases); % get the detrended elevation of all the bases
            kp_magnitude_matrix_1 = smooth_detrended_elev_current(sig_kps_cells)...
                - kp_elev_bases ;  % get the magnitude of knickpoints (lip - base)
            
            
            sig_kps_chi = current_chi_resampled(sig_kps_cells); % get the chi coordinate for each knickpoint
            
            k = 1:(length(sig_kps_chi));
            k = fliplr(k) ; % cycle through tributary from downstream to upstream
            lumpcell_lips = zeros(length(sig_kps_chi),1);  % make matrix to store the cells where we lump knickpoints
            lumpcell_base = zeros(length(sig_kps_chi),1);  % make matrix to store the cells where we lump bases
            
            p =k(1:length(k)-1); % scan from downstream to upstream
            for t = 1:length(p)
                m=p(t);  % remove the last index term (1)
                if sig_kps_chi(m-1) - sig_kps_chi(m) < KZP_parameters.lumping_search_distance && kp_elev_bases(m) < kp_elev_bases(m-1)
                    % if the knickpoints are close together (less than 150 chi)
                    % and the base of the downstream knickpoint is at a lower elevation than that of the base of the upstream knickpoint
                    
                    lumpcell_lips(m) = NaN;  % set the downstream kp to NaN
                    lumpcell_base(m-1) = NaN  ;   % set the base for the upstream kp to NaN
                end
                %( we want to compare the upstream knickpoint elev to
                %the downstream base elev to get the total elev drop)
                
            end
            
        end
        
        %% if lumping knickzones by distance upstream
        if KZP_parameters.distance_upstream_lump_option == 1
            
            kp_elev_bases = smooth_detrended_elev_current(sig_kps_bases); % get the detrended elevation of all the bases
            kp_magnitude_matrix_1 = smooth_detrended_elev_current(sig_kps_cells)...
                - kp_elev_bases ;  % get the magnitude of knickpoints (lip - base)
            
            % find the distance upstream for the potential knickzone boundaries
            kp_cells_d = [];
            
            kp_cells_chi_equispaced_filtering = current_chi_resampled(sig_kps_cells); % chi coordinates for knickpoint lips
            % reference the original chi dataset, find where the values are closest to
            % the knickpoint chi coordinates above
            
            for k = 1:length(kp_cells_chi_equispaced_filtering)
                [~,idx] = min(abs(current_chi-kp_cells_chi_equispaced_filtering(k))) ;  % find cell with the minimum difference (usually pretty small)
                kp_cells_d(k) = idx; % these are the cells in the original data which contain the knickpoints
            end
            
            sig_kps_distance = current_distance_trib(kp_cells_d); % get the chi coordinate for each knickpoint
            
            k = 1:(length(sig_kps_distance));
            k = fliplr(k) ; % cycle through tributary from downstream to upstream
            lumpcell_lips = zeros(length(sig_kps_distance),1);  % make matrix to store the cells where we lump knickpoints
            lumpcell_base = zeros(length(sig_kps_distance),1);  % make matrix to store the cells where we lump bases
            
            p =k(1:length(k)-1); % scan from downstream to upstream
            for t = 1:length(p)
                m=p(t);  % remove the last index term (1)
                if sig_kps_distance(m-1) - sig_kps_distance(m) < lumping_distance_upstream && kp_elev_bases(m) < kp_elev_bases(m-1)
                    % if the knickpoints are close together (less than 150 chi)
                    % and the base of the downstream knickpoint is at a lower elevation than that of the base of the upstream knickpoint
                    
                    lumpcell_lips(m) = NaN;  % set the downstream kp to NaN
                    lumpcell_base(m-1) = NaN  ;   % set the base for the upstream kp to NaN
                end
                %( we want to compare the upstream knickpoint elev to
                %the downstream base elev to get the total elev drop)
                
            end
            
        end
        
        
        % remove the knickpoints that got lumped (nan's)
        
        idx_kps_tribs_2= ~isnan(lumpcell_lips);
        sig_kps_cells_lumped = sig_kps_cells(idx_kps_tribs_2);   % remove the knickpoint lip locations of knickpoints that got lumped
        idx_bases_tribs_2= ~isnan(lumpcell_base);
        sig_kps_bases_lumped  = sig_kps_bases(idx_bases_tribs_2);  % remove respective bases of knickpoints that got lumped
        
        % calculate the magnitude of the knickpoints after lumping (knickpoints
        % that lumped together will have a larger magnitude now)
        kp_magnitude_matrix_lumped_equispaced = smooth_detrended_elev_current(sig_kps_cells_lumped) - smooth_detrended_elev_current(sig_kps_bases_lumped);
        
        % plot option for debugging
        %figure(2)
        %clf
        %plot(current_chi_resampled,current_elev_resampled)
        %hold on
        %plot(current_chi_resampled(sig_kps_cells),current_elev_resampled(sig_kps_cells),'r.','Markersize',12)
        %plot(current_chi_resampled(sig_kps_cells_lumped),current_elev_resampled(sig_kps_cells_lumped),'k.')
        
        
        
        %% Apply 3rd filter, another minimum knickpoint size filter to get rid of knickpoints that didn't grow big enough after lumping
        if KZP_parameters.kp_magnitude_filter_option == 1
            min_kp_size2 = KZP_parameters.min_kp_size2_magnitude;
            for k = 1:length(sig_kps_cells_lumped); % for how many potential knickpoints we still have
                if kp_magnitude_matrix_lumped_equispaced(k) < KZP_parameters.min_kp_size2_magnitude % specify a 2nd minimum knickpoint size if you'd like (after summing closely spaced knickpoints together)
                    % THIS IS OPTIONAL BUT RECOMMENDED (BASICALLY IF YOU'D LIKE TO REMOVE SMALL KPS THAT DIDN"T LUMP UP INTO BIG KPS
                    kp_magnitude_matrix_lumped_equispaced(k) = NaN;  % remove the stored magnitude for knickpoints that are too small
                    sig_kps_bases_lumped(k) =NaN; % remove the stored base position for knickpoints that are too small
                    sig_kps_cells_lumped(k) = NaN; % remove the stored lips for knickpoints that are too small
                end
            end
            
            % remove knickpoints that aren't big enough even after lumping close
            % proximity knickpoints
            idx_kps_tribs_3= ~isnan(sig_kps_cells_lumped);
            sig_kps_cells_lumped_filtered = sig_kps_cells_lumped(idx_kps_tribs_3);   % index out knickpoints that after lumping still didn't get big enough
            
            idx_bases_tribs_3= ~isnan(sig_kps_bases_lumped);
            sig_kps_bases_lumped_filtered = sig_kps_bases_lumped(idx_bases_tribs_3);  % index out respective bases to those knickpoints
            
            idx_magnitude_matrix_2 = ~isnan(kp_magnitude_matrix_lumped_equispaced);
            kp_magnitude_matrix_lumped_filtered  = kp_magnitude_matrix_lumped_equispaced(idx_magnitude_matrix_2);  % magnitudes of those knickponts that didn't get big enough
        end
        
        % if you want to filter by knickzone relief rather than magnitude
        if KZP_parameters.kp_relief_filter_option == 1
            min_kp_size2 = KZP_parameters.min_kp_size2_relief;
            current_rel = current_elev_resampled(sig_kps_cells_lumped) - current_elev_resampled(sig_kps_bases_lumped);
            
            for k = 1:length(sig_kps_cells_lumped); % for how many potential knickpoints we still have
                if current_rel(k) < min_kp_size2_relief % specify a 2nd minimum knickpoint size if you'd like (after summing closely spaced knickpoints together)
                    % THIS IS OPTIONAL BUT RECOMMENDED (BASICALLY IF YOU'D LIKE TO REMOVE SMALL KPS THAT DIDN"T LUMP UP INTO BIG KPS
                    kp_magnitude_matrix_lumped_equispaced(k) = NaN;  % remove the stored magnitude for knickpoints that are too small
                    sig_kps_bases_lumped(k) =NaN; % remove the stored base position for knickpoints that are too small
                    sig_kps_cells_lumped(k) = NaN; % remove the stored lips for knickpoints that are too small
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
        
        % plot option for debugging
        %figure(2)
        %clf
        %plot(current_chi_resampled,current_elev_resampled)
        %hold on
        %plot(current_chi_resampled(sig_kps_cells_lumped_filtered),current_elev_resampled(sig_kps_cells_lumped_filtered),'g.')
        %plot(current_chi_resampled(sig_kps_bases_lumped_filtered),current_elev_resampled(sig_kps_bases_lumped_filtered),'r.')
        
        %% Final Filter (minimum steepness filter) OPTIONAL!
        
        % If the steepness across the knickpoint doesn't exceed some value you can
        % get rid of that knickpoint here.  It is not traditional "steepness",
        % we're looking at:  detrended_elev_drop/chi, where steepness would be elev/chi
        
        % elev/distance upstream is the slope of the knickzone
        
        % find the distance upstream for the potential knickzone boundaries
        kp_cells_dist = [];
        kp_bases_dist = [];
        
        kp_cells_chi_equispaced_slope = current_chi_resampled(sig_kps_cells_lumped_filtered); % chi coordinates for knickpoint lips
        % reference the original chi dataset, find where the values are closest to
        % the knickpoint chi coordinates above
        
        for k = 1:length(kp_cells_chi_equispaced_slope)
            [~,idx] = min(abs(current_chi-kp_cells_chi_equispaced_slope(k))) ;  % find cell with the minimum difference (usually pretty small)
            kp_cells_dist(k) = idx; % these are the cells in the original data which contain the knickpoints
        end
        
        kp_bases_chi_equispaced_slope = current_chi_resampled(sig_kps_bases_lumped_filtered); % chi coordinates for knickpoint lips
        % reference the original chi dataset, find where the values are closest to
        % the knickpoint chi coordinates above
        
        for k = 1:length(kp_bases_chi_equispaced_slope)
            [~,idx] = min(abs(current_chi-kp_bases_chi_equispaced_slope(k))) ;  % find cell with the minimum difference (usually pretty small)
            kp_bases_dist(k) = idx; % these are the cells in the original data which contain the knickpoints
        end
        
        
        kp_length_dist = current_distance_trib(kp_cells_dist)- current_distance_trib(kp_bases_dist);
        % calculate length of knickpoint
        kp_height = current_elev_trib(kp_cells_dist) - current_elev_trib(kp_bases_dist);
        
        kp_face_slope_d = atand(kp_height./kp_length_dist) ; % divide drop in elevation by change in chi
        % this gives us a slope in "change in detrended elevation/change in
        % chi", basically the steepness of the reach relative to the
        % average steepness of the stream
        
        % THis is an optional filter that would filter out
        % knickpoints with face slopes that are too gradual
        
        
        for k = 1:length(kp_length_dist) % for the number of potential knickpoitns we still have
            if kp_face_slope_d(k) < KZP_parameters.min_kp_slope;
                
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
        kp_magnitude_matrix_final{j}  = kp_magnitude_matrix_lumped_filtered(idx_magnitude_matrix_3);  % respective knickpoint magnitudes
        
        idx_kp_chi_length = ~isnan(kp_length_dist);
        kp_chi_length_final{j} = kp_length_dist(idx_kp_chi_length );                       % respective knickpoint lengths (in units of chi)
        
        idx_kp_face_slope_dist_space_detrended  = ~isnan(kp_face_slope_d);
        kp_face_slope_chi_final{j}= kp_face_slope_d(idx_kp_face_slope_dist_space_detrended);   % respective knickpoint slopes  detrended elev/chi
        
        % All the filtering is done! We have our knickpoints identified!
        % Time to re-reference knickpoints back to orignal data
        
        % find the chi values for each knickpoint in the equispaced data, use that
        % to find the corresponding cell with the same chi value in the original
        % data.  then once you have all those cells that contain the knickpoints,
        % you can extract all the rest of the attributes for each knickponit
        % (upstream da, distance, northing easting ect..)
        kp_cells = [];
        kp_cells_bases=[];
        
        kp_cells_chi_equispaced = current_chi_resampled(sig_kps_cells_final); % chi coordinates for knickpoint lips
        kp_cells_bases_chi_equispaced = current_chi_resampled(sig_kps_bases_final); % chi coordinates for knickpoint bases
        
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
        
        % calculate knickpoint relief (elevation drop)
        current_kp_lip_elev = current_elev(kp_cells);
        current_kp_bases_elev = current_elev(kp_cells_bases);
        knickpoint_relief{j} = current_kp_lip_elev - current_kp_bases_elev;
        
        
        % knickpoint lips
        kp_up_da{j} = current_upDA(kp_cells);  % drainage area upstream of knickpoint
        kp_distance{j} = current_distance(kp_cells); % knickpoint distance upstream
        kp_easting{j} = current_easting(kp_cells);
        kp_northing{j} = current_northing(kp_cells);
        kp_chi{j} = current_chi(kp_cells); % chi coordinate
        kp_elev{j} = current_elev(kp_cells); % elevation of knickpoint
        
        % knickpoint bases
        kp_up_da_bases{j} = current_upDA(kp_cells_bases);  % drainage area upstream of knickpoint
        kp_distance_bases{j} = current_distance(kp_cells_bases); % knickpoint distance upstream
        kp_easting_bases{j} = current_easting(kp_cells_bases);
        kp_northing_bases{j} = current_northing(kp_cells_bases);
        kp_chi_bases{j} = current_chi(kp_cells_bases); % chi coordinate
        kp_elevation_bases{j} = current_elev(kp_cells_bases); % elevation of knickpoint
        
        
        % Save all the important attributes in cell arrays, each cell representing the current tributary analyzed
        % Knickpoint lips
        kp_up_da_stored_lips{i} = kp_up_da;
        kp_distance_stored_lips{i} = kp_distance;
        kp_easting_stored_lips{i} = kp_easting;
        kp_northing_stored_lips{i} = kp_northing;
        kp_chi_stored_lips{i} = kp_chi;
        kp_elev_stored_lips{i} = kp_elev;
        
        % Knickpoint Bases
        kp_up_da_stored_bases{i} = kp_up_da_bases;
        kp_distance_stored_bases{i} = kp_distance_bases;
        kp_easting_stored_bases{i} = kp_easting_bases;
        kp_northing_stored_bases{i} = kp_northing_bases;
        kp_chi_stored_bases{i} = kp_chi_bases;
        kp_elev_stored_bases{i} = kp_elevation_bases;
        
        % Knickpoint measurements
        kp_magnitude_matrix_stored{i}=kp_magnitude_matrix_final;
        kp_chi_length_stored{i}=kp_chi_length_final;
        kp_face_slope_chi_stored{i} =kp_face_slope_chi_final;
        kp_relief_stored{i} = knickpoint_relief; %(store knickpoint relief measure)
        
        % save slope of current tributary
        beta_all_tribs{j} = beta_current_trib;  % for each tributary
        beta_all_tribs_all_basins{i} = beta_all_tribs;  % for each basin
    end
end

fprintf('\n');

% All the knickpoints are saved in the above cell arrays
% Organize all knickpoints into a lists (unpack the cell arrays and
% organize these into columns in a matrix)

% So take all the knickpoints that are stored in the different cells of
% the cell arrays and move them all into 1 list
number_of_basins2 = length(kp_elev_stored_lips);

for i = 1:number_of_basins2 ;% % for all stream basins
    %for all the knickpoint lips
    current_basin_elevation = kp_elev_stored_lips{i};
    current_basin_chi = kp_chi_stored_lips{i};
    current_basin_easting = kp_easting_stored_lips{i};
    current_basin_northing = kp_northing_stored_lips{i};
    current_basin_upstream_DA = kp_up_da_stored_lips{i};
    current_basin_distance = kp_distance_stored_lips{i};
    
    % measurements between lips and bases
    current_basin_magnitude = kp_magnitude_matrix_stored{i};
    current_basin_kp_slope_detrend = kp_face_slope_chi_stored{i};
    current_basin_kp_length_distance = kp_chi_length_stored{i};
    current_basin_kp_relief = kp_relief_stored{i}; %(unpack stored relief values)
    
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
    current_basin_all_kps_slope_detrend = [];
    current_basin_all_kps_length_distance = [];
    current_basin_all_kps_relief = [];
    
    %bases
    current_basin_all_elevation_bases = [];
    current_basin_all_chi_bases = [];
    current_basin_all_easting_bases = [];
    current_basin_all_northing_bases = [];
    current_basin_all_upstream_DA_bases = [];
    current_basin_all_distance_bases = [];
    
    
    for j = 1:length(current_basin_elevation);  % for each tributary, vertically concatinate all the knickpoints from each tributary within the basin
        %lips
        current_basin_all_kps_elevation = vertcat(current_basin_all_kps_elevation,current_basin_elevation{j});
        current_basin_all_kps_chi = vertcat(current_basin_all_kps_chi,current_basin_chi{j});
        current_basin_all_kps_easting = vertcat(current_basin_all_kps_easting,current_basin_easting{j});
        current_basin_all_kps_northing = vertcat(current_basin_all_kps_northing,current_basin_northing{j});
        current_basin_all_kps_upstream_DA = vertcat(current_basin_all_kps_upstream_DA,current_basin_upstream_DA{j});
        current_basin_all_kps_distance = vertcat(current_basin_all_kps_distance,current_basin_distance{j});
        %measurements
        current_basin_all_kps_magnitude = vertcat(current_basin_all_kps_magnitude,current_basin_magnitude{j});
        current_basin_all_kps_slope_detrend = vertcat(current_basin_all_kps_slope_detrend,current_basin_kp_slope_detrend{j});
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
    all_basins_all_kps_slope_detrend{i} = current_basin_all_kps_slope_detrend;
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
all_kp_slope_detrend = [];
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
    all_kp_easting = vertcat(all_kp_easting,all_basins_all_kps_easting{i}); %9
    all_kp_northing = vertcat(all_kp_northing,all_basins_all_kps_northing{i}); %10
    all_kp_upstream_DA = vertcat(all_kp_upstream_DA,all_basins_all_kps_upstream_DA{i}); %11
    all_kp_distance = vertcat(all_kp_distance,all_basins_all_kps_distance{i}); %12
    % measurements
    all_kp_magnitude = vertcat(all_kp_magnitude,all_basins_all_kps_magnitude{i}); %7
    all_kp_slope_detrend = vertcat(all_kp_slope_detrend,all_basins_all_kps_slope_detrend{i}); % 15
    all_kp_length_distance = vertcat(all_kp_length_distance,all_basins_all_kps_length_distance{i}); %16
    all_kp_relief = vertcat(all_kp_relief,all_basins_all_kps_relief{i}); %8
    % bases
    all_kp_elevations_bases = vertcat(all_kp_elevations_bases,all_basins_all_kps_elevation_bases{i}); %6
    all_kp_chi_bases = vertcat(all_kp_chi_bases,all_basins_all_kps_chi_bases{i}); %5
    all_kp_easting_bases = vertcat(all_kp_easting_bases,all_basins_all_kps_easting_bases{i}); %9
    all_kp_northing_bases = vertcat(all_kp_northing_bases,all_basins_all_kps_northing_bases{i}); %10
    all_kp_upstream_DA_bases = vertcat(all_kp_upstream_DA_bases,all_basins_all_kps_upstream_DA_bases{i}); %11
    all_kp_distance_bases = vertcat(all_kp_distance_bases,all_basins_all_kps_distance_bases{i}); %12
    
end

%
% find the tributary that each knickpoint came from
% this is tricky because some tributaries have multiple knickpoints, some
% have none

for i = 1:number_of_basins2 % % for all stream basins
    current_basin_elevation = kp_elev_stored_lips{i};
    beta_list_current_basin = beta_all_tribs_all_basins{i};
    for j= 1:length(current_basin_elevation);
        trib_number= j;
        trib_id_matrix = [];
        beta_matrix = [];
        number_of_kps = length(current_basin_elevation{j}) ; % organize cell array so each tributary records the same id number  (how many knickpoints were there in that tributary)
        current_beta = beta_list_current_basin{j}; % get the beta used for that tributary
        for k = 1:number_of_kps ;%of cells as knickpoints, but the value in the cell is just "j" the tributary number
            
            trib_id_matrix(k) = j;
            beta_matrix(k) = current_beta; % create a matrix that has the used beta listed as many times as the number of knickpoints found
            
        end
        trib_array{j} = transpose(trib_id_matrix);  % cell array for basin containing each trib, but each trib just has trib number saved for as many knickpoints as they have
        beta_array{j} = transpose(beta_matrix);
    end
    all_trib_ids_each_basin = [];
    all_beta_each_basin = [];
    
    for j= 1:length(current_basin_elevation);
        all_trib_ids_each_basin = vertcat(all_trib_ids_each_basin,trib_array{j}); % concatinate each trib id within the basin
        all_beta_each_basin = vertcat(all_beta_each_basin,beta_array{j}); %same with beta
    end
    all_trib_ids_all_basins{i} = all_trib_ids_each_basin;  % store each trib id list for each basin in a different cell
    all_beta_all_basins{i}= all_beta_each_basin;
    
    
    % get the stream id number for each knickpoint too (which stream
    % did this come from)
    stream_id_matrix= [];
    
    for k =1:length(all_trib_ids_all_basins{i}); % for the number of knickpoints found in that stream
        stream_id_matrix(k) = i ;  % we want to record the current stream too
    end
    stream_id_array{i} = transpose(stream_id_matrix);  % save the stream id matrix
end

all_trib_ids = [];
stream_id = [];
all_beta_list = [];
% now we need to put them together for each basin

for i = 1:number_of_basins2 % % for all stream basins
    all_trib_ids = vertcat(all_trib_ids,all_trib_ids_all_basins{i});  % vertcally concatinate all the trib ids from each basin stored in the cell array above
    stream_id = vertcat(stream_id,stream_id_array{i});
    all_beta_list = vertcat(all_beta_list,all_beta_all_basins{i});
end

%
total_number_of_kps = length(all_kp_elevations);
% Total  knickpoint number
kp_id_number = 1:1:total_number_of_kps;
kp_id_number = transpose(kp_id_number);    % kp id # (1)

all_trib_ids;   % this goes in the database  tributary id # (3)
stream_id;     % this too  stream id # (2)

all_beta_list; % beta used for that tributary (slope of blf) (4)
% ^^ this is slightly different from beta because beta is fixed to the mouth of the stream, this is simply the linear regression

% parameters used  (last 4 columns in database)
sgolay_window_size = zeros(total_number_of_kps,1);
max_kp_search_dist = zeros(total_number_of_kps,1);
min_kp_size_pre_combination = zeros(total_number_of_kps,1);
min_kp_size_post_combination = zeros(total_number_of_kps,1);
min_kp_slope_detrended = zeros(total_number_of_kps,1);
min_trib_size_cells = zeros(total_number_of_kps,1);

for k = 1 : total_number_of_kps
    sgolay_window_size(k) = KZP_parameters.smoothing_window;   % (17)
    max_kp_search_dist(k) = KZP_parameters.lumping_search_distance ;   % (18)
    min_kp_size_pre_combination(k)= KZP_parameters.min_kp_size1 ;    % (19)
    min_kp_size_post_combination(k)= KZP_parameters.min_kp_size2;     % (20)
    min_kp_slope_detrended (k) = KZP_parameters.min_kp_slope ;   % (21)
    min_trib_size_cells(k) = KZP_parameters.min_trib_size; % (22)
end
%
% CALCULATE REAL SLOPE (Y/X in landscape)
% Slope raw
dy = all_kp_elevations - all_kp_elevations_bases; % (m elevation drop)
dx = all_kp_distance - all_kp_distance_bases; % (m distance upstream)
kp_slope_raw = dy./dx; %13
% slope normalized for drainage area
kp_slope_detrend_elv_dist = all_kp_magnitude./dx; % 14

% list of database measurments and parameters
% 1: knickpoint #
% 2: stream id #
% 3: tributary id #
% 4: slope of stream (bfl slope) elev/chi
% 5: chi
% 6: elevation (m)
% 7: knickpoint magnitude (detrended elevation drop (m))
% 8: knickpoint relief (m)
% 9: easting meters utm
% 10: northing meters utm
% 11: upstream drainage area m2
% 12: knickpoint distance upstream (m)
% 13: knickpoint slope (elev drop/distance upstream)
% 14: elevation drop (detrended elevation drop/distance upstream)
% 15: knickpoint slope (detrended elevation drop/chi)
% 16: knickpoint length (m)
% 17: sgolay smoothing window size
% 18: knickpoint lumping search window size
% 19: minimu knickkpoint size pre-lumping
% 20: minimum knickpoint size post-lumping (final minimum knickpoitn size)
% 21: minimum steepness anomoly
% 22: minimum stream size for analysis (cells)

% make the database

% knickpoint lips
Knickpoint_Database_Lips = horzcat (kp_id_number,stream_id,all_trib_ids, ...
    all_beta_list, all_kp_chi,all_kp_elevations,all_kp_magnitude, all_kp_relief, ...
    all_kp_easting, all_kp_northing ,all_kp_upstream_DA ,all_kp_distance, ...
    kp_slope_raw, kp_slope_detrend_elv_dist, all_kp_slope_detrend, ...
    all_kp_length_distance, sgolay_window_size, max_kp_search_dist, ...
    min_kp_size_pre_combination, min_kp_size_post_combination, ...
    min_kp_slope_detrended,min_trib_size_cells);
% how do you label each column in the database?

%knickpoint bases (database containing all information about the bottom of
%the knickpoints)
Knickpoint_Database_Bases = horzcat (kp_id_number,stream_id,all_trib_ids,...
    all_beta_list, all_kp_chi_bases,all_kp_elevations_bases, ...
    all_kp_magnitude, all_kp_relief, all_kp_easting_bases, all_kp_northing_bases, ...
    all_kp_upstream_DA_bases, all_kp_distance_bases, kp_slope_raw, ...
    kp_slope_detrend_elv_dist, all_kp_slope_detrend, ...
    all_kp_length_distance, sgolay_window_size, max_kp_search_dist, ...
    min_kp_size_pre_combination, min_kp_size_post_combination, ...
    min_kp_slope_detrended, min_trib_size_cells);

% write database to csv file with header
hdr = {'1kp_id', '2stream_id', '3trib_id', '4sl_str', '5chi', '6elev_m', ...
    '7kp_magnt', '8kp_Rel','9Easting_m', '10North_m', '11DA_m2', '12kp_DFM', ...
    '13kp_slp', '14kp_slp_dt', '15slp_dt/chi', '16kp_len_m', '17sgol_smv', ...
    '18kp_lm_ws', '19kp_pr_lu', '20kp_pt_lu', '21kp_stp_a', ...
    '22strea_sz'};
txt = sprintf('%s, ',hdr{:});
txt(end)='';
kp_lips_fn = strcat(KZP_parameters.DEM_basename_nodir, '_kp_lips.csv');
dlmwrite(kp_lips_fn,txt,'');
dlmwrite(kp_lips_fn,Knickpoint_Database_Lips,'-append','delimiter',',', 'precision', '%.5f');

kp_bases_fn = strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases.csv');
dlmwrite(kp_bases_fn,txt,'');
dlmwrite(kp_bases_fn,Knickpoint_Database_Bases,'-append','delimiter',',', 'precision', '%.5f');

%verify if projection file exists
if exist('projection.prj', 'file') ~= 2
    eval([KZP_parameters.gdalsrsinfo_cmd, ' -o wkt ', KZP_parameters.DEM_fname, '> projection.prj']);
end
% Create .vrt file to convert .csv to shapefile or other vector-file
% formats with gdal
kp_lips_shapeout_fn = strcat(KZP_parameters.DEM_basename_nodir, '_kp_lips.shp');
kp_lips_shapeout_fn_out = strcat(KZP_parameters.DEM_basename_nodir, '_kp_lips.out');
kp_lips_shapeout_fn_all = strcat(KZP_parameters.DEM_basename_nodir, '_kp_lips.*');
kp_lips_crt_fn = strcat(KZP_parameters.DEM_basename_nodir, '_kp_lips.crt');
lips_FID = fopen(kp_lips_crt_fn, 'w+');
string2write = sprintf(['<OGRVRTDataSource>\n  <OGRVRTLayer name=\"%s\">\n', ...
    '    <SrcDataSource relativeToVRT=\"1\">%s</SrcDataSource>\n', ...
    '    <GeometryType>wkbPoint</GeometryType>\n', ...
    '    <LayerSRS>WGS84</LayerSRS>\n',...
    '    <GeometryField encoding="PointFromColumns" x="9Easting_m" y="10North_m"/>\n',...
    '    <Field name="1kp_id" type="Integer" width="8"/>\n', ...
    '    <Field name="2stream_id" type="Integer" width="8"/>\n', ...
    '    <Field name="3trib_id" type="Integer" width="8"/>\n', ...
    '    <Field name="4sl_str" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="5chi" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="6elev_m" type="Integer" width="8"/>\n', ...
    '    <Field name="7kp_magnt" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="8kp_Rel" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="9Easting_m" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="10North_m" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="11DA_m2" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="12kp_DFM" type="Integer" width="8"/>\n', ...
    '    <Field name="13kp_slp" type="Real" width="8" precision="7"/>\n',  ...
    '    <Field name="14kp_slp_dt" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="15slp_dt/chi" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="16kp_len_m" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="17sgol_smv" type="Integer" width="8"/>\n', ...
    '    <Field name="18kp_lm_ws" type="Real" width="8" precision="7"/>\n',  ...
    '    <Field name="19kp_pr_lu" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="20kp_pt_lu" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="21kp_stp_a" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="22strea_sz" type="Real" width="8" precision="7"/>\n', ...
    '  </OGRVRTLayer>\n', '</OGRVRTDataSource>\n'], strcat(KZP_parameters.DEM_basename_nodir, '_kp_lips'), ...
    kp_lips_fn);
fwrite(lips_FID, string2write);
fclose(lips_FID);
if exist(strcat(sprintf('%s%s', KZP_parameters.KZP_shapefile_dirname, KZP_parameters.dir_sep), kp_lips_shapeout_fn), 'file') ~= 2
    eval([KZP_parameters.gdalsrsinfo_cmd, ' -o wkt ', KZP_parameters.DEM_fname, '> projection.prj']);
    eval([KZP_parameters.ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj -f "ESRI Shapefile" ', ...
        kp_lips_shapeout_fn, ' ', kp_lips_crt_fn, ' 2> ', kp_lips_shapeout_fn_out]);
end
if exist(kp_lips_shapeout_fn_all, 'file') ~= 2
    eval([KZP_parameters.mv_cmd, ' ', kp_lips_shapeout_fn_all, ' ', sprintf('%s%s', KZP_parameters.KZP_shapefile_dirname, KZP_parameters.dir_sep)]);
end

kp_bases_shapeout_fn = strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases.shp');
kp_bases_shapeout_fn_out = strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases.out');
kp_bases_shapeout_fn_all = strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases.*');
kp_bases_crt_fn = strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases.crt');
bases_FID = fopen(kp_bases_crt_fn, 'w+');
string2write = sprintf(['<OGRVRTDataSource>\n  <OGRVRTLayer name=\"%s\">\n', ...
    '    <SrcDataSource relativeToVRT=\"1\">%s</SrcDataSource>\n', ...
    '    <GeometryType>wkbPoint</GeometryType>\n', ...
    '    <LayerSRS>WGS84</LayerSRS>\n',...
    '    <GeometryField encoding="PointFromColumns" x="9Easting_m" y="10North_m"/>\n',...
    '    <Field name="1kp_id" type="Integer" width="8"/>\n', ...
    '    <Field name="2stream_id" type="Integer" width="8"/>\n', ...
    '    <Field name="3trib_id" type="Integer" width="8"/>\n', ...
    '    <Field name="4sl_str" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="5chi" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="6elev_m" type="Integer" width="8"/>\n', ...
    '    <Field name="7kp_magnt" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="8kp_Rel" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="9Easting_m" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="10North_m" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="11DA_m2" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="12kp_DFM" type="Integer" width="8"/>\n', ...
    '    <Field name="13kp_slp" type="Real" width="8" precision="7"/>\n',  ...
    '    <Field name="14kp_slp_dt" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="15slp_dt/chi" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="16kp_len_m" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="17sgol_smv" type="Integer" width="8"/>\n', ...
    '    <Field name="18kp_lm_ws" type="Real" width="8" precision="7"/>\n',  ...
    '    <Field name="19kp_pr_lu" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="20kp_pt_lu" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="21kp_stp_a" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="22strea_sz" type="Real" width="8" precision="7"/>\n', ...
    '  </OGRVRTLayer>\n', '</OGRVRTDataSource>\n'], ...
    strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases'), kp_bases_fn);
fwrite(bases_FID, string2write);
fclose(bases_FID);
if exist(strcat(sprintf('%s%s', KZP_parameters.KZP_shapefile_dirname, KZP_parameters.dir_sep), kp_bases_shapeout_fn), 'file') ~= 2
    eval([KZP_parameters.gdalsrsinfo_cmd, ' -o wkt ', KZP_parameters.DEM_fname, '> projection.prj']);
    eval([KZP_parameters.ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj -f "ESRI Shapefile" ', ...
        kp_bases_shapeout_fn, ' ', kp_bases_crt_fn, ' 2> ', kp_bases_shapeout_fn_out]);
end
if exist(kp_bases_shapeout_fn_all, 'file') ~= 2
    eval([KZP_parameters.mv_cmd, ' ', kp_bases_shapeout_fn_all, ' ', sprintf('%s%s', KZP_parameters.KZP_shapefile_dirname, KZP_parameters.dir_sep)]);
end
