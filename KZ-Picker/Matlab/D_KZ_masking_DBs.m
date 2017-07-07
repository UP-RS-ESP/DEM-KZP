%% Use selections from data cursor to mask draiange basins of interest

% data cursor selections include the index of each drainage basin the user
% wants to analyze.  

% We will mask the full DEM, StreamOBJ, and FlowAccum to calculate
% slope-area and chi plots for each individual basin

number_of_basins = length(basin_index);
% record # of basins from the selections made using the data cursor

AOI_resolution = AOI_DEM.refmat(2,1);
% record the resolution of the DEM (used to covert minimum drainage area to
% pixels)

% preallocate cell array to store generated streamOBJ
AOI_STR_largest = cell(1,number_of_basins);
AOI_chiplot = cell(1,number_of_basins);
% index values for each basin


for i = 1:number_of_basins
    current_basin_index = basin_index(i);
    basin_num_str = num2str(i);
    fprintf(1,'step 3 of 4: Masking Drainage Basins (may take a few minutes if a large DEM) BASIN#: %s\n', basin_num_str)
    
    AOI_dbasins_index = AOI_dbasins;
    idx_current_basin = find(AOI_dbasins_index.Z ~= current_basin_index); % index out your drainge basin
    AOI_dbasins_index.Z(idx_current_basin) = NaN;
    
    % indexes DEM to include just the basin of interest
    AOI_mg_current = AOI_mg;  % makes a copy of the minimum gradient dem
    AOI_mg_current.Z(idx_current_basin) = NaN;
    
    % index FIL DEM to include just basin of interest
    AOI_FIL_current = AOI_FIL; % makes a copy of fil
    AOI_FIL_current.Z(idx_current_basin) = NaN; % indexes out current fil for just basin
        
    % need to recalculate flowDIR for only the basin of interest
    % flowDIR has a different method of storing values than other rasters
    % so you need to recalculate it but only for the basins of interest
    AOI_FD_current = FLOWobj(AOI_FIL_current,'preprocess','carve');  
    
    % generate a flowaccumulation grid for basin of interest
    AOI_FAC_current = flowacc(AOI_FD_current); 
    
    minApix = Min_DA_threshold/(AOI_resolution*AOI_resolution);
    minApix = ceil(minApix);  % convert area threshold to meters specified above
    
    AOI_FAC_current_w = AOI_FAC_current > minApix;  
    % masks flow accum grid for each basin (above minimum drainage area threshold)
    
    % create streamOBJ for basin of interest
    AOI_STR_current_w = STREAMobj(AOI_FD_current,AOI_FAC_current_w);
    
    % make sure streamOBJ doesn't have two outlets (would cause chiplot to
    % crash), so pick largest connected component
    AOI_STR_largest{i} = klargestconncomps(AOI_STR_current_w);
    
    
    %% Generate Slope-Area Plot and Chi-Plot Data
    
   
    %% Generate Chi plot
if theta_bf_option == 1 % if selected not to use fixed ref concavity 
    % (allow code to find best fit m/n)
    S_chiplot = ...
    chiplot(AOI_STR_largest{i},...
    AOI_mg_current,AOI_FAC_current, ...
    'a0', 1, 'fitto', 'all', 'plot', false); % use concavity that best collapses tributaries
else
    S_chiplot = ...
    chiplot(AOI_STR_largest{i},...
    AOI_mg_current,AOI_FAC_current, ...
    'a0', 1,'mn',-theta_ref, 'plot', false); % use reference concavity 
end
    AOI_chiplot{i} = S_chiplot; % saves chiplot structure array for each basin
    
    % MUST STORE THE VALUES BELOW FOR EACH BASIN!!!!
    theta_bf(i) = S_chiplot.mn; % saves m/n ratio used for each basin
    beta_bf(i) = S_chiplot.beta; % saves regression value for each basin
    ks_chiplot(i) = S_chiplot.ks; % saves ks calculated for each basin
    % use a reference drainage area of 1 (slope = ks) 
    % this value changes the scaling of your chi axis 
    % but doesn't affect the steepness or concavity of the profile
    
    % calculate slope area relation and
     % generate slope area plot 
    h = figure;set(h, 'Visible', 'off')  % create a slope/area figure for each stream)
    AOI_slopearea_current = slopearea(AOI_STR_largest{i},...
    AOI_mg_current,AOI_FAC_current, ... % define trunk stream and DEM
    'areabins', 100, ... % nr of log-bins (100)
    'areabinlocs', 'median', ... % log-area bins use median area of each bin
    'gradaggfun', 'median', ... % log-slope bins use median slope of each bin
    'fitmethod', 'ls', ... % least squares fit method
    'streamgradient', 'robust', ... % robust-fitting method for gradient
    'theta', -theta_ref, ... % fix slope area to reference concavity
    'plot', true);
    grid on
    
    
    % create string for ks from slopearea plot
    % MUST STORE THE VALUE BELOW FOR EACH BASIN!!!!
    ksn_SA_store(i) = AOI_slopearea_current.ks; 
    % save ks for slope area calculation (will be ksn because uses fixed ref. concavity)
    ks_string = [' ksn SA = ', num2str(AOI_slopearea_current.ks)];
    text(Min_DA_threshold,min(AOI_slopearea_current.g),ks_string,'HorizontalAlignment','right');
    % convert to string and plot on figure outputted from slopearea
    % function
    
    % create string for concavity from slopearea plot
    theta_string = [' theta SA = ', num2str(AOI_slopearea_current.theta)];
    text(Min_DA_threshold,min(AOI_slopearea_current.g),theta_string,'HorizontalAlignment','left');
    % save reference concavity and add to figure outputted from slopearea
    % function
    
    % create string for ks from chiplot
    ks_string_chiplot = [' ks chiplot = ', num2str(S_chiplot.ks)];
    text(max(AOI_slopearea_current.a),max(AOI_slopearea_current.g),ks_string_chiplot,'HorizontalAlignment','right');
    % save ks from chiplot (ksn if fixed concavity was used 'specified in
    % parameters script') and add to fig outputted from slopearea function
    
    % create string for concavity from chiplot
    theta_string_chiplot = [' m/n chiplot = ', num2str(S_chiplot.mn)];
    text(max(AOI_slopearea_current.a),max(AOI_slopearea_current.g),theta_string_chiplot,'HorizontalAlignment','left');
    % save m/n from chiplot (either fixed or best fit m/n depending on what
    % was specified in parameters script) and add to fig outputted from slopearea function
     
    %Generate name for slope-area figure
    name_string = ['SA_plot_basin_', num2str(basin_index(i))];
    % need to add path to directory where we store slope-area figures
    
    saveas(h,[pwd '/SA_plots/',name_string,'.png']);   
    % save figure in directory
    clf
end
    
    