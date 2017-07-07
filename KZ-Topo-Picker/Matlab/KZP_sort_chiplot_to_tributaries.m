function [elev_stored_trimmed, chi_stored_trimmed, eastingUTM11_stored_trimmed,...
    northingUTM11_stored_trimmed,upstream_DA_stored_trimmed,distance_stored_trimmed] = ...
    KZP_sort_chiplot_to_tributaries(current_basin_chiplot,Min_trib_size)

% organize AOI_chiplot structure array so that cell arrays contain tributary attributes
% elevation, chi, N, E, distance, and DA information (for each tributary)
% Extract chi, elev, dist, N,E, DA information for each tributary within the current basin
idx_confluences = ~isnan(current_basin_chiplot.chi);
% sets cells that are NaN's to 0   (the confluences are where cells within .chi = NaN)

confluences = find(idx_confluences == 0) ;
% FINDS NAN VALUES (CONFLUENCES)

confluences = [1;confluences];
% GETS CELL LOCATIONS OF CONFLUENCES
% Need to add the first cell in the structure array to the NaN's (start of trunkstream)

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

% Below are needed if looping through multiple basins where
% there are different number of tributaries.  The cell arrays
% must be emptied because if not, basins with more tributaries
% will not have all of their cells rewritten in the the next
% loop when looking at a basin with fewer tributaries
chi_stored = {}; % erase data from previous loop
elev_stored = {}; % erase data from previous loop
eastingUTM11_stored = {}; % erase data from previous loop
northingUTM11_stored={} ;% erase data from previous loop
upstream_DA_stored={} ;% erase data from previous loop
distance_stored = {}; % erase data from previous loop

% ^^ above cell arrays  will store information for each trib.
% within a basin (allows you to plot and construct long
% profiles for each individual tributary within a basin)

%  This starts from the headwaters of the trunk stream (cell 1), measures to the next NaN value (confluence)
for j = 1:(length(confluences)-1);
    tributaries{j} = [confluences(j):1:confluences(j+1)] ;
    % This stores the cells between the each confluence
    %(from the headwaters to the mouth of the tributary). The trunk stream
    % is included. The values are stored in tributaries{j}. Each cell in
    % tributaries{j} contains information for 1 tributary
    % (including ts, which is just the first cell in "tributaries")
    
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

% now, 'chi_stored' contains the chi information in a different
% cell for each tributary 'j' in the basin

% get rid of tributaries that are too small (less than
% 'Min_trib_size' in parameters script)

% find the tributaries that are too small and mark them so we
% can remove them. This avoids analyzing tributaries that are
% smaller than the smoothing window size or tributaries that
% only contain a few datapoints.

for j = 1:length(tributaries);
    % Mark cells for short tributaries in the stored elevation cell array (same
    % cells as chi, .x, .y, ect..
    if length(chi_stored{j}) < Min_trib_size ;
        % measure the length of the chi vector for each tributary.
        % if that vector is too small we ignore that stream segment
        chi_stored{j} = -9999 ;
        elev_stored{j} = -9999;
        eastingUTM11_stored{j} = -9999; % mark as -9999
        northingUTM11_stored{j} = -9999;
        upstream_DA_stored{j} = -9999;
        distance_stored{j} = -9999;
    end
end
% now index out tributaries that are too small
for j = 1:length(tributaries);
    if chi_stored{j} == -9999;
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