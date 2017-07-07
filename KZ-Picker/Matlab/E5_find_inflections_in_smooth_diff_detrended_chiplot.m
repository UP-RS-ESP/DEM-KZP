function [kps_cells_current_trib, base_cells_current_trib, kps_lips_length, kps_base_length,stream_mouth_cell] = ...
    E5_find_inflections_in_smooth_diff_detrended_chiplot(diff_smooth_detrended_elev_current)

local_max=[]; % create vectors to store found local minima and local maxima 
local_min=[];

% find inflections (zeros in the differentiated, smoothed detrended curve)
stream_mouth_cell = length(diff_smooth_detrended_elev_current);
% cell at the mouth of the stream

for k = 1:stream_mouth_cell-1
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

% find the cells within the current tributary that are potential knickpoints
kps_cells_current_trib = find(local_max== 1)'; % cells containing knickpoint lips
base_cells_current_trib = find(local_min== 1)'; % cells containing knickpoint bases

kps_lips_length = length(kps_cells_current_trib);  % find how many concave down inflections are found
kps_base_length = length(base_cells_current_trib); % find how many concave up inflections are found
