function [AOI_ridgecrest_values] = ...
    DEM_ridgevalues(AOI_ridgecrest_MS, AOI_values)
% extract ridge values from structure AOI_ridgecrest_MS
%
nob = size(AOI_ridgecrest_MS,2);
AOI_ridgecrest_values = struct;

for i = 1:nob
    index = AOI_ridgecrest_MS(i).pindex;
    AOI_ridgecrest_values(i).values = AOI_values.Z(index);
    AOI_ridgecrest_values(i).values_norm = (AOI_ridgecrest_values(i).values-min(AOI_ridgecrest_values(i).values)) ./ (max(AOI_ridgecrest_values(i).values) - min(AOI_ridgecrest_values(i).values));
end
