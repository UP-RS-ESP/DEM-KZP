function [AOI_ridgecrest_MS] = ...
    DEM_ridgecrest(AOI_dbasins_regionprops, AOI_dbasins, AOI_dbasins_outlet)
%extract ridgecrests from all drainage basins
%AOI_ridgecrest contains ridgecrests from all drainage basins in:
%[BasinNr imagesize]
%

fprintf(1,'\textract ridge crests for all basins\n');

nob = length(AOI_dbasins_regionprops);
bsnnr_unique = unique(AOI_dbasins.Z); bsnnr_unique(1) = [];

for i = 1:nob
    AOI_dbasins_cp = AOI_dbasins;
    idx = find(AOI_dbasins.Z == bsnnr_unique(i));
    AOI_dbasins_cp.Z = logical(zeros(size(AOI_dbasins.Z), 'int8'));
    AOI_dbasins_cp.Z(idx) = 1;
    AOI_dbasins_outlet_cp = AOI_dbasins_outlet(i);
    [AOI_dbasins_outlet_cp_x, AOI_dbasins_outlet_cp_y] = ind2sub(size(AOI_dbasins_cp.Z), AOI_dbasins_outlet_cp);
    try
        AOI_dbasins_trace = bwtraceboundary(AOI_dbasins_cp.Z, [AOI_dbasins_outlet_cp_x AOI_dbasins_outlet_cp_y], 'W');
    catch error
        try
            AOI_dbasins_trace = bwtraceboundary(AOI_dbasins_cp.Z, [AOI_dbasins_outlet_cp_x AOI_dbasins_outlet_cp_y], 'N');
        catch error
            try
                AOI_dbasins_trace = bwtraceboundary(AOI_dbasins_cp.Z, [AOI_dbasins_outlet_cp_x AOI_dbasins_outlet_cp_y], 'S');
            catch error
                AOI_dbasins_trace = bwtraceboundary(AOI_dbasins_cp.Z, [AOI_dbasins_outlet_cp_x AOI_dbasins_outlet_cp_y], 'E');
            end
        end
    end
    if isempty(AOI_dbasins_trace)
      fprintf(1,'  AOI_dbasins_trace is empty for basin id: %d (iteration %d)\n', bsnnr_unique(i), i);
      continue
    end
    [x,y] = sub2coord(AOI_dbasins_cp,AOI_dbasins_trace(:,1),AOI_dbasins_trace(:,2));
    [IX] = coord2ind(AOI_dbasins, x, y);
    m = [x y];
    %convert x y to pixel indicators
    distance_xy = sqrt(sum((diff(m,1,1).^2), 2));
    distance_xy(1:length(distance_xy)+1) = [0; distance_xy];
    distance_xy_cumsum = cumsum(distance_xy);
    AOI_ridgecrest_MS(i).Geometry='Point';
    AOI_ridgecrest_MS(i).x = x;
    AOI_ridgecrest_MS(i).y = y;
    AOI_ridgecrest_MS(i).area = AOI_dbasins_regionprops(i).Area;
    AOI_ridgecrest_MS(i).pindex = IX;
    AOI_ridgecrest_MS(i).distance = distance_xy_cumsum;
    %AOI_ridgecrest_MS(i).distance_cum = distance_xy_cumsum;
    AOI_ridgecrest_MS(i).distance_norm = distance_xy_cumsum./max(distance_xy_cumsum);
    AOI_ridgecrest_MS(i).ID = bsnnr_unique(i);
end
