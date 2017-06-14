function DEM_FAC_pt_export(shapefile_dirname, AOI_DEM, AOI_FAC, AOI_FAC_pt_threshold, AOI_FD, ...
    AOI_mg, dir_sep, DEM_basename, gdalsrsinfo_cmd, ogr2ogr_cmd, remove_cmd, mv_cmd)
%
% write point shapefile that contains drainage area, dem elevation, stream
% order, river slope. All points are written with drainage area > AOI_FAC_pt_threshold
%

AOI_FAC_m2 = AOI_FAC .* (AOI_FAC.cellsize^2);

idx_nan = find(AOI_FAC_m2.Z < AOI_FAC_pt_threshold);
AOI_FAC_m2.Z(idx_nan) = NaN;
clear idx_nan
AOI_STO = streamorder(AOI_FD, ~isnan(AOI_FAC_m2));
AOI_rivers_pt = STREAMobj(AOI_FD, ~isnan(AOI_FAC_m2));

AOI_rivers_slope = gradient(AOI_rivers_pt,AOI_mg,'unit','tangent');

idx_fac = find(~isnan(AOI_FAC_m2.Z));
[x,y] = ind2coord(AOI_FAC_m2, idx_fac);
AOI_FAC_m2_pt = AOI_FAC_m2.Z(idx_fac);
AOI_DEM_pt = AOI_DEM.Z(idx_fac);
AOI_STO_pt = AOI_STO.Z(idx_fac);

p = mappoint(x, y, 'DA_m2', double(AOI_FAC_m2_pt), 'DEM', double(AOI_DEM_pt), ...
    'STROrder', double(AOI_STO_pt), 'slp', double(AOI_rivers_slope));

shapeout_fn = sprintf('%s%s%s_PT_%s.shp', shapefile_dirname, dir_sep, ...
    DEM_basename, num2str(AOI_FAC_pt_threshold,'%1.0G'));
shapeout_fn_out = sprintf('%s%s%s_PT_%s.out', shapefile_dirname, dir_sep, ...
    DEM_basename, num2str(AOI_FAC_pt_threshold,'%1.0G'));
shapeout_all_fn = sprintf('%s%s%s_PT_%s.*', shapefile_dirname, dir_sep, ...
    DEM_basename, num2str(AOI_FAC_pt_threshold,'%1.0G'));
k = strfind(shapeout_fn, '+'); shapeout_fn(k) = [];
shapeout_all_fn(k) = []; shapeout_fn_out(k) = [];
shapeout_fn_prj = sprintf('%s%s%s_%s_PT_proj.shp', ...
    shapefile_dirname, dir_sep, DEM_basename, num2str(AOI_FAC_pt_threshold,'%1.0e'));
if exist(shapeout_fn_prj, 'file') ~= 2
    fprintf(1,'\twriting shapefile: %s\n', shapeout_fn);
    if ~isempty(p)
        shapewrite(p, shapeout_fn);
        %Because shapewrite doesn't add projection information, we have to add
        %these manually via ogr2ogr or ArcMAP (or something similar)
        eval([ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj ', ...
            shapeout_fn_prj, ' ', shapeout_fn, ' > ', shapeout_fn_out]);
        eval([remove_cmd, ' ', shapeout_all_fn]);
    end
end
