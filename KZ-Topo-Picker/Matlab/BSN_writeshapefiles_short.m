function BSN_writeshapefiles(dir_sep, DEM_basename, DEM_basename_nodir, DEM_fname, AOI_STR_MS, ...
    AOI_STR_MS_area1, AOI_STR_MS_area2, AOI_STR_MS_areaonly, ...
    str_area1, str_area2, gdalsrsinfo_cmd, ogr2ogr_cmd, remove_cmd, mv_cmd)

shapeout_fn = sprintf('shapefiles%s%s_all_MS.shp', dir_sep, DEM_basename);
shapeout_fn_ogr = sprintf('shapefiles%s%s_all_MS.out', dir_sep, DEM_basename);
shapeout_all_fn = sprintf('shapefiles%s%s_all_MS.*', dir_sep, DEM_basename);
shapeout_fn_prj = sprintf('shapefiles%s%s_all_MS_proj.shp', dir_sep, DEM_basename);
if exist(shapeout_fn_prj, 'file') ~= 2
    fprintf(1,'\twriting shapefile: %s\n', shapeout_fn);
    shapewrite(AOI_STR_MS,shapeout_fn);
    %Because shapewrite doesn't add projection information, we have to add
    %these manually via ogr2ogr or ArcMAP (or something similar)
    eval([gdalsrsinfo_cmd, ' -o wkt ', DEM_fname, '> projection.prj']);
    eval([ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj ', ...
        shapeout_fn_prj, ' ', shapeout_fn, ' > ', shapeout_fn_ogr]);
    eval([remove_cmd, ' ', shapeout_all_fn]);
end

shapeout_fn = sprintf('shapefiles%s%s_MS_%s.shp', dir_sep, ...
    DEM_basename, num2str(str_area1,'%1.0G'));
shapeout_fn_out = sprintf('shapefiles%s%s_MS_%s.out', dir_sep, ...
    DEM_basename, num2str(str_area1,'%1.0G'));
shapeout_all_fn = sprintf('shapefiles%s%s_MS_%s.*', dir_sep, ...
    DEM_basename, num2str(str_area1,'%1.0G'));
k = strfind(shapeout_fn, '+'); shapeout_fn(k) = [];
shapeout_all_fn(k) = []; shapeout_fn_out(k) = [];
shapeout_fn_prj = sprintf('shapefiles%s%s_%s_MS_proj.shp', ...
    dir_sep, DEM_basename, num2str(str_area1,'%1.0e'));
if exist(shapeout_fn_prj, 'file') ~= 2
    fprintf(1,'\twriting shapefile: %s\n', shapeout_fn);
    if isstruct(AOI_STR_MS_area1)
        if ~isempty(AOI_STR_MS_area1)
            shapewrite(AOI_STR_MS_area1, shapeout_fn);
            %Because shapewrite doesn't add projection information, we have to add
            %these manually via ogr2ogr or ArcMAP (or something similar)
            eval([gdalsrsinfo_cmd, ' -o wkt ', DEM_fname, '> projection.prj']);
            eval([ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj ', ...
                shapeout_fn_prj, ' ', shapeout_fn, ' > ', shapeout_fn_out]);
            eval([remove_cmd, ' ', shapeout_all_fn]);
        end
    end
end

shapeout_fn = sprintf('shapefiles%s%s_MS_%s_areaonly.shp', dir_sep, ...
    DEM_basename, num2str(str_area2,'%1.0G'));
shapeout_fn_out = sprintf('shapefiles%s%s_MS_%s_areaonly.out', dir_sep, ...
    DEM_basename, num2str(str_area2,'%1.0G'));
shapeout_all_fn = sprintf('shapefiles%s%s_MS_%s_areaonly.*', dir_sep, ...
    DEM_basename, num2str(str_area2,'%1.0G'));
k = strfind(shapeout_fn, '+'); shapeout_fn(k) = [];
shapeout_all_fn(k) = []; shapeout_fn_out(k) = [];
shapeout_fn_prj = sprintf('shapefiles%s%s_%s_MS_areaonly_proj.shp', ...
    dir_sep, DEM_basename, num2str(str_area1,'%1.0e'));
if exist(shapeout_fn_prj, 'file') ~= 2
    fprintf(1,'\twriting shapefile: %s\n', shapeout_fn);
    if isstruct(AOI_STR_MS_areaonly)
        if ~isempty(AOI_STR_MS_areaonly)
            shapewrite(AOI_STR_MS_areaonly, shapeout_fn);
            %Because shapewrite doesn't add projection information, we have to add
            %these manually via ogr2ogr or ArcMAP (or something similar)
            eval([gdalsrsinfo_cmd, ' -o wkt ', DEM_fname, '> projection.prj']);
            eval([ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj ', ...
                shapeout_fn_prj, ' ', shapeout_fn, ' > ', shapeout_fn_out]);
            eval([remove_cmd, ' ', shapeout_all_fn]);
        end
    end
end

shapeout_fn = sprintf('shapefiles%s%s_MS_%s.shp', dir_sep, ...
    DEM_basename, num2str(str_area2,'%1.0G'));
shapeout_fn_out = sprintf('shapefiles%s%s_MS_%s.out', dir_sep, ...
    DEM_basename, num2str(str_area2,'%1.0G'));
shapeout_all_fn = sprintf('shapefiles%s%s_MS_%s.*', dir_sep, ...
    DEM_basename, num2str(str_area2,'%1.0G'));
k = strfind(shapeout_fn, '+'); shapeout_fn(k) = [];
shapeout_all_fn(k) = []; shapeout_fn_out(k) = [];
shapeout_fn_prj = sprintf('shapefiles%s%s_%s_MS_proj.shp', dir_sep, ...
    DEM_basename, num2str(str_area2,'%1.0e'));
if exist(shapeout_fn_prj, 'file') ~= 2
    fprintf(1,'\twriting shapefile: %s\n', shapeout_fn);
    if isstruct(AOI_STR_MS_area2)
        if ~isempty(AOI_STR_MS_area2)
            shapewrite(AOI_STR_MS_area2, shapeout_fn);
            %Because shapewrite doesn't add projection information, we have to add
            %these manually via ogr2ogr or ArcMAP (or something similar)
            eval([gdalsrsinfo_cmd, ' -o wkt ', DEM_fname, '> projection.prj']);
            eval([ogr2ogr_cmd, ' -s_srs projection.prj -t_srs projection.prj ', ...
                shapeout_fn_prj, ' ', shapeout_fn, ' >', shapeout_fn_out]);
            eval([remove_cmd, ' ', shapeout_all_fn]);
        end
    end
end
