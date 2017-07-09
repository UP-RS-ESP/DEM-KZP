%save all knickpoints with additional information to a csv and shapefile
nr_of_kps = length(kz_num_stored);
nr_of_basins = length(kz_basin_id_kps);

for i = 1:nr_of_basins
    if ~isempty(kz_easting_all_tribs_lips{i})
        % for each basin
        
        % get current basin ID number
        current_basin_id = i;
        
        
        % save database files
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
        
        %EDIT: shortened version:
        % 1: knickpoint #
        % 2: easting coordinate of knickzone lip/base meters utm
        % 3: northing coordinate of knickzone lip/base meters utm
        % 4: upstream drainage area coordinate of knickzone lip/base m2
        % 5: knickzone distance upstream coordinate of knickzone lip/base (m)
        % 6: chi coordinate of knickzone lip/base
        % 7: elevation coordinate of knickzone lip/base (m)
        % 8: knickzone magnitude
        % 9: knickzone relief (m)
        % 10: knickzone length (distance usptream)
        % 11: knickzone slope (elev drop/distance upstream)
        
        kz_db_lips{i} = horzcat(kz_id_nr{i}, kz_easting_all_tribs_lips{i}, ...
            kz_northing_all_tribs_lips{i}, kz_up_da_all_tribs_lips{i}, ...
            kz_distance_all_tribs{i}, kz_chi_all_tribs_lips{i}, ...
            kz_elev_all_tribs_lips{i}, kz_magnitude_all_tribs{i}, ...
            kz_relief_all_tribs{i},  kz_length_all_tribs{i}, kz_face_slope_all_tribs{i}  );
        
        
        kz_db_bases{i} = horzcat(kz_id_nr{i}, kz_easting_all_tribs_bases{i}, ...
            kz_northing_all_tribs_bases{i}, kz_up_da_all_tribs_bases{i}, ...
            kz_distance_all_tribs{i}, kz_chi_all_tribs_bases{i}, ...
            kz_elev_all_tribs_lips{i}, kz_magnitude_all_tribs{i}, ...
            kz_relief_all_tribs{i}, kz_length_all_tribs{i}, kz_face_slope_all_tribs{i}  );
    end
end

kz_db_lips_all = {cat(1, kz_db_lips{:})};
kz_db_bases_all = {cat(1, kz_db_bases{:})};

%% write to csv:
hdr = {'1kp_id', '2kz_E', '3kz_N', '4kz_DA', '5kz_dist', '6kz_chi', ...
    '7kz_elev', '8kz_mag', '9kz_rel', '10kz_lgt', '11kz_slp'};
txt = sprintf('%s, ',hdr{:});
txt(end)='';
kp_lips_fn = strcat(KZP_parameters.DEM_basename_nodir, '_kp_lips.csv');
dlmwrite(kp_lips_fn,txt,'');
dlmwrite(kp_lips_fn,kz_db_lips_all,'-append','delimiter',',', 'precision', '%.5f');

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
    '    <GeometryField encoding="PointFromColumns" x="2kz_E" y="3kz_N"/>\n',...
    '    <Field name="1kp_id" type="Integer" width="8"/>\n', ...
    '    <Field name="2kz_E" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="3kz_N" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="4kz_DA" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="5kz_dist" type="Integer" width="8"/>\n', ...
    '    <Field name="6kz_chi" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="7kz_elev" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="8kz_mag" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="9kz_rel" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="10kz_lgt" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="11kz_slp" type="Real" width="8" precision="7"/>\n', ...
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

kp_bases_fn = strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases.csv');
dlmwrite(kp_bases_fn,txt,'');
dlmwrite(kp_bases_fn,kz_db_bases_all,'-append','delimiter',',', 'precision', '%.5f');

kp_bases_shapeout_fn = strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases.shp');
kp_bases_shapeout_fn_out = strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases.out');
kp_bases_shapeout_fn_all = strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases.*');
kp_bases_crt_fn = strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases.crt');
bases_FID = fopen(kp_bases_crt_fn, 'w+');
string2write = sprintf(['<OGRVRTDataSource>\n  <OGRVRTLayer name=\"%s\">\n', ...
    '    <SrcDataSource relativeToVRT=\"1\">%s</SrcDataSource>\n', ...
    '    <GeometryType>wkbPoint</GeometryType>\n', ...
    '    <LayerSRS>WGS84</LayerSRS>\n',...
    '    <GeometryField encoding="PointFromColumns" x="2kz_E" y="3kz_N"/>\n',...
    '    <Field name="1kp_id" type="Integer" width="8"/>\n', ...
    '    <Field name="2kz_E" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="3kz_N" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="4kz_DA" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="5kz_dist" type="Integer" width="8"/>\n', ...
    '    <Field name="6kz_chi" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="7kz_elev" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="8kz_mag" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="9kz_rel" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="10kz_lgt" type="Real" width="8" precision="7"/>\n', ...
    '    <Field name="11kz_slp" type="Real" width="8" precision="7"/>\n', ...
    '  </OGRVRTLayer>\n', '</OGRVRTDataSource>\n'], strcat(KZP_parameters.DEM_basename_nodir, '_kp_bases'), ...
    kp_bases_fn);
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
