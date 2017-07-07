%% comparison between algorithm outputs and calibration knickzones

% USE ONLY CALIBRATION LIPS! (no measure G.. calibrate only lips position)
% this script will run if specified in the parameters file

% load calibration knickzones located using profiler toolbar (Whipple et
% al., 2007)  

% load in each algorithm output with various parameter combination (saved
% .csv files from last script 'part_1')

% compare position of algorithm knickzones with the position of the
% calibration knickzones (located manually)

% compare measurements of algorithm knickzone height with knickpoint
% dimensions measured from calibration knickzones 

%% Use formula from Orlandini et al., 2011, also used by Clubb et al., 2014

% method relies on counting number of true positives (TP), false positives
% (FP), and false negatives (FN). 

% A true positive, would be a algorithm kz lies within an allowable error
% circle drawn around the calibration kz.

% A false positive, would be an algorithm kz. that does not lie within an
% allowable error circle drawn around the calibration kzs. 

% A false negative, would be an calibration kz. that does not have an
% algorithm knickzone within its allowable error circle. 

% reliability = #TP/(#TP + #FP)    Score of 1 is perfect

% sensitivity = #TP/(#TP + #FN)    Score of 1 is perfect

% WE need to weight these formulas by a measure of how accurate each
% method measures the calibration knickzone height

% If this is not done, the calibration proceedure will favor selecting more
% knickzones just to produce less false positives and negatives, but the 
% height of these knickzones will not match the height of the calibration 
% knickzones.

% FOR THIS REASON, we calculate the geometric error, which is the %
% difference in the algorithm knickzone height and the calibration
% knickzone height (only for TP knickzone lips)

% Geometric_error = 1 - abs(kp_relief_for_TP_alg_kps -
% kp_relief_for_closest_calib_kps)/(kp_relief_for_closest_calib_kps+kp_relief_for_TP_alg_kps)

% ideally, Geometric_error = 1 
% meaning: (relief of calibration kps = relief of algorithm kps)

% G = Average(geometric_error for all TP_alg_kps found with Param. comb)

% We can rate the algorithm's performance with the average of R, S, and G
% (this assumes = weighting, but could be adjusted)

% Max_R_S_G = ((#TP/(#TP + #FP)) + (#TP/(#TP + #FN)) + G)/3
               % Reliability       %sensitivity    %Geometric accuracy
               
% Find where these are maximized

% Uses error tolerance value from parameters script

% Spatial error tolerance should be something that is significantly less 
% than the length of average knickzones, and something that reflects the 
% uncertainty in knickzone position that comes from manual knickzone
% identification using the Whipple et al., 2007 method (profiler toolbar)

% Also, this value should be small enough that none of the circles drawn
% around each calibration point overlap. i.e. small enough so that one 
% algorithm knickzone boundary cannot create a true positive for more than one
% calibration knickzone boundary.  This may depend on your DEM resolution
% and the size of your landscape.


%% Step 1: Load calibration knickpoints

% calibration knickpoints are located using the profiler toolbar and
% procedures outlined in Wobus et al., 2006. I selected all of the
% knickpoint lips and combined these into a point-shapefile.  Then I
% selected all of the corresponding bases and combined these into a
% different point-shapefile.  This way, the base and lip information is
% stored in different shapefiles.

% make sure the FID for the knickpoint base is the same for its
% corresponding knickpoint lip (from ArcMap)


%% load attribute tables for calibration knickzone lips ONLY!
format long g

% set current directory as 'oldfolder'
oldfolder = pwd;

% load attribute table for the calibration knickzone lips. ROW 1
% is assumed to be the header! so starts reading at ROW 2!!!!!!
[Easting_calib_lips, Northing_calib_lips, numb_calib_kps] = ...
    F1_1_calib_load_calib_kz_lips_only(KZ_lips_calib_fname, KZ_calib_easting_column_num,...
    KZ_calib_northing_column_num);

% user needs calibration dataset to have 3 columns containing KZ northing,
% easting, and relief for each KZ boundary (lips only)



%% Reference calibration knickpoints to the streamobj file
% reference using geographic position (easting and northing)

% use function (will output cells in the chiplot structure array that
% contain the calibration knickzone lips and hte calibration knickzone
% bases
[calib_lips_cells, numb_calib_kz_ref] = ...
    F2_1_calib_calib_kz_reference2streamOBJ_only_calib_lips(current_basin_streamobj, Easting_calib_lips,...
    Northing_calib_lips, current_basin_chiplot,calibration_snapping_tolerance);

% ^^ function will plot calibration knickzones pre-refrencing to streamobj
% and then plot calibration knickzones after referencing them to the
% streamOBJ


% Get northing, easting, elevation, and distance upstream information for 
% each calibration knickpoint lip (on the chiplot network)

% northing
stream_loc_north = current_basin_chiplot.y; 
calib_lips_northing_str_loc = stream_loc_north(calib_lips_cells);

% easting
stream_loc_east = current_basin_chiplot.x;
calib_lips_easting_str_loc = stream_loc_east(calib_lips_cells); 

% elevation
stream_loc_elevation = current_basin_chiplot.elev;
calib_lips_elevation_str_loc = stream_loc_elevation(calib_lips_cells);

% distance from stream outlet
stream_loc_distance = current_basin_chiplot.distance;
calib_lips_distance_str_loc = stream_loc_distance(calib_lips_cells);

% plot longitudinal profile to double check to make sure KPS are properly
% referenced to streamOBJ
% position (after referencing them to the streamOBJ)

han1 = figure;set(han1, 'Visible', 'off');
hold on
plotdz(current_basin_streamobj,AOI_DEM,'color','k')
% create long-profile
% add calibration knickzone lips only
for i = 1:length(calib_lips_distance_str_loc)
    plot(calib_lips_distance_str_loc(i),calib_lips_elevation_str_loc (i),'r.','MarkerSize',14)
end

title('Calibration KZ Long-Profile Post-Referencing to StreamOBJ')

saveas(han1,[pwd '/calib_figure_outputs/','Calib_KZ_long_prof_post-ref.png']);   
% save figure in directory
clf 
% if this looks good, then you can start running the calibration

%% Load each set of algorithm kps and characterize fit with calibration kps

% FIRST, we have to perform this with the knickpoint lips

% we'll use the formula outlined in the beginning of the script to
% calculate reliability, sensitivity, and geometric accuracy

% change the directory to where you have stored all of the .csv files with
% the knickpoint lips positions found using the different parameter
% combinations

fprintf(1,'step 4 of 4: comparing results of different parameter combinations')

% find out # of different parameter combination outputs you made (from
% previous script)
Number_of_param_simulations = length(Knickpoint_Database_Lips_all_params);

% pre-allocate the size of your matrix where you'll store values that can
% be used to record error (Reliability, Sensitivity, G, ect..)
Error_stats_lips = zeros(Number_of_param_simulations,13);
% multiply by 2 because also iterating through the knickzone bases as well

% for each parameter combination run in the previous calibration script
% (part 1)... (for the knickzone lips)
for i = 1:Number_of_param_simulations

    % load file containing algorithm knickzone lip positions, height, and
    % parameter values used to make those selections
    [algorithm_kp_easting, algorithm_kp_northing,...
    algorithm_kp_relief,algorithm_kp_sgolay_smooth_current,...
    algorithm_kp_chi_lump_current,algorithm_kp_minkp_1_pre_current,...
    algorithm_kp_minkp_2_post_current,num_alg_kps_current] = ...
    F3_calib_load_algorithm_output_table(Knickpoint_Database_Lips_all_params,i);
    
    %% Find True Positives and False Positives
    
    % to do this, cycle through each algorithm knickzone lip for the 
    % current iteration and see if there is a calibration knickpoint within
    % its allowable error radius (specified in parameters script)
    
    % if yes, record that alg_kp# with a TP
    % if no, record that alg_kp# with a FP
    
   [TP, FP, FN,Confirmed_kp_pair] = ...
    F4_calib_alg_tp_or_fp_and_fn(error_radius, num_alg_kps_current,...
    numb_calib_kz_ref,calib_lips_northing_str_loc,...
    algorithm_kp_northing,calib_lips_easting_str_loc,algorithm_kp_easting);
    % outputs number of true positives, false positives, and false
    % negatives, and the index values of these confirmed knickzone
    % boundaries (the cells containing algorithm knickzones that were close
    % to a calibration knickzone)

    
%% Don't measure geometric accuracy because we don't have the relief measurement

    %% Calculate Sensitivity, Reliability, Adj. Sens., Adj. Reli. and Sum of Sens and Relia.
    % (we have G already from the function above)
    
    R = TP/(TP + FP); %reliability
    S = TP/(TP + FN); %sensitivity
    
    AVG_RS = (R+S)/2; % combined  (just R and S)   
    
    % Calculate the difference in number of calibrtion and algorithm kps
    diff_nkps = num_alg_kps_current - numb_calib_kz_ref;
    %% Save statistics
    Error_stats_param = horzcat(R,S,AVG_RS,TP,...
        FP,FN,num_alg_kps_current,numb_calib_kz_ref,diff_nkps,algorithm_kp_sgolay_smooth_current,...
        algorithm_kp_chi_lump_current,algorithm_kp_minkp_2_post_current,...
        algorithm_kp_minkp_1_pre_current);
     % save in a matrix for each calibration file read   
    Error_stats_lips(i,:)=Error_stats_param;
end

  
% write database to csv file with header
hdr = {'1Rel', '2Sens', '3AVG_RS', '4TP', '5FP', '6FN',...
    '7n_al_kp', '8n_cal_kp', '9diff_kp','10Sgol_w','11chi_L','12min_kp2','13min_kp1'};
txt = sprintf('%s, ',hdr{:});
txt(end)='';

% write the filename for the calibration results of knickzone lips
kp_lips_fn = strcat('calibration_database_kz_lips.csv');

% add header
cd([oldfolder,'/calib_database_outputs/'])
fid = fopen(kp_lips_fn, 'w');
fprintf(fid, txt);
fclose(fid);

% add data starting with row 2 (row 1 is headers)
dlmwrite([oldfolder '/calib_database_outputs/',kp_lips_fn],Error_stats_lips,'-append','delimiter',',', 'precision', '%.5f','roffset',1);

% accuracy of lip selections
avg_R_S_lips = Error_stats_lips(:,3);

% find the index (parameter combination) that produced the best result (for
% lips and bases)
[Best_Fit_RS, idx_BF_RS] = max(avg_R_S_lips);

best_ft_rs_params = Error_stats_lips(idx_BF_RS,:);

%% MAKE PLOTS
% plot spread of RS and RSG (accuracy) from using different parameter 
% values: for S-golay window size, lumping window size, min-kz1, and min-kz2

% plot distribution of results from changing smoothing window size
han2 = figure;set(han2, 'Visible', 'off');

hold on
plot(Error_stats_lips(:,10),avg_R_S_lips,'k.');
plot(best_ft_rs_params(:,10),avg_R_S_lips(idx_BF_RS),'rs');
xlabel('smoothing window size (cells)')
ylabel('Average of R and S only')
grid on

saveas(han2,[oldfolder '/calib_figure_outputs/','SG_window_Accuracy.png']);
clf

% plot distribution of results from changing lumping distance window size
han3 = figure;set(han3, 'Visible', 'off');

hold on
plot(Error_stats_lips(:,11),avg_R_S_lips,'k.');
plot(best_ft_rs_params(:,11),avg_R_S_lips(idx_BF_RS),'rs');
xlabel('lumping window size (m)')
ylabel('Average of R and S only')
grid on

saveas(han3,[oldfolder '/calib_figure_outputs/','Lumping_distance_accuracy.png']);
clf

% plot distribution of results from changing minimum knickzone height 2
han4 = figure;set(han4, 'Visible', 'off');

hold on
plot(Error_stats_lips(:,12),avg_R_S_lips,'k.');
plot(best_ft_rs_params(:,12),avg_R_S_lips(idx_BF_RS),'rs');
xlabel('min kz height post lumping (m)')
ylabel('Average of R and S only')
grid on

saveas(han4,[oldfolder '/calib_figure_outputs/','min_kz_height_post_lumping.png']);
clf

% plot distribution of results from changing minimum knickzone height 1
han5 = figure;set(han5, 'Visible', 'off');

hold on
plot(Error_stats_lips(:,13),avg_R_S_lips,'k.');
plot(best_ft_rs_params(:,13),avg_R_S_lips(idx_BF_RS),'rs');
xlabel('min kz height pre lumping (m)')
ylabel('Average of R and S only')
grid on

saveas(han5,[oldfolder '/calib_figure_outputs/','min_kz_height_pre_lumping.png']);
clf
%%  Plot Longitudinal profiles and maps of best fit algorithm results

% need to load the BF algorithm results

% load the lips infomation for the best fit result
% get results that maximized RS and RSG (lips)
Best_Fit_RS_lips = Knickpoint_Database_Lips_all_params{idx_BF_RS};

% Make map with BF algorithm KZs and calibration KZs (RS maximized only)

han6 = figure;set(han6, 'Visible', 'off')
hold on
grid on
% plot streamOBJ
plot(current_basin_streamobj,'k')
% plot calibration knickzone positions
h1 = plot(calib_lips_easting_str_loc,calib_lips_northing_str_loc,'b.','MarkerSize',15);

% plot best fit algorithm knickzone positions
h2 = plot(Best_Fit_RS_lips(:,12),Best_Fit_RS_lips(:,13),'r.','MarkerSize',11);
title('Best fit parameters maximize R & S')
xlabel('easting')
ylabel('northing')
legend([h1 h2],{'calibration lips','Alg. BF lips'},'Location','southoutside')

saveas(han6,[oldfolder '/calib_figure_outputs/','BF_map_RS_only.png']);
clf   

% Plot longitudinal profiles AND display the best fit parameter values and
% resulting R,S,and G

han7 = figure;set(han7, 'Visible', 'off')
hold on
% plot streamOBJ
plotdz(current_basin_streamobj,AOI_DEM,'color','k')
grid on
% plot calibration knickzone positions
h1 = plot(calib_lips_distance_str_loc,calib_lips_elevation_str_loc,'b.','MarkerSize',15);

% plot best fit algorithm knickzone positions
h2 = plot(Best_Fit_RS_lips(:,15),Best_Fit_RS_lips(:,9),'r.','MarkerSize',11);
title('Best fit parameters maximize R & S')

% add text with parameters used and accuracy of results (retrieve from BF
% index and the table storing the parameter calibration results)
R_stat = sprintf('R: %0.3f',Error_stats_lips(idx_BF_RS,1));
    str_R = sprintf('%s', R_stat);
    
S_stat = sprintf('S: %0.3f',Error_stats_lips(idx_BF_RS,2));
    str_S = sprintf('%s', S_stat);
    
SW_stat = sprintf('BF SW: %0.3f',Error_stats_lips(idx_BF_RS,10));
    str_SW = sprintf('%s', SW_stat);

CL_stat = sprintf('BF CL: %0.3f',Error_stats_lips(idx_BF_RS,11));
    str_CL = sprintf('%s', CL_stat);
    
Mkp2_stat = sprintf('BF Mkz2: %0.3f',Error_stats_lips(idx_BF_RS,12));
    str_Mkp2 = sprintf('%s', Mkp2_stat);
    
Mkp1_stat = sprintf('BF Mkz1: %0.3f',Error_stats_lips(idx_BF_RS,13));
    str_Mkp1 = sprintf('%s', Mkp1_stat);
    
descr_maxRS = {str_R;
    str_S;
    str_SW;
    str_CL;
    str_Mkp2;
    str_Mkp1};

% place the text in upper lefthand corner of longitudinal profile
text(0 + (max(current_basin_chiplot.distance))/25,max(stream_loc_elevation),descr_maxRS,'VerticalAlignment','top');
legend([h1 h2],{'calibration lips','Alg. BF lips'},'Location','southoutside')

saveas(han7,[oldfolder '/calib_figure_outputs/','BF_long_profile_RS_only.png']);
clf  

% Write best fit results table to file and store (so user could plot in
% arcmap

hdr = {'1kp_id', '2stream_id', '3ksn_SA', '4ks_chi', '5m_n', '6trib_id', ...
                '7trib_ks', '8chi', '9Elev', '10kz_mag', '11kz_rel', '12kz_E', ...
                '13kz_N', '14kz_DA', '15kz_dFm', '16kz_slp', '17kz_leng', ...
                '18sg_wz', '19sg_po', '20lwsz', '21mkz_1', ...
                '22mkz_2','23mkz_slp','24m_trib','25min_DA'};
            txt = sprintf('%s, ',hdr{:}); % change from cell array to text
            txt(end)='';
     
kp_lips_fn = strcat('kp_lips_best_fit_RS.csv');
%give filename

cd([oldfolder '\calib_lips\']);
fid = fopen(kp_lips_fn, 'w');
fprintf(fid, txt);
fclose(fid);

% write .csv file containing knickzone lip information for
% current parameters
dlmwrite([oldfolder '\calib_lips\',kp_lips_fn],Best_Fit_RS_lips,'-append','delimiter',',', 'precision', '%.5f','roffset',1);

cd(oldfolder)