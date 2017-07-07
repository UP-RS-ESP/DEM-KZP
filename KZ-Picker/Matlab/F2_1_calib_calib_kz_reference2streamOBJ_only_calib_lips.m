function [calib_lips_cells, numb_calib_kz_ref] = ...
    F2_1_calib_calib_kz_reference2streamOBJ_only_calib_lips(current_basin_streamobj, Easting_calib_lips,...
    Northing_calib_lips, current_basin_chiplot,calibration_snapping_tolerance)


% inputs: need streamOBJ & chiplot structure arrays. Need easting and
% northing position of calibration lips and bases (before referencing to
% streamobj

% outputs: cells within the chiplot structure array that contain the
% calibration knickzone lips and bases. number of calibration knickzone
% lips and bases

% this function also plots map of streamOBJ with calibration knickzones
% before referencing them to streamOBJ and plots a map of the streamOBJ
% with calibration knickzones after referencing them to the streamOBJ

%% Reference calibration knickpoints to the streamobj file
% reference using geographic position (easting and northing)
stream_loc_east = current_basin_chiplot.x;
stream_loc_north = current_basin_chiplot.y; 
% we use the chiplot struct array to get position of calibration knickzones
% in stream network

calib_idx = []; % store index of calibration pts closer to streamOBJ than 
% 'calibration_snapping_tolerance'

n = 1;
for i = 1:length(Easting_calib_lips) % for all calibration knickzone lips
    [m,idx] = min((abs(stream_loc_east-Easting_calib_lips(i))+ abs(stream_loc_north-Northing_calib_lips(i))));
    % minimize distance between streamobj and knickpoint position  
    if m < calibration_snapping_tolerance % remove calib. kz that are further from streamOBJ than snapping distance
    calib_lips_cells(n) = idx; 
    % store the index in the streamOBJ array that contains the calibration
    % knickzone cell (only store this for knickzone lips that were closer
    % than the snapping distance! -> use index (n)
    calib_idx(n) = i; 
   % save index in the calibration knickzone lip vector for the calibration pts 
   % that were close enough to streamobj
    n = n+1; % increase the counter used to store the index values
    end
    % find cells on streamobj that contain the calibration knickpoints
end

calib_lips_easting_str_loc = stream_loc_east(calib_lips_cells); 
% get the easting of the calibration knickpoint on the streamOBJ 

calib_lips_northing_str_loc = stream_loc_north(calib_lips_cells);
% get the northing of the calibration knickpoint on the streamOBJ 

% plot the referenced calibration knickzones and the un-referenced 
% calibration knickzones on the streamOBJ

% recalculate the # of calibration knickzone boundaries (now that some were
% deleted because they were too far away from the streamobj)
numb_calib_kz_ref = length(calib_lips_cells);


% use this plot to see if kncikpoints moved a lot during the referencing
% process
han = figure;set(han, 'Visible', 'off');
hold on
plot(current_basin_streamobj,'k')
xlabel('Easting')
ylabel('Northing')
% plot the northing/easting position of the calibration knickzone lips and
% bases before referencing them to the streamOBJ
for i = 1:length(Easting_calib_lips)
    h1=plot(Easting_calib_lips(i),Northing_calib_lips(i),'k.','MarkerSize',16);
end

% plot the referenced northing/easting position of the calibration
% knickzone lips and bases

for i = 1:length(calib_lips_easting_str_loc)
    h2=plot(calib_lips_easting_str_loc(i),calib_lips_northing_str_loc(i),'r.','MarkerSize',12);
end

legend([h1 h2],{'lips pre-referencing','lips post-referencing'},'Location','southoutside')
title('Calibration KZ Position Pre- and Post-Referencing to StreamOBJ')

saveas(han,[pwd '/calib_figure_outputs/','Calib_KZ_sOBJ_post-ref.png']);   
% save figure in directory
clf 
