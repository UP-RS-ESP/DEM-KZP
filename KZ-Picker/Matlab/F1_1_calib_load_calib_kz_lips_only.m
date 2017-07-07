function [Easting_calib_lips, Northing_calib_lips, numb_calib_kps] = ...
    F1_1_calib_load_calib_kz_lips_only(KZ_lips_calib_fname, KZ_calib_easting_column_num,...
    KZ_calib_northing_column_num)

% ^^ inputs: filename for calibration lips and calibration bases, columns
% within that file containing the kz easting, northing, and relief

% ^^ outputs: vectors containing easting,northing, and relief of
% calibration knickzone lips and bases.

KZ_lips_calib= csvread(KZ_lips_calib_fname,1);
% load attribute table for the calibration knickzone bases and lips. ROW 1
% is assumed to be the header! so starts reading at ROW 2.


% index out the easting coordinates of kz lips and bases
Easting_calib_lips = KZ_lips_calib(:,KZ_calib_easting_column_num);
% easting

% index out the northing coordinates of kz lips and bases
Northing_calib_lips = KZ_lips_calib(:,KZ_calib_northing_column_num);
% northing

% record the number of calibration knickzones
numb_calib_kps= length(Easting_calib_lips);
% number of calibration knickpoints