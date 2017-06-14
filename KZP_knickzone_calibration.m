%% (3) Optional KZP Calibration
%
fprintf(1,'KZP identifying knickzones calibration step 1 of 1: calibration for %s\n', KZP_parameters.DEM_fname);
if KZP_parameters.Calibration_option == 0
    fprintf(1,'no calibration, continue\n');
    return
elseif KZP_parameters.Calibration_option == 1
    KZP_calibration_1;
    if KZP_parameters.do_you_have_calibration_KZ_bases == 1
        KZP_calib_KZ_calibration_comparison;
    else
        KZP_calib_KZ_calibration_comparison_only_calib_lips
    end
end
