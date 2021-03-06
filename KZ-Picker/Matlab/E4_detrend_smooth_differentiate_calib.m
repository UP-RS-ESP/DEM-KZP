function [smooth_detrended_elev_current,diff_smooth_detrended_elev_current] = ...
    E4_detrend_smooth_differentiate_calib(current_chi_resampled, detrended_elev_current,...
    l,SG_smoothing_calib,sgolayfilt_order,smoothing_option)

if smoothing_option == 1 % if smoothing is turned on (turn off for low res DEMs)
                                        
if license('test', 'curve_fitting_toolbox') == 1
smooth_detrended_elev_current = smooth(current_chi_resampled,detrended_elev_current,SG_smoothing_calib(l),'sgolay');  % currently we're using an sgolay filter with 50 cell smoothing
else
smooth_detrended_elev_current = sgolayfilt(detrended_elev_current, sgolayfilt_order, SG_smoothing_calib(l));
end

else
smooth_detrended_elev_current = detrended_elev_current; % if no smoothing just rename
end

% Differentiate
diff_smooth_detrended_elev_current = diff(smooth_detrended_elev_current); % differentiate
% used to find inflections in smoothed, detrended, chi elev
% data