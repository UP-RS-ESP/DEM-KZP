function [beta_current_trib, detrended_elev_current] = ...
    E3_reference_profile_construct_detrend(current_elev_resampled,...
    current_chi_resampled,current_stepsize_mean)

% regression in chi/elev space taken from differencing upstream-most
% point and mouth of stream
mx_b = (max(current_elev_resampled)-min(current_elev_resampled))/(max(current_chi_resampled) - min(current_chi_resampled)); 
% slope between headwater point and mouth point

blf_y_current_tribs_f = []; % store reference profile points
blf_y_current_tribs = []; % store reference profile points

n =[];
for n = 1:length(current_chi_resampled)
blf_y_current_tribs(n,1) = min(current_elev_resampled)+(mx_b*current_stepsize_mean*(n-1)); 
% get the elev coordinates for the trendline
end
blf_y_current_tribs_f = flipud(blf_y_current_tribs);
% need to flip vector

% save slope of reference profile in chi elev space
beta_current_trib = mx_b;

detrended_elev_current = current_elev_resampled - blf_y_current_tribs_f;  
% perform detrending on resampled data