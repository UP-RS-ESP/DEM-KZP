function [current_chi_resampled, current_elev_resampled,current_stepsize_mean] = ...
    E2_chi_interpolate_stepsize(current_chi_trib, current_elev_trib)

                % get chi stepsize
                for d = 2:length(current_chi_trib)
                    current_stepsize(d) = current_chi_trib(d) - current_chi_trib(d-1);
                end
                % takes the difference between each successive point and stores them in
                % current_stepsize
                current_stepsize(1) = 0; 
                % add cell to the begining because we lost a cell in that loop

                %the next line calculates averaged stepsize for chi:
                current_stepsize_mean = mean(abs(current_stepsize));

                % make a new vector with average stepsize for resampling/interpolation
                current_chi_resampled = min(current_chi_trib):current_stepsize_mean:max(current_chi_trib);

                %perform 1D interpolation
                current_elev_resampled = interp1(current_chi_trib,current_elev_trib,current_chi_resampled,'linear');

                % need to flip these and transpose them because they come out upside down
                current_chi_resampled = flipud(current_chi_resampled');
                current_elev_resampled = flipud(current_elev_resampled');