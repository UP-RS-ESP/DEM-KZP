function [AOI_ridgecrest_MS, AOI_ridgecrest_MS_Dy_all, AOI_ridgecrest_MS_yi_all, ridgecrest_stepsize] = ...
    DEM_ridgecrest_analyze(AOI_ridgecrest_MS, AOI_ridgecrest_DEM, DEM_HYD_MAT_fname)
% function [AOI_ridgecrest_MS, AOI_ridgecrest_MS_Dy_all, AOI_ridgecrest_MS_yi_all, ridgecrest_stepsize] = ...
%    DEM_ridgecrest_analyze(AOI_ridgecrest_MS, AOI_ridgecrest_DEM, DEM_HYD_MAT_fname)
%
% Analyze ridgecrest profile
%

if exist('AOI_ridgecrest_MS', 'var') ~= 1
    load(DEM_HYD_MAT_fname, 'AOI_ridgecrest_MS', 'AOI_ridgecrest_DEM', 'AOI_ridgecrest_diffusionf_curv');
end

% get maximum numbers of points and use this as stepsize for normalized
% vector
x_length = NaN(length(AOI_ridgecrest_MS),1);
for i = 1:length(AOI_ridgecrest_MS)
    x = AOI_ridgecrest_MS(i).distance_norm; 
    x_length(i) = length(x);
end
ridgecrest_stepsize = max(x_length);
if ridgecrest_stepsize > 1e6
    ridgecrest_stepsize = round(ridgecrest_stepsize,-5);
elseif ridgecrest_stepsize > 1e5
    ridgecrest_stepsize = round(ridgecrest_stepsize,-4);
elseif ridgecrest_stepsize > 1e4
    ridgecrest_stepsize = round(ridgecrest_stepsize,-3);
elseif ridgecrest_stepsize > 1e3
    ridgecrest_stepsize = round(ridgecrest_stepsize,-2);
elseif ridgecrest_stepsize > 1e2
    ridgecrest_stepsize = round(ridgecrest_stepsize,-1);
end
ridgecrest_stepsize = 2/ridgecrest_stepsize;
x1 = [-1:ridgecrest_stepsize:1];
x1_norm = (x1 - min(x1)) ./ (max(x1) - min(x1));

% normalize and calculate difference to parabola and cosh*2, cosh*4
y1_parabola = -x1.^2; y1_parabola_norm = y1_parabola - min(y1_parabola) ./ (max(y1_parabola) - min(y1_parabola));
y1_cosh = cosh(x1.*2); y1_cosh_norm = (y1_cosh - min(y1_cosh)) ./ (max(y1_cosh) - min(y1_cosh));
y1_cosh_norm = -y1_cosh_norm; y1_cosh2_norm = y1_cosh_norm + 1;
y1_cosh = cosh(x1.*4); y1_cosh_norm = (y1_cosh - min(y1_cosh)) ./ (max(y1_cosh) - min(y1_cosh));
y1_cosh_norm = -y1_cosh_norm; y1_cosh4_norm = y1_cosh_norm + 1;

for i = 1:length(AOI_ridgecrest_MS)
    x = AOI_ridgecrest_MS(i).distance_norm; y = AOI_ridgecrest_DEM(i).values_norm;
    if isempty(x) || isempty(y)
        AOI_ridgecrest_MS(i).Dy_parab = NaN(1, length(x1_norm));
        AOI_ridgecrest_MS(i).Dy_cosh2 = NaN(1, length(x1_norm));
        AOI_ridgecrest_MS(i).Dy_cosh4 = NaN(1, length(x1_norm));
        AOI_ridgecrest_MS(i).elev_norm = NaN(1, length(x1_norm));
        continue
    end
	yi = double(interp1(x, y, x1_norm));
    AOI_ridgecrest_MS(i).elev = y;
    AOI_ridgecrest_MS(i).elev_norm = yi;
    AOI_ridgecrest_MS(i).Dy_parab = yi-y1_parabola_norm;
    AOI_ridgecrest_MS(i).Dy_cosh2 = yi-y1_cosh2_norm;
    AOI_ridgecrest_MS(i).Dy_cosh4 = yi-y1_cosh4_norm;
    clear yi
end

%merge all profiles into one variable and get average
AOI_ridgecrest_MS_Dy_parab_all = [AOI_ridgecrest_MS.Dy_parab];
AOI_ridgecrest_MS_Dy_parab_all = reshape(AOI_ridgecrest_MS_Dy_parab_all, ...
    length(x1_norm), length(AOI_ridgecrest_MS));
AOI_ridgecrest_MS_Dy_parab_all_mean = nanmean(AOI_ridgecrest_MS_Dy_parab_all');
AOI_ridgecrest_MS_Dy_parab_all_std = nanstd(AOI_ridgecrest_MS_Dy_parab_all');

AOI_ridgecrest_MS_Dy_cosh2_all = [AOI_ridgecrest_MS.Dy_cosh2];
AOI_ridgecrest_MS_Dy_cosh2_all = reshape(AOI_ridgecrest_MS_Dy_cosh2_all, ...
    length(x1_norm), length(AOI_ridgecrest_MS));
AOI_ridgecrest_MS_Dy_cosh2_all_mean = nanmean(AOI_ridgecrest_MS_Dy_cosh2_all');
AOI_ridgecrest_MS_Dy_cosh2_all_std = nanstd(AOI_ridgecrest_MS_Dy_cosh2_all');

AOI_ridgecrest_MS_Dy_cosh4_all = [AOI_ridgecrest_MS.Dy_cosh4];
AOI_ridgecrest_MS_Dy_cosh4_all = reshape(AOI_ridgecrest_MS_Dy_cosh4_all, ...
    length(x1_norm), length(AOI_ridgecrest_MS));
AOI_ridgecrest_MS_Dy_cosh4_all_mean = nanmean(AOI_ridgecrest_MS_Dy_cosh4_all');
AOI_ridgecrest_MS_Dy_cosh4_all_std = nanstd(AOI_ridgecrest_MS_Dy_cosh4_all');

AOI_ridgecrest_MS_yi_all = [AOI_ridgecrest_MS.elev_norm];
AOI_ridgecrest_MS_yi_all = reshape(AOI_ridgecrest_MS_yi_all', ...
    length(x1_norm), length(AOI_ridgecrest_MS));
AOI_ridgecrest_MS_yi_all_mean = nanmean(AOI_ridgecrest_MS_yi_all');
AOI_ridgecrest_MS_yi_all_std = nanstd(AOI_ridgecrest_MS_yi_all');

%calculate delta (residual) for mean profile and put into structure
for i = 1:length(AOI_ridgecrest_MS)
    AOI_ridgecrest_MS(i).Dy_mean = AOI_ridgecrest_MS(i).elev_norm - AOI_ridgecrest_MS_yi_all_mean;
end
AOI_ridgecrest_MS_Dy_all = [AOI_ridgecrest_MS.Dy_mean];
AOI_ridgecrest_MS_Dy_all = reshape(AOI_ridgecrest_MS_Dy_all', ...
    length(x1_norm), length(AOI_ridgecrest_MS));
AOI_ridgecrest_MS_Dy_all_std = nanstd(AOI_ridgecrest_MS_Dy_all');

% calculate standard deviation and find profiles locations that are >2*std away
% Mean catchment height
for i = 1:length(AOI_ridgecrest_MS)
    %AOI_ridgecrest_MS(i).Dy_1std = AOI_ridgecrest_MS_Dy_all_std;
    %AOI_ridgecrest_MS(i).Dy_2std = (2*AOI_ridgecrest_MS_Dy_all_std);
    %find points that lie outside the error envelope:
    idx1 = find((AOI_ridgecrest_MS(i).elev_norm > AOI_ridgecrest_MS_yi_all_mean+AOI_ridgecrest_MS_yi_all_std) | ...
        (AOI_ridgecrest_MS(i).elev_norm < AOI_ridgecrest_MS_yi_all_mean-AOI_ridgecrest_MS_yi_all_std));
    idx2 = find((AOI_ridgecrest_MS(i).elev_norm > AOI_ridgecrest_MS_yi_all_mean+(2.*AOI_ridgecrest_MS_yi_all_std)) | ...
        (AOI_ridgecrest_MS(i).elev_norm < AOI_ridgecrest_MS_yi_all_mean-(2.*AOI_ridgecrest_MS_yi_all_std)));
    % now find corresponding coordinates/distance along profile
    % distance stored in x1_norm
    
    % convert/transform indices to distance vector
    idx1_distance = round((idx1./length(x1)) .* length(AOI_ridgecrest_MS(i).x));
    idx2remove = find(isinf(idx1_distance) | idx1_distance == 0);
    if ~isempty(idx2remove)
        idx1_distance(idx2remove) = [];
        idx1(idx2remove) = [];
    end
    clear idx2remove
    idx2_distance = round((idx2./length(x1)) * length(AOI_ridgecrest_MS(i).x));
    idx2remove = find(isinf(idx2_distance) | idx2_distance == 0);
    if ~isempty(idx2remove)
        idx2_distance(idx2remove) = [];
        idx2(idx2remove) = [];
    end
    AOI_ridgecrest_MS(i).Dy_mean_lg_1std = idx1_distance;
    AOI_ridgecrest_MS(i).distance_Dy_mean_lg_1std = AOI_ridgecrest_MS(i).distance(idx1_distance)';
    flag = zeros(1, length(AOI_ridgecrest_MS(i).distance));
    flag(idx1_distance) = 1;
    AOI_ridgecrest_MS(i).Dy_mean_lg_1std_flag = flag;
    AOI_ridgecrest_MS(i).Dyi_mean_lg_1std = idx1;
    
    AOI_ridgecrest_MS(i).Dy_mean_lg_2std = idx2_distance;
    AOI_ridgecrest_MS(i).distance_Dy_mean_lg_2std = AOI_ridgecrest_MS(i).distance(idx2_distance)';
    flag = zeros(1, length(AOI_ridgecrest_MS(i).distance));
    flag(idx2_distance) = 1;
    AOI_ridgecrest_MS(i).Dy_mean_lg_2std_flag = flag;
    AOI_ridgecrest_MS(i).Dyi_mean_lg_2std = idx2;

    %DISTANCE TO PARABOLA:
    idx1 = find((AOI_ridgecrest_MS(i).Dy_parab > AOI_ridgecrest_MS_Dy_parab_all_mean+AOI_ridgecrest_MS_Dy_parab_all_std) | ...
        (AOI_ridgecrest_MS(i).Dy_parab < AOI_ridgecrest_MS_Dy_parab_all_mean-AOI_ridgecrest_MS_Dy_parab_all_std));
    idx2 = find((AOI_ridgecrest_MS(i).Dy_parab > AOI_ridgecrest_MS_Dy_parab_all_mean+(2.*AOI_ridgecrest_MS_Dy_parab_all_std)) | ...
        (AOI_ridgecrest_MS(i).Dy_parab < AOI_ridgecrest_MS_Dy_parab_all_mean-(2.*AOI_ridgecrest_MS_Dy_parab_all_std)));
    % convert/transform indices to distance vector
    idx1_distance = round((idx1./length(x1)) * length(AOI_ridgecrest_MS(i).x));
    idx2remove = find(isinf(idx1_distance) | idx1_distance == 0);
    if ~isempty(idx2remove)
        idx1_distance(idx2remove) = [];
        idx1(idx2remove) = [];
    end
    clear idx2remove
    idx2_distance = round((idx2./length(x1)) * length(AOI_ridgecrest_MS(i).x));
    idx2remove = find(isinf(idx2_distance) | idx2_distance == 0);
    if ~isempty(idx2remove)
        idx2_distance(idx2remove) = [];
        idx2(idx2remove) = [];
    end
    AOI_ridgecrest_MS(i).Dy_parab_lg_1std = idx1_distance;
    AOI_ridgecrest_MS(i).distance_Dy_parab_lg_1std = AOI_ridgecrest_MS(i).distance(idx1_distance)';
    flag = zeros(1, length(AOI_ridgecrest_MS(i).distance));
    flag(idx1_distance) = 1;
    AOI_ridgecrest_MS(i).Dy_parab_lg_1std_flag = flag;
	AOI_ridgecrest_MS(i).Dyi_parab_lg_1std = idx1;
    AOI_ridgecrest_MS(i).Dy_parab_lg_2std = idx2_distance;
    AOI_ridgecrest_MS(i).distance_Dy_parab_lg_2std = AOI_ridgecrest_MS(i).distance(idx2_distance)';
    flag = zeros(1, length(AOI_ridgecrest_MS(i).distance));
    flag(idx2_distance) = 1;
    AOI_ridgecrest_MS(i).Dy_parab_lg_2std_flag = flag;
	AOI_ridgecrest_MS(i).Dyi_parab_lg_2std = idx2;

    %DISTANCE TO cosh.*2:
    idx1 = find((AOI_ridgecrest_MS(i).Dy_cosh2 > AOI_ridgecrest_MS_Dy_cosh2_all_mean+AOI_ridgecrest_MS_Dy_cosh2_all_std) | ...
        (AOI_ridgecrest_MS(i).Dy_cosh2 < AOI_ridgecrest_MS_Dy_cosh2_all_mean-AOI_ridgecrest_MS_Dy_cosh2_all_std));
    idx2 = find((AOI_ridgecrest_MS(i).Dy_cosh2 > AOI_ridgecrest_MS_Dy_cosh2_all_mean+(2.*AOI_ridgecrest_MS_Dy_cosh2_all_std)) | ...
        (AOI_ridgecrest_MS(i).Dy_cosh2 < AOI_ridgecrest_MS_Dy_cosh2_all_mean-(2.*AOI_ridgecrest_MS_Dy_cosh2_all_std)));
    % convert/transform indices to distance vector
    idx1_distance = round((idx1./length(x1)) * length(AOI_ridgecrest_MS(i).x));
    idx2remove = find(isinf(idx1_distance) | idx1_distance == 0);
    if ~isempty(idx2remove)
        idx1_distance(idx2remove) = [];
        idx1(idx2remove) = [];
    end
    clear idx2remove
    idx2_distance = round((idx2./length(x1)) * length(AOI_ridgecrest_MS(i).x));
    idx2remove = find(isinf(idx2_distance) | idx2_distance == 0);
    if ~isempty(idx2remove)
        idx2_distance(idx2remove) = [];
        idx2(idx2remove) = [];
    end
    AOI_ridgecrest_MS(i).Dy_cosh2_lg_1std = idx1_distance;
    AOI_ridgecrest_MS(i).distance_Dy_cosh2_lg_1std = AOI_ridgecrest_MS(i).distance(idx1_distance)';
    flag = zeros(1, length(AOI_ridgecrest_MS(i).distance));
    flag(idx1_distance) = 1;
    AOI_ridgecrest_MS(i).Dy_cosh2_lg_1std_flag = flag;
    AOI_ridgecrest_MS(i).Dyi_cosh2_lg_1std = idx1;
    AOI_ridgecrest_MS(i).Dy_cosh2_lg_2std = idx2_distance;
    AOI_ridgecrest_MS(i).distance_Dy_cosh2_lg_2std = AOI_ridgecrest_MS(i).distance(idx2_distance)';
    flag = zeros(1, length(AOI_ridgecrest_MS(i).distance));
    flag(idx2_distance) = 1;
    AOI_ridgecrest_MS(i).Dy_cosh2_lg_2std_flag = flag;
    AOI_ridgecrest_MS(i).Dyi_cosh2_lg_2std = idx2;
    
    %DISTANCE TO cosh.*4:
    idx1 = find((AOI_ridgecrest_MS(i).Dy_cosh4 > AOI_ridgecrest_MS_Dy_cosh4_all_mean+AOI_ridgecrest_MS_Dy_cosh4_all_std) | ...
        (AOI_ridgecrest_MS(i).Dy_cosh4 < AOI_ridgecrest_MS_Dy_cosh4_all_mean-AOI_ridgecrest_MS_Dy_cosh4_all_std));
    idx2 = find((AOI_ridgecrest_MS(i).Dy_cosh4 > AOI_ridgecrest_MS_Dy_cosh4_all_mean+(2.*AOI_ridgecrest_MS_Dy_cosh4_all_std)) | ...
        (AOI_ridgecrest_MS(i).Dy_cosh4 < AOI_ridgecrest_MS_Dy_cosh4_all_mean-(2.*AOI_ridgecrest_MS_Dy_cosh4_all_std)));
    % convert/transform indices to distance vector
    idx1_distance = round((idx1./length(x1)) * length(AOI_ridgecrest_MS(i).x));
    idx2remove = find(isinf(idx1_distance) | idx1_distance == 0);
    if ~isempty(idx2remove)
        idx1_distance(idx2remove) = [];
        idx1(idx2remove) = [];
    end
    clear idx2remove
    idx2_distance = round((idx2./length(x1)) * length(AOI_ridgecrest_MS(i).x));
    idx2remove = find(isinf(idx2_distance) | idx2_distance == 0);
    if ~isempty(idx2remove)
        idx2_distance(idx2remove) = [];
        idx2(idx2remove) = [];
    end
    AOI_ridgecrest_MS(i).Dy_cosh4_lg_1std = idx1_distance;
    AOI_ridgecrest_MS(i).distance_Dy_cosh4_lg_1std = AOI_ridgecrest_MS(i).distance(idx1_distance)';
    flag = zeros(1, length(AOI_ridgecrest_MS(i).distance));
    flag(idx1_distance) = 1;
    AOI_ridgecrest_MS(i).Dy_cosh4_lg_1std_flag = flag;
    AOI_ridgecrest_MS(i).Dyi_cosh4_lg_1std = idx1;
    AOI_ridgecrest_MS(i).Dy_cosh4_lg_2std = idx2_distance;
    AOI_ridgecrest_MS(i).distance_Dy_cosh4_lg_2std = AOI_ridgecrest_MS(i).distance(idx2_distance)';
    flag = zeros(1, length(AOI_ridgecrest_MS(i).distance));
    flag(idx2_distance) = 1;
    AOI_ridgecrest_MS(i).Dy_cosh4_lg_2std_flag = flag;
	AOI_ridgecrest_MS(i).Dyi_cosh4_lg_2std = idx2;
end


%%
%Plot map view of 1-std.dev. and 2-std.dev.
% figure
% grid on
% hold on 
% for i = 1:length(AOI_ridgecrest_MS)
% scatter(AOI_ridgecrest_MS(i).x(AOI_ridgecrest_MS(i).Dy_lg_1std), ...
%     AOI_ridgecrest_MS(i).y(AOI_ridgecrest_MS(i).Dy_lg_1std), ...
%     abs(AOI_ridgecrest_MS(i).Dy(AOI_ridgecrest_MS(i).Dy_lg_1std).*1000), ...
%     AOI_ridgecrest_MS(i).Dy(AOI_ridgecrest_MS(i).Dy_lg_1std))
% end
% 
% 
% i = 1;
% figure
% plot(AOI_ridgecrest_MS(i).y, AOI_ridgecrest_MS(i).x, '.')
% hold on
% grid on
% plot(AOI_ridgecrest_MS(i).y(AOI_ridgecrest_MS(i).Dy_lg_1std),...
%     AOI_ridgecrest_MS(i).x(AOI_ridgecrest_MS(i).Dy_lg_1std), 'bx')
% plot(AOI_ridgecrest_MS(i).y(AOI_ridgecrest_MS(i).Dy_lg_2std),...
%     AOI_ridgecrest_MS(i).x(AOI_ridgecrest_MS(i).Dy_lg_2std), 'rs')
% 
% figure
% clf
% subplot(3,1,1,'align')
% plot(x1_norm, AOI_ridgecrest_MS(i).elev_norm); hold on; grid
% plot(x1_norm, y1_norm);
% plot(x1_norm, AOI_ridgecrest_MS_yi_all_mean);
% legend('elev_norm', 'y1\_norm', 'AOI\_ridgecrest\_MS\_yi\_all\_mean');
% 
% subplot(3,1,2,'align')
% h1 = plot(x1_norm, AOI_ridgecrest_MS_Dy_all_mean, 'k-', 'Linewidth', 2); hold on; grid on
% h2 = plot(x1_norm, AOI_ridgecrest_MS_Dy_all_mean+AOI_ridgecrest_MS_Dy_all_std, 'k-', 'Linewidth', 1); 
% plot(x1_norm, AOI_ridgecrest_MS_Dy_all_mean-AOI_ridgecrest_MS_Dy_all_std, 'k-', 'Linewidth', 1); 
% h3 = plot(x1_norm, AOI_ridgecrest_MS(i).Dy_mean, 'r-', 'Linewidth', 2);
% legend([h1 h2 h3], 'Ridgecrest Dy all mean', 'Ridgecrest Dy all mean +/- std', 'current Dy mean')
% 
% subplot(3,1,3,'align')
% plot(x1_norm, AOI_ridgecrest_MS(i).Dy_mean, 'k-', 'Linewidth', 2); hold on; grid
% plot(x1_norm, AOI_ridgecrest_MS(i).Dy_mean+AOI_ridgecrest_MS(i).Dy_1std, 'k-', 'Linewidth', 1);
% plot(x1_norm, AOI_ridgecrest_MS(i).Dy_mean-AOI_ridgecrest_MS(i).Dy_1std, 'k-', 'Linewidth', 1);
% plot(x1_norm, AOI_ridgecrest_MS(i).Dy_mean+AOI_ridgecrest_MS(i).Dy_2std, 'r-', 'Linewidth', 1);
% plot(x1_norm, AOI_ridgecrest_MS(i).Dy_mean-AOI_ridgecrest_MS(i).Dy_2std, 'r-', 'Linewidth', 1);
% plot(x1_norm(AOI_ridgecrest_MS(i).Dy_lg_1std),...
%     AOI_ridgecrest_MS(i).Dy_mean(AOI_ridgecrest_MS(i).Dy_lg_1std), 'bx');
% plot(x1_norm(AOI_ridgecrest_MS(i).Dy_lg_2std),...
%     AOI_ridgecrest_MS(i).Dy_mean(AOI_ridgecrest_MS(i).Dy_lg_2std), 'rs');
% 
% % ydelta = yi-y1_norm;
% % plot(x1_norm,ydelta); hold on; grid
% % plot([min(x1_norm) max(x1_norm)], [nanmean(yi-y1_norm) nanmean(yi-y1_norm)], 'k-', 'Linewidth', 2);
% % plot([min(x1_norm) max(x1_norm)], [nanmean(yi-y1_norm)+nanstd(yi-y1_norm) nanmean(yi-y1_norm)+nanstd(yi-y1_norm)], 'k-.', 'Linewidth', 1);
% % plot([min(x1_norm) max(x1_norm)], [nanmean(yi-y1_norm)-nanstd(yi-y1_norm) nanmean(yi-y1_norm)-nanstd(yi-y1_norm)], 'k-.', 'Linewidth', 1);


% figure;
% clf
% subplot(2,2,1,'align')
% plot(x1_norm, AOI_ridgecrest_MS_yi_all_mean, 'k-', 'Linewidth', 2)
% hold on; grid on
% plot(x1_norm, AOI_ridgecrest_MS_yi_all_mean+AOI_ridgecrest_MS_yi_all_std, 'k-', 'Linewidth', 1)
% plot(x1_norm, AOI_ridgecrest_MS_yi_all_mean-AOI_ridgecrest_MS_yi_all_std, 'k-', 'Linewidth', 1)
% plot(x1_norm, AOI_ridgecrest_MS(i).elev_norm, 'r-', 'Linewidth', 2)
% plot(x1_norm(idx1),  AOI_ridgecrest_MS(i).elev_norm(idx1), 'rs')
% plot(x1_norm(idx2),  AOI_ridgecrest_MS(i).elev_norm(idx2), 'bs')
% title('\Delta Mean profile', 'Fontsize', 14)
% 
% subplot(2,2,2,'align')
% plot(x1_norm, AOI_ridgecrest_MS_Dy_parab_all_mean, 'k-', 'Linewidth', 2)
% hold on; grid on
% plot(x1_norm, AOI_ridgecrest_MS_Dy_parab_all_mean+AOI_ridgecrest_MS_Dy_parab_all_std, 'k-', 'Linewidth', 1)
% plot(x1_norm, AOI_ridgecrest_MS_Dy_parab_all_mean-AOI_ridgecrest_MS_Dy_parab_all_std, 'k-', 'Linewidth', 1)
% plot(x1_norm, AOI_ridgecrest_MS(i).Dy_parab, 'r-', 'Linewidth', 2)
% plot(x1_norm(idx1),  AOI_ridgecrest_MS(i).Dy_parab(idx1), 'rs')
% plot(x1_norm(idx2),  AOI_ridgecrest_MS(i).Dy_parab(idx2), 'bs')
% title('\Delta Parabola profile', 'Fontsize', 14)
% 
% subplot(2,2,3,'align')
% plot(x1_norm, AOI_ridgecrest_MS_Dy_cosh2_all_mean, 'k-', 'Linewidth', 2)
% hold on; grid on
% plot(x1_norm, AOI_ridgecrest_MS_Dy_cosh2_all_mean+AOI_ridgecrest_MS_Dy_cosh2_all_std, 'k-', 'Linewidth', 1)
% plot(x1_norm, AOI_ridgecrest_MS_Dy_cosh2_all_mean-AOI_ridgecrest_MS_Dy_cosh2_all_std, 'k-', 'Linewidth', 1)
% plot(x1_norm, AOI_ridgecrest_MS(i).Dy_cosh2, 'r-', 'Linewidth', 2)
% plot(x1_norm(idx1),  AOI_ridgecrest_MS(i).Dy_cosh2(idx1), 'rs')
% plot(x1_norm(idx2),  AOI_ridgecrest_MS(i).Dy_cosh2(idx2), 'bs')
% title('\Delta Cosh*2 profile', 'Fontsize', 14)
% 
% subplot(2,2,4,'align')
% plot(x1_norm, AOI_ridgecrest_MS_Dy_cosh4_all_mean, 'k-', 'Linewidth', 2)
% hold on; grid on
% plot(x1_norm, AOI_ridgecrest_MS_Dy_cosh4_all_mean+AOI_ridgecrest_MS_Dy_cosh4_all_std, 'k-', 'Linewidth', 1)
% plot(x1_norm, AOI_ridgecrest_MS_Dy_cosh4_all_mean-AOI_ridgecrest_MS_Dy_cosh4_all_std, 'k-', 'Linewidth', 1)
% plot(x1_norm, AOI_ridgecrest_MS(i).Dy_cosh4, 'r-', 'Linewidth', 2)
% plot(x1_norm(idx1),  AOI_ridgecrest_MS(i).Dy_cosh4(idx1), 'rs')
% plot(x1_norm(idx2),  AOI_ridgecrest_MS(i).Dy_cosh4(idx2), 'bs')
% title('\Delta Cosh*4 profile', 'Fontsize', 14)
