function AOI_DEM_curvdiff = DEM_curv_diffusion(AOI_DEM)

[fx, fy] = gradient(AOI_DEM.Z, AOI_DEM.cellsize);
AOI_DEM_slope = sqrt(fx.^2 + fy.^2); %calculate slope
fx = fx ./ (AOI_DEM_slope + 0.000001);
fy = fy ./ (AOI_DEM_slope + 0.000001);

[fx1, ~] = gradient(fx, AOI_DEM.cellsize);
[~, fy2] = gradient(fy, AOI_DEM.cellsize);
AOI_DEM_curvdiff = fx1 + fy2;
%diffusionp(fx1 + fy2, 'pm2',10, 1, 0.1, 0.05);

