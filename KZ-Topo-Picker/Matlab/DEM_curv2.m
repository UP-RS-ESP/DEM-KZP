function AOI_DEM_curv2 = DEM_curv2(AOI_DEM, AOI_DEM_curv2_fname)

[fx, fy] = gradient(double(AOI_DEM.Z), AOI_DEM.cellsize);
AOI_DEM_slope = sqrt(fx.^2 + fy.^2); %calculate slope
fx = fx ./ (AOI_DEM_slope + 0.000001);
fy = fy ./ (AOI_DEM_slope + 0.000001);

[fx1, ~] = gradient(fx, AOI_DEM.cellsize);
[~, fy2] = gradient(fy, AOI_DEM.cellsize);
AOI_DEM_curv = fx1 + fy2;
AOI_DEM_curv2 = AOI_DEM;
AOI_DEM_curv2.Z = AOI_DEM_curv;
AOI_DEM_curv2.name = 'curvature';

if exist(AOI_DEM_curv2_fname, 'file') ~= 2
    GRIDobj2geotiff(AOI_DEM_curv2,AOI_DEM_curv2_fname);
end
