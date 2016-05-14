function AOI_DEM_diffusionf = ...
    DEM_diffusionf(AOI_DEM, difkernelWidth, difSSquared, ...
    difTimeIncrement, difFilterI, difMethod, AOI_DEM_diff_fname, REGEN)
% alternative way of calculating curvature by first smooth the DEM before
% curvature calculation

%smoothing operations adapted from GeoNet
%(https://sites.google.com/site/geonethome/home) based on Guy Gilboa,
%Catte et al. (1992), Perona and Malik (1990), Passalacqua et al. (2010a,
%2010b).
fprintf('\tperform diffusion filtering on DEM\n');

[fx, fy] = gradient(AOI_DEM.Z, AOI_DEM.cellsize);
AOI_DEM_slope = sqrt(fx.^2 + fy.^2); %calculate slope
difedgeT = prctile(abs(AOI_DEM_slope(:)),90);

[Ny,Nx]=size(AOI_DEM.Z);
pixelSize = AOI_DEM.cellsize;
demArray = double(AOI_DEM.Z);

for i=1:difFilterI
    originalDemArray = demArray;
    halfdifkernelWidth=(difkernelWidth-1)/2;
    % Make a ramp array with 5 rows each containing [-2, -1, 0, 1, 2]
    x=ones(difkernelWidth,1)*(-halfdifkernelWidth:halfdifkernelWidth);
    y=x';
    gaussianFilter=exp(-(x.^2+y.^2)/(2*difSSquared));  % 2D Gaussian
    gaussianFilter=gaussianFilter/sum(sum(gaussianFilter));      % Normalize
    
    xL=mean(demArray(:,1:halfdifkernelWidth)')';
    xR=mean(demArray(:,Nx-halfdifkernelWidth+1:Nx)')';
    eI=[xL*ones(1,halfdifkernelWidth) demArray xR*ones(1,halfdifkernelWidth)];
    xU=mean(eI(1:halfdifkernelWidth,:));
    xD=mean(eI(Ny-halfdifkernelWidth+1:Ny,:));
    eI=[ones(halfdifkernelWidth,1)*xU; eI; ones(halfdifkernelWidth,1)*xD];
    AOI_DEM_smooth=conv2(eI,gaussianFilter,'valid');
    
    % Calculate gradient in all directions (N,S,E,W) by simple differencing
    % - with repeat padding in each direction.
    % This step will propagate NaNs one pixel inward in each dirn.
    In=([AOI_DEM_smooth(1,:); AOI_DEM_smooth(1:Ny-1,:)]-AOI_DEM_smooth)./pixelSize;
    Is=([AOI_DEM_smooth(2:Ny,:); AOI_DEM_smooth(Ny,:)]-AOI_DEM_smooth)./pixelSize;
    Ie=([AOI_DEM_smooth(:,2:Nx) AOI_DEM_smooth(:,Nx)]-AOI_DEM_smooth)./pixelSize;
    Iw=([AOI_DEM_smooth(:,1) AOI_DEM_smooth(:,1:Nx-1)]-AOI_DEM_smooth)./pixelSize;
    In(isnan(In))=0;
    Is(isnan(Is))=0;
    Ie(isnan(Ie))=0;
    Iw(isnan(Iw))=0;
    
    % Calculate diffusion coefficients in all dirns according to diffusionMethod
    if strcmp(difMethod,'linear')
        Cn=difedgeT;
        Cs=difedgeT;
        Ce=difedgeT;
        Cw=difedgeT;
    elseif strcmp(difMethod,'PeronaMalik1')
        Cn=exp(-(abs(In)/difedgeT).^2);
        Cs=exp(-(abs(Is)/difedgeT).^2);
        Ce=exp(-(abs(Ie)/difedgeT).^2);
        Cw=exp(-(abs(Iw)/difedgeT).^2);
    elseif strcmp(difMethod, 'PeronaMalik2')
        Cn=1./(1+(abs(In)/difedgeT).^2);
        Cs=1./(1+(abs(Is)/difedgeT).^2);
        Ce=1./(1+(abs(Ie)/difedgeT).^2);
        Cw=1./(1+(abs(Iw)/difedgeT).^2);
    elseif strcmp(difMethod, 'Tukey')
        [lx ly] = find(In < difedgeT);
        chr = sparse(lx, ly, 1, length(In), length(In(1,:)));
        Cn = .5 * ((1 - (In/difedgeT).^2).^2) .*chr;
        [lx ly] = find(Is < difedgeT);
        chr = sparse(lx, ly, 1, length(Is), length(Is(1,:)));
        Cs = .5 * ((1 - (Is/difedgeT).^2).^2) .*chr;
        [lx ly] = find(Ie < difedgeT);
        chr = sparse(lx, ly, 1, length(Ie), length(Ie(1,:)));
        Ce = .5 * ((1 - (Ie/difedgeT).^2).^2) .*chr;
        [lx ly] = find(Iw < difedgeT);
        chr = sparse(lx, ly, 1, length(Iw), length(Iw(1,:)));
        Cw = .5 * ((1 - (Iw/difedgeT).^2).^2) .*chr;
    elseif strcmp(difMethod, 'rampPreserving')
        k=difedgeT(1); theta=difedgeT(2); j=sqrt(-1);
        Cn=exp(j*theta)./(1+(imag(In)/(k*theta)).^2);
        Cs=exp(j*theta)./(1+(imag(Is)/(k*theta)).^2);
        Ce=exp(j*theta)./(1+(imag(Ie)/(k*theta)).^2);
        Cw=exp(j*theta)./(1+(imag(Iw)/(k*theta)).^2);
    else
        error(['Unknown smoothing method "' difMethod '"']);
    end
    
    if (difSSquared>0)
        % Calculate real gradients (not smoothed) - with repeat padding in each
        % direction.  This step will propagate NaNs one pixel inward in each dirn.
        In=([originalDemArray(1,:); originalDemArray(1:Ny-1,:)]-originalDemArray)./pixelSize;
        Is=([originalDemArray(2:Ny,:); originalDemArray(Ny,:)]-originalDemArray)./pixelSize;
        Ie=([originalDemArray(:,2:Nx) originalDemArray(:,Nx)]-originalDemArray)./pixelSize;
        Iw=([originalDemArray(:,1) originalDemArray(:,1:Nx-1)]-originalDemArray)./pixelSize;
        In(isnan(In))=0;
        Is(isnan(Is))=0;
        Ie(isnan(Ie))=0;
        Iw(isnan(Iw))=0;
        demArray=originalDemArray;
    end
    demArray=demArray+difTimeIncrement*(Cn.*In + Cs.*Is + Ce.*Ie + Cw.*Iw);
end
AOI_DEM_diffusionf = AOI_DEM;
AOI_DEM_diffusionf.Z = demArray;
AOI_DEM_diffusionf.name = 'Diffusion-filtered DEM';

if exist(AOI_DEM_diff_fname, 'file') ~= 2 || REGEN == 1
    GRIDobj2geotiff(AOI_DEM_diffusionf,AOI_DEM_diff_fname);
end

%diffusionp(fx1 + fy2, 'pm2',10, 1, 0.2, 0.05);

