FILENAME = 'tiegcm_init.nc';
original = 'TGCM.tiegcm1.95_dres.p_whi2008_weimer_imf.nc';
varname = 'VN_NM';
% Read var
varMyInit = ncread(FILENAME, varname);
varCompare = ncread(original, varname);
level = size(varCompare, 3)-1;%1;
% Take layer of lowest altitude at first recorded time
varMyInitSurf = varMyInit(:,:,level,1);
varCompareSurf = varCompare(:,:,level,1);
% Read coordinates
longitude = ncread(FILENAME, 'lon');
latitude = ncread(FILENAME, 'lat');
% Make coordinates matrices to be able to use surf
[x,y] = meshgrid(longitude, latitude);

% Plot surface plot
subplot(2,1,1)
surf(x,y,double(varMyInitSurf'), 'linestyle', 'none')
view(2);
title('Empirical Init');
colorbar
axis('tight')

subplot(2,1,2)
surf(x,y,double(varCompareSurf'), 'linestyle', 'none')
view(2);
title('Original Init');
colorbar
axis('tight')