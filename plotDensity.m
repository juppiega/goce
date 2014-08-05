FILENAME = 'JUHO.tiegcm1.95.scntr_mareqx_smin_001.nc';
% Read density
densArray = ncread(FILENAME, 'DEN');
% Take layer of lowest altitude at first recorded time
densSurf = densArray(:,:,29,5);
% Read coordinates
longitude = ncread(FILENAME, 'lon');
latitude = ncread(FILENAME, 'lat');
% Make coordinates matrices to be able to use surf
[x,y] = meshgrid(longitude, latitude);
% Plot surface plot
surf(x,y,double(densSurf'))
view(2);