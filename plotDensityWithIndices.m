FILENAME = 'JUHO.tiegcm1.95.scntr_mareqx_smin_001.nc';
% Read density
densArray = ncread(FILENAME, 'DEN');
% Take layer of lowest altitude at first recorded time
densSurf = densArray(:,:,1,1);
% Remove fist column of density matrix, as it includes impossible values
densSurf(:,1) = [];
% Read coordinates
longitude = ncread(FILENAME, 'lon');
latitude = ncread(FILENAME, 'lat');
% Now, plot using row and column numbers, taking account that density
% matrix now has one column less
[x,y] = meshgrid(1:length(longitude), 1:length(latitude) - 1);
% plot surface plot
surf(x,y,double(densSurf'))
view(2);