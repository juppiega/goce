FILENAME = 'tiegcm_dres.sec_goce_spinup.nc';
varname = 'TE';
% Read variable
var = ncread(FILENAME, varname);
if strcmpi(varname,'OP') % Take log, if OP
    var = log10(var);
    varname = 'log_{10} [O+]';
end

levels = ncread(FILENAME, 'ilev'); % Load level coordinate.
levelToPlot = 6.75; % Plot this level.
levelInd = find(abs(levels - levelToPlot) < eps(100)); % Pressure level index.
timeInd = size(var, 4); % Time index
varLevelSurf = var(:,:,levelInd,timeInd); % Take a slice at a specific level and time.

% Read coordinates
longitude = ncread(FILENAME, 'lon');
latitude = ncread(FILENAME, 'lat');
% Make coordinates matrices to be able to use surf
[x,y] = meshgrid(longitude, latitude);

% Plot surface plot
figure('color', 'w')
% Convert single -> double and transpose to match the coordinate grids.
surf(x,y,double(varLevelSurf'), 'linestyle', 'none')
view(2); % "Top-down" view
title([varname, ' @ level ', num2str(levelToPlot)]);
colorbar
axis('tight')
xlabel('lon')
ylabel('lat')

% Plot solar time slices.

solarTimes = [0, 6, 12, 18]; % Plot these LTs
mtime = ncread(FILENAME, 'mtime');
utNow = double(mtime(2,timeInd)) + double(mtime(3,timeInd))/60; % UT in hours.
alt = ncread(FILENAME, 'Z');
% Unfortunately, the top most layer needs to be removed, since its values could be undefined.
alt(:,:,end,:) = []; 
var(:,:,end,:) = [];
figure('color','w','units','normalized','outerposition',[0 0 1 1]); % Full screen figure.

for i = 1:length(solarTimes)
     % Longitude of requested solar time.
    lstLon = 15 * (solarTimes(i) - utNow); 
    lstLon(lstLon<-180) = lstLon(lstLon<-180) + 360;
    lstLon(lstLon>180) = lstLon(lstLon>180) - 360;
    [~,lonInd] = min(abs(longitude - lstLon)); % lonInd = nearest model longitude to the requested solar time.
    
    altSlice = alt(lonInd,:,:,timeInd); % Slice geopotential grid at specific longitude and time.
    altSlice = reshape(altSlice,size(alt,2),size(alt,3)); % Reshape from 1 x nlat x nlev-1 to nlat x nlev-1
    altSlice = flipud(altSlice') / 1E5; % Transpose and flip upside down for plotting. Convert to km.
    
    % Same as above.
    varSlice = var(lonInd,:,:,timeInd);
    varSlice = reshape(varSlice,size(alt,2),size(alt,3));
    varSlice = flipud(varSlice');
    
    [latGrid,~] = meshgrid(latitude, 1:size(alt,3)); % nlev-1 x nlat latitude-grid for plotting.
    
    % Plot.
    subplot(2,2,i);
    surf(latGrid, altSlice, double(varSlice), 'linestyle', 'none', 'edgecolor', 'none')
    view(2);
    title([varname, ' @ ', num2str(solarTimes(i)), ' LT']);
    grid off;
    xlim([min(latitude) max(latitude)])
    colorbar;
    xlabel('Geographic Latitude', 'fontsize', 14)
    ylabel('Geopotential Height [km]', 'fontsize', 14)
    hold all;
end
hold off;