FILENAME = 'tiegcm_dres.sec_goce_spinup.nc';
original = 'TGCM.tiegcm1.95_dres.p_whi2008_weimer_imf.nc';
varname = 'OP';
% Read var
varMyInit = ncread(FILENAME, varname);
if strcmpi(varname,'OP')
    varMyInit = log10(varMyInit);
    varname = 'log_{10} [O+]';
end
%varCompare = ncread(original, varname);
%level = size(varCompare, 3)-1;%1;

varMyInitSurf = varMyInit(:,:,end-1,end);
%varCompareSurf = varCompare(:,:,level,1);
% Read coordinates
longitude = ncread(FILENAME, 'lon');
latitude = ncread(FILENAME, 'lat');
% Make coordinates matrices to be able to use surf
[x,y] = meshgrid(longitude, latitude);

% Plot surface plot
%subplot(2,1,1)
figure('color', 'w')
surf(x,y,double(varMyInitSurf'), 'linestyle', 'none')
view(2);
title(varname);
colorbar
axis('tight')
xlabel('lon')
ylabel('lat')

% subplot(2,1,2)
% surf(x,y,double(varCompareSurf'), 'linestyle', 'none')
% view(2);
% title('Original Init');
% colorbar
% axis('tight')

solarTimes = [0, 6, 12, 18];
mtime = ncread(FILENAME, 'mtime');
utNow = double(mtime(2,end)) + double(mtime(3,end))/60;
alt = ncread(FILENAME, 'Z');
alt(:,:,end,:) = [];
varMyInit(:,:,end,:) = [];
figure('color','w','units','normalized','outerposition',[0 0 1 1]);

for i = 1:length(solarTimes)
    lstLon = 15 * (solarTimes(i) - utNow);
    lstLon(lstLon<-180) = lstLon(lstLon<-180) + 360;
    lstLon(lstLon>180) = lstLon(lstLon>180) - 360;
    [~,lonInd] = min(abs(longitude - lstLon));
    
    altSlice = alt(lonInd,:,:,end);
    altSlice = reshape(altSlice,size(alt,2),size(alt,3));
    altSlice = flipud(altSlice') / 1E5;
    
    fieldSlice = varMyInit(lonInd,:,:,end);
    fieldSlice = reshape(fieldSlice,size(alt,2),size(alt,3));
    fieldSlice = flipud(fieldSlice');
    
    [latGrid,~] = meshgrid(latitude, 1:size(alt,3));
    
    subplot(2,2,i);
    surf(latGrid, altSlice, double(fieldSlice), 'linestyle', 'none', 'edgecolor', 'none')
    view(2);
    title([varname, ' @', num2str(solarTimes(i)), ' LT']);
    grid off;
    xlim([min(latitude) max(latitude)])
    %ylim([100 600])
    colorbar;
    xlabel('Geographic Latitude', 'fontsize', 14)
    ylabel('Geopotential Height [km]', 'fontsize', 14)
    hold all;
end
hold off;