function plotDensityLatitudeTimeSurf(firstDatenum, magneticLatitude, timestamps10s, regriddedLatitude, regriddedTime, regriddedGoceDensity, regriddedMsisDensity, timeOfDay)
% plotDensityLatitudeTimeSurf(averagedDensity, averagedLatitude, timestamps

persistent colormapFigHandle

if ~isempty(strfind(lower(timeOfDay), 'morning'))
    colormapFigHandle = figure('units','normalized','outerposition',[0 0 1 1]);
    goceDensitySubplot = 1;
    msisDensitySubplot = 3;
else
    figure(colormapFigHandle);
    goceDensitySubplot = 2;
    msisDensitySubplot = 4;
end
secondsInDay = 60 * 60 * 24;
[minLat, maxLat] = findInterpolationLimits(magneticLatitude);
indicesToRemove = findMatrixIndicesInDatagap(regriddedTime, timestamps10s);
regriddedTime(indicesToRemove) = nan(1);
regriddedLatitude(indicesToRemove) = nan(1);
regriddedGoceDensity(indicesToRemove) = nan(1);
regriddedMsisDensity(indicesToRemove) = nan(1);

regriddedTime = regriddedTime / secondsInDay + firstDatenum;
referenceDay = datestr(min(regriddedTime(:)), 'mmmm dd, yyyy');
regriddedTime = regriddedTime - datenum(referenceDay, 'mmmm dd, yyyy');

subplot(2,2,goceDensitySubplot)
surf(regriddedTime, regriddedLatitude, regriddedGoceDensity, 'EdgeColor', 'None')

xlim([min(regriddedTime(:)) max(regriddedTime(:))]);
ylim([minLat maxLat]);
caxis([min(regriddedGoceDensity(:)) max(regriddedGoceDensity(:))])
colorbar('Location', 'NorthOutside');
view(2);
ylabel('Geomagnetic latitude (°)')
title(['Goce ', timeOfDay,' density'])

subplot(2,2,msisDensitySubplot)
surf(regriddedTime, regriddedLatitude, regriddedMsisDensity, 'EdgeColor', 'None')

xlim([min(regriddedTime(:)) max(regriddedTime(:))]);
ylim([minLat maxLat]);
caxis([min(regriddedGoceDensity(:)) max(regriddedGoceDensity(:))])
colorbar('Location', 'NorthOutside');
view(2);
xlabel(['Days since the beginning of ', referenceDay])
ylabel('Geomagnetic latitude (°)')
title(['Msis ', timeOfDay,' density'])

% if ~isempty(strfind(lower(timeOfDay), 'evening'))
%     tightfig(colormapFigHandle);
% end

end