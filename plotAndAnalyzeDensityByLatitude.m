function [results] = plotAndAnalyzeDensityByLatitude(firstDatenum, ae, timestamps1min, aeIntegral, timestampsAeInt, timestamps1minFixed, correctedDensity, msisDensity, jbDensity, dtmDensity, tiegcmDensity, aeProxyDensity,...
    timestamps10s, magneticLatitude, i, winds, timeOfDay, plotFigures, results)
% plotCorrectedDensityLatitudes(ae, timestamps1min, correctedDensity, timestamps10s, latitude, timestampsDatenum, computeLatitudes);

if ~isempty(strfind(lower(timeOfDay), 'morning'))
    goceV = winds.morningGoceV{i};
    goceU = winds.morningGoceU{i};
    hwmV = winds.morningHwmV{i};
    hwmU = winds.morningHwmV{i};
else
    goceV = winds.eveningGoceV{i};
    goceU = winds.eveningGoceU{i};
    hwmV = winds.eveningHwmV{i};
    hwmU = winds.eveningHwmV{i};
end

[limitedTimestamps, limitedLatitude, minAllowedLatitude, maxAllowedLatitude] = giveExactOrbits(timestamps10s, magneticLatitude);

[crossingTimes, goceDensityByLatitude, msisDensityByLatitude, jbDensityByLatitude, dtmDensityByLatitude, aeProxyDensityByLatitude, tiegcmDensityByLatitude, oneDegreeStep] = interpolateAndPlotByLatitude(firstDatenum, aeIntegral, timestampsAeInt, timestamps10s, magneticLatitude, ...
    correctedDensity, msisDensity, jbDensity, dtmDensity, tiegcmDensity, aeProxyDensity, goceV, goceU, hwmV, hwmU, limitedLatitude, limitedTimestamps, minAllowedLatitude, maxAllowedLatitude, plotFigures, timeOfDay);

if plotFigures ~= 0
    analyzeLagByLatitude(timestamps1min, ae, goceDensityByLatitude, crossingTimes, minAllowedLatitude, maxAllowedLatitude, timeOfDay, 'DTM');
end

results = plotAndAnalyzeByHemisphere(firstDatenum, goceDensityByLatitude, ae, crossingTimes, timestamps1min, timeOfDay, oneDegreeStep, plotFigures, 'DTM', results);

end

function plotDensityLatitudeTimeSurf(firstDatenum, aeIntegral, timestamps1min, magneticLatitude, timestamps10s, regriddedLatitude, regriddedTime, regriddedGoceDensity, regriddedMsisDensity, regriddedJbDensity, regriddedDtmDensity, regriddedTiegcmDensity, ...
    regriddedAeProxy,  goceVMatrix, goceUMatrix, hwmVMatrix, hwmUMatrix, timeOfDay)
% plotDensityLatitudeTimeSurf(averagedDensity, averagedLatitude, timestamps

persistent colormapFigHandle
persistent windFigHandle
persistent minDensity
persistent maxDensity

numPlotRows = 2; % 2 or 6
if ~isempty(strfind(lower(timeOfDay), 'morning'))
    colormapFigHandle = figure('Color', 'white', 'units','normalized','outerposition',[0 0 1 1]);
    windFigHandle = figure('Color', 'white', 'units','normalized','outerposition',[0 0 1 1]);
    goceSubplot = 1;
    tiegcmSubplot = 3;
    hwmSubplot = 3;
    if numPlotRows == 6
        msisDensitySubplot = 5;
        jbDensitySubplot = 7;
        dtmDensitySubplot = 9;
        aeProxySubplot = 11;
    end
    minDensity = min(regriddedGoceDensity(:));
    maxDensity = max(regriddedGoceDensity(:));
else
    goceSubplot = 2;
    tiegcmSubplot = 4;
    hwmSubplot = 4;
    if numPlotRows == 6
        msisDensitySubplot = 6;
        jbDensitySubplot = 8;
        dtmDensitySubplot = 10;
        aeProxySubplot = 12;
    end
end
secondsInDay = 60 * 60 * 24;
minDensityTime = min(regriddedTime(:));
[minLat, maxLat] = findInterpolationLimits(magneticLatitude);
indicesToRemove = findMatrixIndicesInDatagap(regriddedTime, timestamps10s) | (regriddedTime < secondsInDay + minDensityTime);
regriddedTime(indicesToRemove) = nan(1);

regriddedTime = regriddedTime / secondsInDay + firstDatenum;
referenceDay = datestr(min(regriddedTime(:)), 'mmmm dd, yyyy');
regriddedTime = regriddedTime - datenum(referenceDay, 'mmmm dd, yyyy');

timestamps1min = timestamps1min / secondsInDay + firstDatenum - datenum(referenceDay, 'mmmm dd, yyyy');

minDensityTime = min(regriddedTime(:));
maxDensityTime = max(regriddedTime(:));
firstDayIndices = regriddedTime < minDensityTime + 1;
firstDayGoce = regriddedGoceDensity(firstDayIndices);
firstDayMsis = regriddedMsisDensity(firstDayIndices);
firstDayAe = regriddedAeProxy(firstDayIndices);
firstDayJb = regriddedJbDensity(firstDayIndices);
firstDayDtm = regriddedDtmDensity(firstDayIndices);
firstDayTiegcm = regriddedTiegcmDensity(firstDayIndices);
msisMultiplier = mean(firstDayGoce(:) ./ firstDayMsis(:));
aeMultiplier = mean(firstDayGoce(:) ./ firstDayAe(:));
jbMultiplier = mean(firstDayGoce(:) ./ firstDayJb(:));
dtmMultiplier = mean(firstDayGoce(:) ./ firstDayDtm(:));
tiegcmMultiplier = mean(firstDayGoce(:) ./ firstDayTiegcm(:));

plotIndices = (regriddedTime >= minDensityTime + 0.0 & regriddedTime <= maxDensityTime);
[plotRows, ~] = find(plotIndices);
plotRows = unique(plotRows);

regriddedGoceDensity = regriddedGoceDensity(plotRows, :);
regriddedTime = regriddedTime(plotRows, :);
regriddedLatitude = regriddedLatitude(plotRows, :);
regriddedMsisDensity = msisMultiplier * regriddedMsisDensity(plotRows, :);
regriddedAeProxy = aeMultiplier * regriddedAeProxy(plotRows, :);
regriddedJbDensity = jbMultiplier * regriddedJbDensity(plotRows, :);
regriddedDtmDensity = dtmMultiplier * regriddedDtmDensity(plotRows, :);
regriddedTiegcmDensity = regriddedTiegcmDensity(plotRows, :);%tiegcmMultiplier * regriddedTiegcmDensity(plotRows, :);
goceVMatrix = goceVMatrix(plotRows, :);
hwmVMatrix = hwmVMatrix(plotRows, :);
goceUMatrix = goceUMatrix(plotRows, :);
hwmUMatrix = hwmUMatrix(plotRows, :);
% i = round(size(regriddedTime, 1) / 20);
% j = round(size(regriddedTime, 2) / 15);
% goceVMatrix = goceVMatrix(1:i:end, 1:j:end);
% goceUMatrix = goceUMatrix(1:i:end, 1:j:end);
% hwmVMatrix = hwmVMatrix(1:i:end, 1:j:end);
% hwmUMatrix = hwmUMatrix(1:i:end, 1:j:end);

plotHeight = maxDensity;

minDensityTime = minDensityTime + 0.0;
maxDensityTime = maxDensityTime - 0.0;
[~, minAeIndex] = min(abs(timestamps1min - minDensityTime));
[~, maxAeIndex] = min(abs(timestamps1min - maxDensityTime));
aeIndicesToPlot = minAeIndex:maxAeIndex;

figure(colormapFigHandle);

subplotAxesHandle = subplot(numPlotRows,2,goceSubplot);
surf(subplotAxesHandle, regriddedTime, regriddedLatitude, regriddedGoceDensity, 'EdgeColor', 'None')
xlim([minDensityTime maxDensityTime]);
ylim([minLat maxLat]);
caxis([minDensity maxDensity])
view(2);
colorbar('Location', 'EastOutside');
ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
title(['Goce ', timeOfDay,' density'], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
set(gca, 'fontsize', 12)

hold all;
aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
set(aeLineHandle, 'LineWidth', 0.1)
view(2);
set(aeAxesHandle, 'yaxislocation', 'right');
ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
set(aeAxesHandle, 'Color', 'none', 'XTick', []);
hold off;
set(gca, 'fontsize', 12)

subplotAxesHandle = subplot(numPlotRows,2,tiegcmSubplot);
surf(subplotAxesHandle, regriddedTime, regriddedLatitude, regriddedTiegcmDensity, 'EdgeColor', 'None')
xlim([minDensityTime maxDensityTime]);
ylim([minLat maxLat]);
%caxis([minDensity maxDensity])
view(2);
colorbar('Location', 'EastOutside');
ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
title(['TIEGCM ', timeOfDay,' density'], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
set(gca, 'fontsize', 12)

hold all;
aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
set(aeLineHandle, 'LineWidth', 0.1)
view(2);
set(aeAxesHandle, 'yaxislocation', 'right');
ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
set(aeAxesHandle, 'Color', 'none', 'XTick', []);
hold off;
set(gca, 'fontsize', 12)

if numPlotRows == 6
    subplotAxesHandle = subplot(numPlotRows,2,msisDensitySubplot);
    surf(subplotAxesHandle, regriddedTime, regriddedLatitude, regriddedMsisDensity, 'EdgeColor', 'None')

    xlim([minDensityTime maxDensityTime]);
    ylim([minLat maxLat]);
    caxis([minDensity maxDensity])
    colorbar('Location', 'EastOutside');
    view(2);
    ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
    title(['NRLMSISE-00 ', timeOfDay,' density'], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
    set(gca, 'fontsize', 12)

    hold all;
    aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
    aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
    set(aeLineHandle, 'LineWidth', 0.1)
    view(2);
    set(aeAxesHandle, 'yaxislocation', 'right');
    ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
    set(aeAxesHandle, 'Color', 'none', 'XTick', []);
    hold off;
    set(gca, 'fontsize', 12)

    subplotAxesHandle = subplot(numPlotRows,2,jbDensitySubplot);
    surf(subplotAxesHandle, regriddedTime, regriddedLatitude, regriddedJbDensity, 'EdgeColor', 'None')

    xlim([minDensityTime maxDensityTime]);
    ylim([minLat maxLat]);
    caxis([minDensity maxDensity])
    colorbar('Location', 'EastOutside');
    view(2);
    ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
    title(['JB2008 ', timeOfDay,' density'], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
    set(gca, 'fontsize', 12)

    hold all;
    aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
    aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
    set(aeLineHandle, 'LineWidth', 0.1)
    view(2);
    set(aeAxesHandle, 'yaxislocation', 'right');
    ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
    set(aeAxesHandle, 'Color', 'none', 'XTick', []);
    hold off;
    set(gca, 'fontsize', 12)


    subplotAxesHandle = subplot(numPlotRows,2,dtmDensitySubplot);
    surf(subplotAxesHandle, regriddedTime, regriddedLatitude, regriddedDtmDensity, 'EdgeColor', 'None')

    xlim([minDensityTime maxDensityTime]);
    ylim([minLat maxLat]);
    caxis([minDensity maxDensity])
    colorbar('Location', 'EastOutside');
    view(2);
    ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
    title(['DTM-2013 ', timeOfDay,' density'], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
    set(gca, 'fontsize', 12)

    hold all;
    aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
    aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
    set(aeLineHandle, 'LineWidth', 0.1)
    view(2);
    set(aeAxesHandle, 'yaxislocation', 'right');
    ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
    set(aeAxesHandle, 'Color', 'none', 'XTick', []);
    hold off;
    set(gca, 'fontsize', 12)


    subplotAxesHandle = subplot(numPlotRows,2,aeProxySubplot);
    surf(subplotAxesHandle, regriddedTime, regriddedLatitude, regriddedAeProxy, 'EdgeColor', 'None')

    xlim([minDensityTime maxDensityTime]);
    ylim([minLat maxLat]);
    caxis([minDensity maxDensity])
    colorbar('Location', 'EastOutside');
    view(2);
    xlabel(['Days since the UTC beginning of ', referenceDay], 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
    ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
    title(['AE Predicted ', timeOfDay,' density'], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
    set(gca, 'fontsize', 12)

    hold all;
    aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
    aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
    set(aeLineHandle, 'LineWidth', 0.1)
    view(2);
    set(aeAxesHandle, 'yaxislocation', 'right');
    ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
    set(aeAxesHandle, 'Color', 'none', 'XTick', []);
    hold off;
    set(gca, 'fontsize', 12)
end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

figure(windFigHandle);

goceTotWind = sqrt(goceVMatrix.^2 + goceUMatrix.^2);
upDir = repmat([0 0 1], size(goceTotWind(:)), 1);
northDir = repmat([0 1 0], size(goceTotWind(:)), 1);
goceCrWind3d = [goceUMatrix(:), goceVMatrix(:), zeros(size(goceTotWind(:)))];
ind = (dot(upDir, cross(northDir, goceCrWind3d, 2), 2)) > 0;
goceTotWind(ind) = -goceTotWind(ind);

hwmTotWind = sqrt(hwmVMatrix.^2 + hwmUMatrix.^2);
hwmCrWind3d = [hwmUMatrix(:), hwmVMatrix(:), zeros(size(hwmTotWind(:)))];
ind = (dot(upDir, cross(northDir, hwmCrWind3d, 2), 2)) > 0;
hwmTotWind(ind) = -hwmTotWind(ind);

plotHeight = max(goceTotWind(:));
greatestVel = max([max(goceVMatrix(:)), -min(goceVMatrix(:))]);

subplotAxesHandle = subplot(2,2,goceSubplot);

% %quiver(subplotAxesHandle, regriddedTime(1:i:end,1:j:end), regriddedLatitude(1:i:end,1:j:end), goceUMatrix, goceVMatrix, 'o')
% quiver(subplotAxesHandle, goceUMatrix, goceVMatrix, 'o')
% % xlim([minDensityTime maxDensityTime]);
% % ylim([minLat maxLat]);
% ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
% title(['Goce wind ', timeOfDay], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
% set(gca, 'fontsize', 12)

surf(subplotAxesHandle, regriddedTime, regriddedLatitude, goceVMatrix, 'EdgeColor', 'None')
xlim([minDensityTime maxDensityTime]);
ylim([minLat maxLat]);
colorbar('Location', 'EastOutside');
caxis([-greatestVel greatestVel])
view(2);
ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
title(['Goce wind ', timeOfDay], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
set(gca, 'fontsize', 12)

% hold all;
% aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
% aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
% set(aeLineHandle, 'LineWidth', 0.1)
% view(2);
% set(aeAxesHandle, 'yaxislocation', 'right');
% ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
% set(aeAxesHandle, 'Color', 'none', 'XTick', []);
% hold off;
% set(gca, 'fontsize', 12)



subplotAxesHandle = subplot(2,2,hwmSubplot);

% %quiver(subplotAxesHandle, regriddedTime(1:i:end,1:j:end), regriddedLatitude(1:i:end,1:j:end), hwmUMatrix, hwmVMatrix, 'o')
% quiver(subplotAxesHandle, hwmUMatrix, hwmVMatrix, 'o')
% % xlim([minDensityTime maxDensityTime]);
% % ylim([minLat maxLat]);
% ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
% title(['HWM wind ', timeOfDay], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
% set(gca, 'fontsize', 12)

surf(subplotAxesHandle, regriddedTime, regriddedLatitude, hwmVMatrix, 'EdgeColor', 'None')
xlim([minDensityTime maxDensityTime]);
ylim([minLat maxLat]);
colorbar('Location', 'EastOutside');
caxis([-greatestVel greatestVel])
view(2);
ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
title(['HWM07 wind ', timeOfDay], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
set(gca, 'fontsize', 12)
% 
% hold all;
% aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
% aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
% set(aeLineHandle, 'LineWidth', 0.1)
% view(2);
% set(aeAxesHandle, 'yaxislocation', 'right');
% ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
% xlabel(['Days since the UTC beginning of ', referenceDay], 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
% set(aeAxesHandle, 'Color', 'none', 'XTick', []);
% hold off;
% set(gca, 'fontsize', 12)

end

function [crossingTimes, goceDensityByLatitude, msisDensityByLatitude, jbDensityByLatitude, dtmDensityByLatitude, aeProxyDensityByLatitude, tiegcmDensityByLatitude, oneDegreeStep] = interpolateAndPlotByLatitude(firstDatenum, aeIntegral, timestamps1min, timestamps10s, magneticLatitude, ...
    correctedDensity, msisDensity, jbDensity, dtmDensity, tiegcmDensity, aeProxyDensity, goceV, goceU, hwmV, hwmU, limitedLatitude, limitedTimestamps, minAllowedLatitude, maxAllowedLatitude, plotFigures, timeOfDay)
%

oneDegreeStep = minAllowedLatitude:maxAllowedLatitude;
if plotFigures == 0   
    goceInterpolant = scatteredInterpolant(timestamps10s, magneticLatitude, correctedDensity);
    msisInterpolant = scatteredInterpolant(timestamps10s, magneticLatitude, msisDensity);
    jbInterpolant = scatteredInterpolant(timestamps10s, magneticLatitude, jbDensity);
    dtmInterpolant = scatteredInterpolant(timestamps10s, magneticLatitude, dtmDensity);
    aeProxyInterpolant = scatteredInterpolant(timestamps10s, magneticLatitude, aeProxyDensity);
    tiegcmInterpolant = scatteredInterpolant(timestamps10s, magneticLatitude, tiegcmDensity);
    for i = 1:length(oneDegreeStep)
        crossingTimes(:,i) = latitudeCrossingTimes(limitedLatitude, limitedTimestamps, oneDegreeStep(i));    
        goceDensityByLatitude(:,i) = goceInterpolant(crossingTimes(:,i), ones(size(crossingTimes(:,i))) * oneDegreeStep(i));
        msisDensityByLatitude(:,i) = msisInterpolant(crossingTimes(:,i), ones(size(crossingTimes(:,i))) * oneDegreeStep(i));
        jbDensityByLatitude(:,i) = jbInterpolant(crossingTimes(:,i), ones(size(crossingTimes(:,i))) * oneDegreeStep(i));
        dtmDensityByLatitude(:,i) = dtmInterpolant(crossingTimes(:,i), ones(size(crossingTimes(:,i))) * oneDegreeStep(i));
        aeProxyDensityByLatitude(:,i) = aeProxyInterpolant(crossingTimes(:,i), ones(size(crossingTimes(:,i))) * oneDegreeStep(i));
        tiegcmDensityByLatitude(:,i) = tiegcmInterpolant(crossingTimes(:,i), ones(size(crossingTimes(:,i))) * oneDegreeStep(i));
    end
else
    oneQuarterDegreeStep = minAllowedLatitude:0.25:maxAllowedLatitude;
    for i = 1:length(oneQuarterDegreeStep)
        regriddedTime(:,i) = latitudeCrossingTimes(limitedLatitude, limitedTimestamps, oneQuarterDegreeStep(i));    
    end
    
%     goceV = smooth(goceV, 15);
%     goceU = smooth(goceU, 15);
%     hwmV = smooth(hwmV, 15);
%     hwmU = smooth(hwmU, 15);
    
    regriddedGoceDensity = interp1(timestamps10s, correctedDensity, regriddedTime, 'spline');
    regriddedMsisDensity = interp1(timestamps10s, msisDensity, regriddedTime, 'spline');
    regriddedJbDensity = interp1(timestamps10s, jbDensity, regriddedTime, 'spline');
    regriddedDtmDensity = interp1(timestamps10s, dtmDensity, regriddedTime, 'spline');
    regriddedTiegcmDensity = interp1(timestamps10s, tiegcmDensity, regriddedTime, 'spline');
    regriddedAeProxy = interp1(timestamps10s, aeProxyDensity, regriddedTime, 'spline');
    regriddedGoceV = interp1(timestamps10s, goceV, regriddedTime, 'spline');
    regriddedHwmV = interp1(timestamps10s, hwmV, regriddedTime, 'spline');
    regriddedGoceU = interp1(timestamps10s, goceU, regriddedTime, 'spline');
    regriddedHwmU = interp1(timestamps10s, hwmU, regriddedTime, 'spline');
    
    crossingTimes = regriddedTime(:,1:4:end);
    goceDensityByLatitude = regriddedGoceDensity(:,1:4:end);
    msisDensityByLatitude = regriddedMsisDensity(:,1:4:end);
    jbDensityByLatitude = regriddedJbDensity(:,1:4:end);
    tiegcmDensityByLatitude = regriddedTiegcmDensity(:,1:4:end);
    dtmDensityByLatitude = regriddedDtmDensity(:,1:4:end);
    aeProxyDensityByLatitude = regriddedAeProxy(:,1:4:end);
    
    
    numOfOrbits = length(regriddedTime(:,1));
    numOfValuesInOrbit = length(regriddedTime(1,:));
    for i = 1:numOfValuesInOrbit
        timeThisLatitude = regriddedTime(:,i);
        goceDensityThisLatitude = regriddedGoceDensity(:,i);
        msisDensityThisLatitude = regriddedMsisDensity(:,i);
        jbDensityThisLatitude = regriddedJbDensity(:,i);
        dtmDensityThisLatitude = regriddedDtmDensity(:,i);
        tiegcmDensityThisLatitude = regriddedTiegcmDensity(:,i);
        aeProxyThisLatitude = regriddedAeProxy(:,i);
        goceVThisLatitude = regriddedGoceV(:,i);
        hwmVThisLatitude = regriddedHwmV(:,i);
        goceUThisLatitude = regriddedGoceU(:,i);
        hwmUThisLatitude = regriddedHwmU(:,i);

        tInterp = interp1(1:numOfOrbits, timeThisLatitude, 1:1/20:numOfOrbits);
        interpolatedGoceDensity = interp1(timeThisLatitude, goceDensityThisLatitude, tInterp, 'spline');
        interpolatedMsisDensity = interp1(timeThisLatitude, msisDensityThisLatitude, tInterp, 'spline');
        interpolatedJbDensity = interp1(timeThisLatitude, jbDensityThisLatitude, tInterp, 'spline');
        interpolatedDtmDensity = interp1(timeThisLatitude, dtmDensityThisLatitude, tInterp, 'spline');
        interpolatedTiegcmDensity = interp1(timeThisLatitude, tiegcmDensityThisLatitude, tInterp, 'spline');
        interpolatedAeProxy = interp1(timeThisLatitude, aeProxyThisLatitude, tInterp, 'spline');
        interpolatedGoceV = interp1(timeThisLatitude, goceVThisLatitude, tInterp, 'spline');
        interpolatedHwmV = interp1(timeThisLatitude, hwmVThisLatitude, tInterp, 'spline'); 
        interpolatedGoceU = interp1(timeThisLatitude, goceUThisLatitude, tInterp, 'spline');
        interpolatedHwmU = interp1(timeThisLatitude, hwmUThisLatitude, tInterp, 'spline');

        latitudeMatrix(:,i) = ones(length(tInterp), 1) * oneQuarterDegreeStep(i); 
        goceDensityMatrix(:,i) = interpolatedGoceDensity;
        msisDensityMatrix(:,i) = interpolatedMsisDensity;
        jbDensityMatrix(:,i) = interpolatedJbDensity;
        dtmDensityMatrix(:,i) = interpolatedDtmDensity;
        tiegcmDensityMatrix(:,i) = interpolatedTiegcmDensity;
        aeProxyDensityMatrix(:,i) = interpolatedAeProxy;
        goceVMatrix(:,i) = interpolatedGoceV;
        hwmVMatrix(:,i) = interpolatedHwmV;
        goceUMatrix(:,i) = interpolatedGoceU;
        hwmUMatrix(:,i) = interpolatedHwmU;
        timeMatrix(:,i) = tInterp;
    end

    plotDensityLatitudeTimeSurf(firstDatenum, aeIntegral, timestamps1min, magneticLatitude, timestamps10s, latitudeMatrix, ...
        timeMatrix, goceDensityMatrix, msisDensityMatrix, jbDensityMatrix, dtmDensityMatrix, tiegcmDensityMatrix, aeProxyDensityMatrix, goceVMatrix, goceUMatrix, hwmVMatrix, hwmUMatrix, timeOfDay);
end

end

function analyzeLagByLatitude(timestamps1min, ae, densityByLatitude, crossingTimes, minAllowedLatitude, maxAllowedLatitude, timeOfDay, densityType)
% writeAndPlotPeakAnalysis(timestamps, ae, splinedDensity, computeLatitutes)

parfor i = 1:length(minAllowedLatitude:maxAllowedLatitude)
    interpolatedDensity(:,i) = interp1(crossingTimes(:,i), densityByLatitude(:,i), timestamps1min, 'nearest', 'extrap');
end

[timelag, errInLag, allLatitudesAeLag, northernInterval, step] = analyzeNorthernLatitudes(minAllowedLatitude, maxAllowedLatitude, ...
    interpolatedDensity, ae);

[timelag, errInLag, intervalLatitudes] = analyzeSouthernLatitudes(minAllowedLatitude, maxAllowedLatitude, ...
    interpolatedDensity, ae, timelag, errInLag, allLatitudesAeLag, step, northernInterval);
     
plotPeakAnalysisFigure(timelag, errInLag, intervalLatitudes, timeOfDay, densityType)

end

function plotPeakAnalysisFigure(timelag, errInLag, intervalLatitudes, timeOfDay, densityType)
%

persistent aeGoceDensityLagByLatHandle
persistent aeMsisDensityLagByLatHandle

if strcmpi(densityType, 'goce'); thisFig = aeGoceDensityLagByLatHandle; else thisFig = aeMsisDensityLagByLatHandle; end

if ~isempty(strfind(lower(timeOfDay), 'morning')); thisFig = figure; subplotNum = 1; else subplotNum = 2; end
figure(thisFig);
subplot(1,2,subplotNum)
timestampsInHours = timelag / 60;
errInHours = errInLag / 60;
herrorbar(timestampsInHours, intervalLatitudes, errInHours, '-s')
xlabel('t / h');
ylabel('Geomagnetic latitude');
title(timeOfDay)
annotation('textbox', [0 0.9 1 0.1], ...
    'String', [densityType, ' Density-AE lag for different latitudes'], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

if strcmpi(densityType, 'goce'); aeGoceDensityLagByLatHandle = thisFig; else aeMsisDensityLagByLatHandle = thisFig; end

end

function [timelag, errInLag, allLatitudesAeLag, northernInterval, step] = analyzeNorthernLatitudes(minAllowedLatitude, maxAllowedLatitude, interpolatedDensity, ae)
%

step = 8;
northernInterval = -step/2:step:maxAllowedLatitude;
northernInterval(ismember(northernInterval, maxAllowedLatitude)) = [];
allLatitudesAeLag = zeros(maxAllowedLatitude - minAllowedLatitude + 1, 1);
catTimelag = 0;
for i = northernInterval
    densityInterval = i - minAllowedLatitude + 1 : i - minAllowedLatitude + 1 + step;
    if ~isempty(find(densityInterval > maxAllowedLatitude - minAllowedLatitude + 1, 1))
        densityInterval = i - minAllowedLatitude + 1 : maxAllowedLatitude - minAllowedLatitude + 1;
        catTimelag = 1;
    end
    
    lagsAllIntervalLats = zeros(length(densityInterval), 1);
    for k = densityInterval
        maxCrossCorr = giveMaxCrossCorrelation(interpolatedDensity(:,k), ae, 60);
        allLatitudesAeLag(k) = maxCrossCorr;
        lagsAllIntervalLats(k - min(densityInterval) + 1) = maxCrossCorr;
    end
    
    if catTimelag > 0
        timelag = [timelag mean(lagsAllIntervalLats)];
        errInLag = [errInLag (std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats)))];  
    else
        timelag(length(i:-step:minAllowedLatitude) + 1) = mean(lagsAllIntervalLats);
        errInLag(length(i:-step:minAllowedLatitude) + 1) = (std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats)));
    end
end

end

function [timelag, errInLag, intervalLatitudes] = analyzeSouthernLatitudes(minAllowedLatitude, maxAllowedLatitude, ...
    interpolatedDensity, ae, timelag, errInLag, allLatitudesAeLag, step, northernInterval)
%

southernInterval = -step/2:-step:minAllowedLatitude;
southernInterval(ismember(southernInterval, minAllowedLatitude)) = [];
for i = southernInterval
    densityInterval = i - minAllowedLatitude + 1 - step : i - minAllowedLatitude + 1;
    if ~isempty(find(densityInterval < 1, 1))
        densityInterval = 1 : i - minAllowedLatitude + 1;
    end
    
    lagsAllIntervalLats = zeros(length(densityInterval), 1);
    for k = densityInterval
        maxCrossCorr = giveMaxCrossCorrelation(interpolatedDensity(:,k), ae, 60);
        allLatitudesAeLag(k) = maxCrossCorr;
        lagsAllIntervalLats(k - min(densityInterval) + 1) = maxCrossCorr;
    end
    
    timelag(length(i:-step:minAllowedLatitude)) = mean(lagsAllIntervalLats);
    errInLag(length(i:-step:minAllowedLatitude)) = std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats));
end

intervalLatitudes = (southernInterval(end) : step : northernInterval(end)) + step/2;
intervalLatitudes(end) = [];
intervalLatitudes = [ round(mean([southernInterval(end) minAllowedLatitude]))...
                      intervalLatitudes                                      ...
                      round(mean([northernInterval(end) maxAllowedLatitude])) ]; 

end

function results = plotAndAnalyzeByHemisphere(firstDatenum, densityByLatitude, ae, crossingTimes, timestamps1min, timeOfDay, oneDegreeStep, plotFigures, densityType, results)
%

northIndices = (oneDegreeStep < 90 & oneDegreeStep > 45);
equatorIndices = (oneDegreeStep < 30 & oneDegreeStep > -30);
southIndices = (oneDegreeStep < -45 & oneDegreeStep > -90);

northernDensity = mean(densityByLatitude(:,northIndices), 2);
equatorDensity = mean(densityByLatitude(:,equatorIndices), 2);
southernDensity = mean(densityByLatitude(:,southIndices), 2);

northernError = std(densityByLatitude(:,northIndices), 0, 2);
equatorError = std(densityByLatitude(:,equatorIndices), 0, 2);
southernError = std(densityByLatitude(:,southIndices), 0, 2);

northTimestamps = mean(crossingTimes(:,northIndices), 2);
equatorTimestamps = mean(crossingTimes(:,equatorIndices), 2);
southTimestamps = mean(crossingTimes(:,southIndices), 2);

[~, northPeakBegin, northPeakEnd] = limitToNearPeak(northernDensity, 'noSmooth', 'median');
[~, equatorPeakBegin, equatorPeakEnd] = limitToNearPeak(equatorDensity, 'noSmooth', 'median');
[~, southPeakBegin, southPeakEnd] = limitToNearPeak(southernDensity, 'noSmooth', 'median');

northCalmMean = mean(northernDensity(1:northPeakBegin));
equatorCalmMean = mean(equatorDensity(1:equatorPeakBegin));
southCalmMean = mean(southernDensity(1:southPeakBegin));

nanIndices = isnan(southCalmMean);
if ~isempty(find(nanIndices, 1))
    a = 1;
end

northAbsDiff = northernDensity - northCalmMean;
equatorAbsDiff = equatorDensity - equatorCalmMean;
southAbsDiff = southernDensity - southCalmMean;

northRelDiff = northAbsDiff / northCalmMean;
equatorRelDiff = equatorAbsDiff / equatorCalmMean;
southRelDiff = southAbsDiff / southCalmMean;

meanNorthAbsDiff = mean(northAbsDiff(northPeakBegin:northPeakEnd));
meanEquatorAbsDiff = mean(equatorAbsDiff(equatorPeakBegin:equatorPeakEnd));
meanSouthAbsDiff = mean(southAbsDiff(southPeakBegin:southPeakEnd));
meanNorthRelDiff = mean(northRelDiff(northPeakBegin:northPeakEnd));
meanEquatorRelDiff = mean(equatorRelDiff(equatorPeakBegin:equatorPeakEnd));
meanSouthRelDiff = mean(southRelDiff(southPeakBegin:southPeakEnd));

maxNorthAbsDiff = max(northAbsDiff);
maxEquatorAbsDiff = max(equatorAbsDiff);
maxSouthAbsDiff = max(southAbsDiff);
maxNorthRelDiff = max(northRelDiff);
maxEquatorRelDiff = max(equatorRelDiff);
maxSouthRelDiff = max(southRelDiff);

% [northEfoldingTimeRisingLimb, equatorEfoldingTimeRisingLimb, southEfoldingTimeRisingLimb, northEfoldingTimeFallingLimb, equatorEfoldingTimeFallingLimb,...
%  southEfoldingTimeFallingLimb, northRisingBegin, equatorRisingBegin, southRisingBegin] = giveEfoldingTimes(northTimestamps, equatorTimestamps, southTimestamps, northernDensity, equatorDensity, southernDensity,timestamps1min);

[northEfoldingTimeRisingLimb, northEfoldingTimeFallingLimb, northRisingBegin] = giveEfoldingTimes(northTimestamps, northernDensity, timestamps1min);
[equatorEfoldingTimeRisingLimb, equatorEfoldingTimeFallingLimb, equatorRisingBegin] = giveEfoldingTimes(equatorTimestamps, equatorDensity, timestamps1min);
[southEfoldingTimeRisingLimb, southEfoldingTimeFallingLimb, southRisingBegin] = giveEfoldingTimes(southTimestamps, southernDensity, timestamps1min);

northDensForXcorr = interp1(northTimestamps, northernDensity, timestamps1min, 'linear', mean([northernDensity(1), northernDensity(end)]));
equatorDensForXcorr = interp1(equatorTimestamps, equatorDensity, timestamps1min, 'linear', mean([equatorDensity(1), equatorDensity(end)]));
southDensForXcorr = interp1(southTimestamps, southernDensity, timestamps1min, 'linear', mean([southernDensity(1), southernDensity(end)]));

northDensForXcorr = limitToDay(northDensForXcorr, timestamps1min, northRisingBegin);
equatorDensForXcorr = limitToDay(equatorDensForXcorr, timestamps1min, equatorRisingBegin);
southDensForXcorr = limitToDay(southDensForXcorr, timestamps1min, southRisingBegin);
northDensForXcorr = deleteDescendingParts(northDensForXcorr);
equatorDensForXcorr = deleteDescendingParts(equatorDensForXcorr);
southDensForXcorr = deleteDescendingParts(southDensForXcorr);

aeSmooth = deleteDescendingParts(smooth(ae, 91));
aeNorth = limitToDay(aeSmooth, timestamps1min, northRisingBegin);
northLag = giveMaxCrossCorrelation(northDensForXcorr, aeNorth, 60);
aeEq = limitToDay(aeSmooth, timestamps1min, equatorRisingBegin);
equatorLag = giveMaxCrossCorrelation(equatorDensForXcorr, aeEq, 60);
aeSouth = limitToDay(aeSmooth, timestamps1min, southRisingBegin);
southLag = giveMaxCrossCorrelation(southDensForXcorr, aeSouth, 60);

% if northLag < eps
%     t = limitToDay(timestamps1min, timestamps1min, northRisingBegin);
%     t = (t-t(1))/3600;
%     plotyy(t, northDensForXcorr, t, aeNorth)
% end

[rowNum, ~] = size(results);
emptyCells = cellfun(@isempty,results);
[~, emptyColPositions] = find(emptyCells);
colNum = min(emptyColPositions);
if rowNum == 2
    colNum = length(results(rowNum,:)) + 1;
    results{1, colNum}     = [densityType, ' NH AE-Dens. ', timeOfDay, ' lag'];
    results{1, colNum + 1} = [densityType, ' EQ AE-Dens. ', timeOfDay, ' lag'];
    results{1, colNum + 2} = [densityType, ' SH AE-Dens. ', timeOfDay, ' lag'];
    
    results{1, colNum + 3} = [densityType, ' Mean NH Abs. Diff ', timeOfDay];
    results{1, colNum + 4} = [densityType, ' Mean EQ Abs. Diff ', timeOfDay];
    results{1, colNum + 5} = [densityType, ' Mean SH Abs. Diff ', timeOfDay];
    results{1, colNum + 6} = [densityType, ' Mean NH Rel. Diff ', timeOfDay];
    results{1, colNum + 7} = [densityType, ' Mean EQ Rel. Diff ', timeOfDay];
    results{1, colNum + 8} = [densityType, ' Mean SH Rel. Diff ', timeOfDay];
    
    results{1, colNum + 9} = [densityType, ' Max NH Abs. Diff ', timeOfDay];
    results{1, colNum + 10} = [densityType, ' Max EQ Abs. Diff ', timeOfDay];
    results{1, colNum + 11} = [densityType, ' Max SH Abs. Diff ', timeOfDay];
    results{1, colNum + 12} = [densityType, ' Max NH Rel. Diff ', timeOfDay];
    results{1, colNum + 13} = [densityType, ' Max EQ Rel. Diff ', timeOfDay];
    results{1, colNum + 14} = [densityType, ' Max SH Rel. Diff ', timeOfDay];
    
    results{1, colNum + 15} = [densityType, ' NH e-fold.t. rising limb ', timeOfDay];
    results{1, colNum + 16} = [densityType, ' EQ e-fold.t. rising limb ', timeOfDay];
    results{1, colNum + 17} = [densityType, ' SH e-fold.t. rising limb ', timeOfDay];
    results{1, colNum + 18} = [densityType, ' NH e-fold.t. falling limb ', timeOfDay];
    results{1, colNum + 19} = [densityType, ' EQ e-fold.t. falling limb ', timeOfDay];
    results{1, colNum + 20} = [densityType, ' SH e-fold.t. falling limb ', timeOfDay];
end
results{rowNum, colNum} =     northLag;
results{rowNum, colNum + 1} = equatorLag;
results{rowNum, colNum + 2} = southLag;

results{rowNum, colNum + 3} = meanNorthAbsDiff;
results{rowNum, colNum + 4} = meanEquatorAbsDiff;
results{rowNum, colNum + 5} = meanSouthAbsDiff;
results{rowNum, colNum + 6} = meanNorthRelDiff;
results{rowNum, colNum + 7} = meanEquatorRelDiff;
results{rowNum, colNum + 8} = meanSouthRelDiff;

results{rowNum, colNum + 9} = maxNorthAbsDiff;
results{rowNum, colNum + 10} = maxEquatorAbsDiff;
results{rowNum, colNum + 11} = maxSouthAbsDiff;
results{rowNum, colNum + 12} = maxNorthRelDiff;
results{rowNum, colNum + 13} = maxEquatorRelDiff;
results{rowNum, colNum + 14} = maxSouthRelDiff;

results{rowNum, colNum + 15} = northEfoldingTimeRisingLimb;
results{rowNum, colNum + 16} = equatorEfoldingTimeRisingLimb;
results{rowNum, colNum + 17} = southEfoldingTimeRisingLimb;
results{rowNum, colNum + 18} = northEfoldingTimeFallingLimb;
results{rowNum, colNum + 19} = equatorEfoldingTimeFallingLimb;
results{rowNum, colNum + 20} = southEfoldingTimeFallingLimb;

if plotFigures ~= 0
    secondsInDay = 24 * 60 * 60;

    northTimestamps = northTimestamps / secondsInDay + firstDatenum;
    equatorTimestamps = equatorTimestamps / secondsInDay + firstDatenum;
    southTimestamps = southTimestamps / secondsInDay + firstDatenum;
    
    figHandle = figure;
    hAx(1) = subplot(3,2,1);
    errorbar(northTimestamps, northAbsDiff, northernError);
    datetick('x', 'dd', 'keeplimits')
    ylabel('Density [10^{-11} kgm^{-3}]')
    title('Absolute difference, NH')
    grid on

    hAx(3) = subplot(3,2,3);
    errorbar(equatorTimestamps, equatorAbsDiff, equatorError, 'g');
    datetick('x', 'dd', 'keeplimits')
    ylabel('Density [10^{-11} kgm^{-3}]')
    title('Absolute difference, equator')
    grid on

    hAx(5) = subplot(3,2,5);
    errorbar(southTimestamps, southAbsDiff, southernError, 'r');
    datetick('x', 'dd', 'keeplimits')
    ylabel('Density [10^{-11} kgm^{-3}]')
    title('Absolute difference, SH')
    grid on

    hAx(2) = subplot(3,2,2);
    errorbar(northTimestamps, northRelDiff * 100, 100 * northernError/northCalmMean);
    datetick('x', 'dd', 'keeplimits')
    ylabel('%')
    title('Relative difference, NH')
    grid on

    hAx(4) = subplot(3,2,4);
    errorbar(equatorTimestamps, equatorRelDiff * 100, 100 * equatorError/equatorCalmMean, 'g');
    datetick('x', 'dd', 'keeplimits')
    ylabel('%')
    title('Relative difference, equator')
    grid on

    hAx(6) = subplot(3,2,6);
    errorbar(southTimestamps, southRelDiff * 100, 100 * southernError/southCalmMean, 'r');
    datetick('x', 'dd', 'keeplimits')
    ylabel('%')
    title('Relative difference, SH')
    grid on
    
    axis(hAx, 'tight')

    annotation('textbox', [0 0.9 1 0.1], ...
    'String', [densityType, ' absolute and relative density changes on hemispheres: ', timeOfDay], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

    %tightfig(figHandle);
end

end

function value = limitToDay(value, timestamps, ind)

[~,dayBeforeInd] = min(abs(timestamps - timestamps(ind) + 86400));
[~,dayAfterInd] = min(abs(timestamps - timestamps(ind) - 2*86400));
value = value(dayBeforeInd:dayAfterInd);

end

function value = deleteDescendingParts(value)

dVal = diff(value);
dVal = [dVal; dVal(end)];
value(dVal < 0) = 0;

end

function [densityIndexTimelag] = giveMaxCrossCorrelation(density, geomIndex, indicesInHour)
% plotCrossCorrelation(averagedDensityNoBg, ae)

maxLag = round(indicesInHour * 10.5);
[correlations, lags] = xcorr(density, geomIndex, maxLag, 'coeff');
correlations = correlations(lags >= 0);
lags = lags(lags >= 0);
densityIndexTimelag = lags(correlations == max(correlations)) / indicesInHour;

end


function hh = herrorbar(x, y, l, u, symbol)
%HERRORBAR Horizontal Error bar plot.
%   HERRORBAR(X,Y,L,R) plots the graph of vector X vs. vector Y with
%   horizontal error bars specified by the vectors L and R. L and R contain the
%   left and right error ranges for each point in X. Each error bar
%   is L(i) + R(i) long and is drawn a distance of L(i) to the right and R(i)
%   to the right the points in (X,Y). The vectors X,Y,L and R must all be
%   the same length. If X,Y,L and R are matrices then each column
%   produces a separate line.
%
%   HERRORBAR(X,Y,E) or HERRORBAR(Y,E) plots X with error bars [X-E X+E].
%   HERRORBAR(...,'LineSpec') uses the color and linestyle specified by
%   the string 'LineSpec'. See PLOT for possibilities.
%
%   H = HERRORBAR(...) returns a vector of line handles.
%
%   Example:
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      herrorbar(x,y,e)
%   draws symmetric horizontal error bars of unit standard deviation.
%
%   This code is based on ERRORBAR provided in MATLAB.   
%
%   See also ERRORBAR

% Copyright (c) 2009, Jos van der Geest
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%   Jos van der Geest
%   email: jos@jasen.nl
%
%   File history:
%   August 2006 (Jos): I have taken back ownership. I like to thank Greg Aloe from
%   The MathWorks who originally introduced this piece of code to the
%   Matlab File Exchange. 
%   September 2003 (Greg Aloe): This code was originally provided by Jos
%   from the newsgroup comp.soft-sys.matlab:
%   http://newsreader.mathworks.com/WebX?50@118.fdnxaJz9btF^1@.eea3ff9
%   After unsuccessfully attempting to contact the orignal author, I
%   decided to take ownership so that others could benefit from finding it
%   on the MATLAB Central File Exchange.

if min(size(x))==1,
    npt = length(x);
    x = x(:);
    y = y(:);
    if nargin > 2,
        if ~isstr(l),
            l = l(:);
        end
        if nargin > 3
            if ~isstr(u)
                u = u(:);
            end
        end
    end
else
    [npt,n] = size(x);
end

if nargin == 3
    if ~isstr(l)
        u = l;
        symbol = '-';
    else
        symbol = l;
        l = y;
        u = y;
        y = x;
        [m,n] = size(y);
        x(:) = (1:npt)'*ones(1,n);;
    end
end

if nargin == 4
    if isstr(u),
        symbol = u;
        u = l;
    else
        symbol = '-';
    end
end

if nargin == 2
    l = y;
    u = y;
    y = x;
    [m,n] = size(y);
    x(:) = (1:npt)'*ones(1,n);;
    symbol = '-';
end

u = abs(u);
l = abs(l);

if isstr(x) | isstr(y) | isstr(u) | isstr(l)
    error('Arguments must be numeric.')
end

if ~isequal(size(x),size(y)) | ~isequal(size(x),size(l)) | ~isequal(size(x),size(u)),
    error('The sizes of X, Y, L and U must be the same.');
end

tee = (max(y(:))-min(y(:)))/100; % make tee .02 x-distance for error bars
% changed from errorbar.m
xl = x - l;
xr = x + u;
ytop = y + tee;
ybot = y - tee;
n = size(y,2);
% end change

% Plot graph and bars
hold_state = ishold;
cax = newplot;
next = lower(get(cax,'NextPlot'));

% build up nan-separated vector for bars
% changed from errorbar.m
xb = zeros(npt*9,n);
xb(1:9:end,:) = xl;
xb(2:9:end,:) = xl;
xb(3:9:end,:) = NaN;
xb(4:9:end,:) = xl;
xb(5:9:end,:) = xr;
xb(6:9:end,:) = NaN;
xb(7:9:end,:) = xr;
xb(8:9:end,:) = xr;
xb(9:9:end,:) = NaN;

yb = zeros(npt*9,n);
yb(1:9:end,:) = ytop;
yb(2:9:end,:) = ybot;
yb(3:9:end,:) = NaN;
yb(4:9:end,:) = y;
yb(5:9:end,:) = y;
yb(6:9:end,:) = NaN;
yb(7:9:end,:) = ytop;
yb(8:9:end,:) = ybot;
yb(9:9:end,:) = NaN;
% end change


[ls,col,mark,msg] = colstyle(symbol); if ~isempty(msg), error(msg); end
symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid

h = plot(xb,yb,esymbol); hold on
h = [h;plot(x,y,symbol)];

if ~hold_state, hold off; end

if nargout>0, hh = h; end

end
