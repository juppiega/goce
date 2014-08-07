function [results] = plotAndAnalyzeDensityByLatitude(firstDatenum, ae, timestamps1min, timestamps1minFixed, correctedDensity, msisDensity,...
    timestamps10s, magneticLatitude, timeOfDay, plotFigures, results)
% plotCorrectedDensityLatitudes(ae, timestamps1min, correctedDensity, timestamps10s, latitude, timestampsDatenum, computeLatitudes);
    
[limitedTimestamps, limitedLatitude, minAllowedLatitude, maxAllowedLatitude] = giveExactOrbits(timestamps10s, magneticLatitude);

[crossingTimes, densityByLatitude, oneDegreeStep] = interpolateAndPlotByLatitude(firstDatenum, ae, timestamps1min, timestamps10s, magneticLatitude, ...
    correctedDensity, msisDensity, limitedLatitude, limitedTimestamps, minAllowedLatitude, maxAllowedLatitude, plotFigures, timeOfDay);

if plotFigures ~= 0
    analyzeLagByLatitude(limitedTimestamps, timestamps1min, timestamps1minFixed, ae, densityByLatitude, crossingTimes, minAllowedLatitude, maxAllowedLatitude, timeOfDay);
end

results = plotAndAnalyzeByHemisphere(firstDatenum, densityByLatitude, ae, crossingTimes, timestamps1min, timeOfDay, oneDegreeStep, plotFigures, results);

end

function plotDensityLatitudeTimeSurf(firstDatenum, ae, timestamps1min, magneticLatitude, timestamps10s, regriddedLatitude, regriddedTime, regriddedGoceDensity, regriddedMsisDensity, timeOfDay)
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

timestamps1min = timestamps1min / secondsInDay + firstDatenum - datenum(referenceDay, 'mmmm dd, yyyy');

minDensityTime = min(regriddedTime(:));
maxDensityTime = max(regriddedTime(:));

minDensity = min(regriddedGoceDensity(:));
maxDensity = max(regriddedGoceDensity(:));

[~, minAeIndex] = min(abs(timestamps1min - minDensityTime));
[~, maxAeIndex] = min(abs(timestamps1min - maxDensityTime));
aeIndicesToPlot = minAeIndex:maxAeIndex;

hAx = subplot(2,2,goceDensitySubplot);
h = plot3(timestamps1min(aeIndicesToPlot), ae(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * maxDensity, 'k');
view(2);
set(hAx, 'yaxislocation', 'right');
ylabel(hAx, 'AE')
hold all;
hAxNew = axis();
surf(regriddedTime, regriddedLatitude, regriddedGoceDensity, 'EdgeColor', 'None')

xlim([minDensityTime maxDensityTime]);
ylim([minLat maxLat]);
caxis([minDensity maxDensity])
colorbar('Location', 'NorthOutside');
ylabel('Geomagnetic latitude (°)')
xlabel(['Days since the UTC beginning of ', referenceDay])
title(['Goce ', timeOfDay,' density'])

% hold on;
% hAx = axis();
% ishandle(hAx)
% h = plot(hAx, timestamps1min(aeIndicesToPlot), ae(aeIndicesToPlot), 'k');% ones(size(aeIndicesToPlot)) * maxDensity,
% set(hAx, 'yaxislocation', 'right');
% ylabel(hAx, 'AE')

% hold all;
% h = plot3(timestamps1min(aeIndicesToPlot), ae(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * maxDensity);
% set(h, 'yaxislocation', 'right');
% hold off;

subplot(2,2,msisDensitySubplot)
surf(regriddedTime, regriddedLatitude, regriddedMsisDensity, 'EdgeColor', 'None')

xlim([minDensityTime maxDensityTime]);
ylim([minLat maxLat]);
caxis([minDensity maxDensity])
colorbar('Location', 'NorthOutside');
view(2);
xlabel(['Days since the UTC beginning of ', referenceDay])
ylabel('Geomagnetic latitude (°)')
title(['Msis ', timeOfDay,' density'])

% hold all;
% h = plot3(timestamps1min(aeIndicesToPlot), ae(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * maxDensity);
% set(h, 'yaxislocation', 'right');
% hold off;

% if ~isempty(strfind(lower(timeOfDay), 'evening'))
%     tightfig(colormapFigHandle);
% end

end

function [limitedTimestamps, limitedLatitude, minAllowedLatitude, maxAllowedLatitude] = giveExactOrbits(timestamps10s, magneticLatitude)
%

[latBeginIndex, latEndIndex] = limitLatitudeToIntegerMultipleOfOrbitalPeriod(magneticLatitude);
limitedLatitude = magneticLatitude(latBeginIndex:latEndIndex);
limitedTimestamps = timestamps10s(latBeginIndex:latEndIndex);
[~, exactOrbitIndices] = splitIntoOrbits(limitedLatitude, limitedTimestamps);
limitedLatitude = limitedLatitude(exactOrbitIndices);
limitedTimestamps = limitedTimestamps(exactOrbitIndices);
[orbits, ~] = splitIntoOrbits(limitedLatitude, limitedTimestamps);

[minAllowedLatitude, maxAllowedLatitude] = findInterpolationLimits(limitedLatitude);
orbitsToDelete = find(abs(limitedLatitude(orbits(:,2)) - limitedLatitude(orbits(:,1))) < ...
    (maxAllowedLatitude - minAllowedLatitude));
newIndices = 1:length(limitedLatitude);
for i = 1:length(orbitsToDelete)
    newIndices = setdiff(newIndices, (orbits(orbitsToDelete(i), 1) : orbits(orbitsToDelete(i), 2)));
end
limitedTimestamps = limitedTimestamps(newIndices);
limitedLatitude = limitedLatitude(newIndices);

end

function [crossingTimes, densityByLatitude, oneDegreeStep] = interpolateAndPlotByLatitude(firstDatenum, ae, timestamps1min, timestamps10s, magneticLatitude, ...
    correctedDensity, msisDensity, limitedLatitude, limitedTimestamps, minAllowedLatitude, maxAllowedLatitude, plotFigures, timeOfDay)
%

oneDegreeStep = minAllowedLatitude:maxAllowedLatitude;
if plotFigures == 0   
    F = scatteredInterpolant(timestamps10s, magneticLatitude, correctedDensity);
    for i = 1:length(oneDegreeStep)
        crossingTimes(:,i) = latitudeCrossingTimes(limitedLatitude, limitedTimestamps, oneDegreeStep(i));    
        densityByLatitude(:,i) = F(crossingTimes(:,i), ones(size(crossingTimes(:,i))) * oneDegreeStep(i));
    end
else
    oneQuarterDegreeStep = minAllowedLatitude:0.25:maxAllowedLatitude;
    for i = 1:length(oneQuarterDegreeStep)
        regriddedTime(:,i) = latitudeCrossingTimes(limitedLatitude, limitedTimestamps, oneQuarterDegreeStep(i));    
    end

    regriddedGoceDensity = interp1(timestamps10s, correctedDensity, regriddedTime, 'spline');
    regriddedMsisDensity = interp1(timestamps10s, msisDensity, regriddedTime, 'spline');
    crossingTimes = regriddedTime(:,1:4:end);
    densityByLatitude = regriddedGoceDensity(:,1:4:end);
    
    numOfOrbits = length(regriddedTime(:,1));
    numOfValuesInOrbit = length(regriddedTime(1,:));
    for i = 1:numOfValuesInOrbit
        timeThisLatitude = regriddedTime(:,i);
        goceDensityThisLatitude = regriddedGoceDensity(:,i);
        msisDensityThisLatitude = regriddedMsisDensity(:,i);

        tInterp = interp1(1:numOfOrbits, timeThisLatitude, 1:1/20:numOfOrbits);
        interpolatedGoceDensity = interp1(timeThisLatitude, goceDensityThisLatitude, tInterp, 'spline');
        interpolatedMsisDensity = interp1(timeThisLatitude, msisDensityThisLatitude, tInterp, 'spline');

        latitudeMatrix(:,i) = ones(length(tInterp), 1) * oneQuarterDegreeStep(i); 
        goceDensityMatrix(:,i) = interpolatedGoceDensity;
        msisDensityMatrix(:,i) = interpolatedMsisDensity;
        timeMatrix(:,i) = tInterp;
    end

    plotDensityLatitudeTimeSurf(firstDatenum, ae, timestamps1min, magneticLatitude, timestamps10s, latitudeMatrix, timeMatrix, goceDensityMatrix, msisDensityMatrix, timeOfDay);
end

end

function analyzeLagByLatitude(limitedTimestamps, timestamps1min, timestamps1minFixed, ae, densityByLatitude, crossingTimes, minAllowedLatitude, maxAllowedLatitude, timeOfDay)
% writeAndPlotPeakAnalysis(timestamps, ae, splinedDensity, computeLatitutes)

% limitedTimestamps = limitedTimestamps(ismember(limitedTimestamps, timestamps1minFixed));
% ae = ae(ismember(timestamps1min, limitedTimestamps));
% ae = removePeriodicBackground(ae, 432, 1, 0); % 432 min =^ 0.3 days
% ae = smooth(ae, 431);
% [ae, ~, ~] = limitToNearPeak(ae, 'noSmooth', 'mean');
% for i = 1:length(densityByLatitude(1,:))
%    samplingFreq = 86400 / mean(diff(crossingTimes(:,i)));
%    densityByLatitude(:,i) = removePeriodicBackground(densityByLatitude(:,i), 0.3, samplingFreq, 0);
% end

parfor i = 1:length(minAllowedLatitude:maxAllowedLatitude)
    interpolatedDensity(:,i) = interp1(crossingTimes(:,i), densityByLatitude(:,i), timestamps1min, 'nearest', 'extrap');
end

[timelag, errInLag, allLatitudesAeLag, northernInterval, step] = analyzeNorthernLatitudes(minAllowedLatitude, maxAllowedLatitude, ...
    interpolatedDensity, ae);

[timelag, errInLag, intervalLatitudes] = analyzeSouthernLatitudes(minAllowedLatitude, maxAllowedLatitude, ...
    interpolatedDensity, ae, timelag, errInLag, allLatitudesAeLag, step, northernInterval);
     
plotPeakAnalysisFigure(timelag, errInLag, intervalLatitudes, timeOfDay)

end

function plotPeakAnalysisFigure(timelag, errInLag, intervalLatitudes, timeOfDay)
%

persistent aeDensityLagByLatHandle

if ~isempty(strfind(lower(timeOfDay), 'morning')); aeDensityLagByLatHandle = figure; subplotNum = 1; else subplotNum = 2; end
figure(aeDensityLagByLatHandle);
subplot(1,2,subplotNum)
timestampsInHours = timelag / 60;
errInHours = errInLag / 60;
herrorbar(timestampsInHours, intervalLatitudes, errInHours, '-s')
xlabel('t / h');
ylabel('Geomagnetic latitude');
title(timeOfDay)
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Density-AE lag for different latitudes', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

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
        maxCrossCorr = giveMaxCrossCorrelation(interpolatedDensity(:,k), ae);
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
        maxCrossCorr = giveMaxCrossCorrelation(interpolatedDensity(:,k), ae);
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

function results = plotAndAnalyzeByHemisphere(firstDatenum, densityByLatitude, ae, crossingTimes, timestamps1min, timeOfDay, oneDegreeStep, plotFigures, results)
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

northRelDiffForXcorr = interp1(northTimestamps, northRelDiff, timestamps1min, 'nearest', 'extrap');
equatorRelDiffForXcorr = interp1(equatorTimestamps, equatorRelDiff, timestamps1min, 'nearest', 'extrap');
southRelDiffForXcorr = interp1(southTimestamps, southRelDiff, timestamps1min, 'nearest', 'extrap');

northLag = giveMaxCrossCorrelation(northRelDiffForXcorr, ae);
equatorLag = giveMaxCrossCorrelation(equatorRelDiffForXcorr, ae);
southLag = giveMaxCrossCorrelation(southRelDiffForXcorr, ae);

[rowNum, ~] = size(results);
emptyCells = cellfun(@isempty,results);
[~, emptyColPositions] = find(emptyCells);
colNum = min(emptyColPositions);
if rowNum == 2
    colNum = length(results(rowNum,:)) + 1;
    results{1, colNum}     = ['NH AE-Dens. ', timeOfDay, ' lag'];
    results{1, colNum + 1} = ['EQ AE-Dens. ', timeOfDay, ' lag'];
    results{1, colNum + 2} = ['SH AE-Dens. ', timeOfDay, ' lag'];
    
    results{1, colNum + 3} = ['Mean NH Abs. Diff ', timeOfDay];
    results{1, colNum + 4} = ['Mean EQ Abs. Diff ', timeOfDay];
    results{1, colNum + 5} = ['Mean SH Abs. Diff ', timeOfDay];
    results{1, colNum + 6} = ['Mean NH Rel. Diff ', timeOfDay];
    results{1, colNum + 7} = ['Mean EQ Rel. Diff ', timeOfDay];
    results{1, colNum + 8} = ['Mean SH Rel. Diff ', timeOfDay];
    
    results{1, colNum + 9} = ['Max NH Abs. Diff ', timeOfDay];
    results{1, colNum + 10} = ['Max EQ Abs. Diff ', timeOfDay];
    results{1, colNum + 11} = ['Max SH Abs. Diff ', timeOfDay];
    results{1, colNum + 12} = ['Max NH Rel. Diff ', timeOfDay];
    results{1, colNum + 13} = ['Max EQ Rel. Diff ', timeOfDay];
    results{1, colNum + 14} = ['Max SH Rel. Diff ', timeOfDay];
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
    'String', ['Absolute and relative density changes on hemispheres: ', timeOfDay], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

    tightfig(figHandle);
end

end

function [densityIndexTimelag] = giveMaxCrossCorrelation(density, geomIndex)
% plotCrossCorrelation(averagedDensityNoBg, ae)

maxLag = 60 * 24;
[correlations, lags] = xcorr(density, geomIndex, maxLag, 'coeff');
correlations = correlations(lags > 0);
lags = lags(lags > 0);
indicesInHour = 60;
densityIndexTimelag = lags(correlations == max(correlations)) / indicesInHour;

end
