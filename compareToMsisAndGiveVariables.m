function [ae, ap, absB, vBz, akasofuEpsilon, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningAeProxy, eveningAeProxy, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, ...
 morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength, firstDatenum] ...
 = compareToMsisAndGiveVariables(threshold, results)

if exist('goceVariables.mat', 'file') == 2
    load('goceVariables.mat')
else
    fprintf(2, '\n%s\n', 'Warning: no goceVariables.mat found in Matlab PATH, now attempting to create new one -> readFiles.m');
    readFiles();
    load('goceVariables.mat')
end

%compareGoceDensityToMsis(measuredDensity, msisDensityVariableAlt, ae, timestampsAeDatenum, timestampsDensityDatenum, results);

intervalsOfInterest = findInterestingIntervals(ae, timestampsAeDatenum, timestamps1minFixed, averagedDensityNoBg, epsilonQualityFlag, timestampsEpsilonDatenum, timestampsDensityDatenum, threshold);

[morningTimestamps10s, morningMagneticLatitude, morningDensityNoBg, morningMsisDensity, eveningTimestamps10s, ...
    morningAeProxy, eveningAeProxy, eveningMagneticLatitude, morningLatitudeForFit, eveningLatitudeForFit, eveningDensityNoBg, eveningMsisDensity, morningIndex, eveningIndex] = ...
    splitBySolarTime(timestamps10sFixed, magneticLatitude, latitude, densityNoBg, msisDensity270km, aeProxyDensity, densityIndex, solarTime);
morningTimestamps10sAll = morningTimestamps10s;
eveningTimestamps10sAll = eveningTimestamps10s;

%plotParametrizationResults(morningFourierGrid, eveningFourierGrid, morningBins, eveningBins);

[ae, ap, absB, vBz, akasofuEpsilon, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningAeProxy, eveningAeProxy, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = sliceToInterestingIntervals(ae, ap, absB, vBz, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningAeProxy, eveningAeProxy, morningTimestamps10s, eveningTimestamps10s, timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed,...
 density3h, morningMagneticLatitude, eveningMagneticLatitude, timestamps10sFixed, timestamps1min, timestampsAbsB, akasofuEpsilon, timestampsDensityDatenum, intervalsOfInterest);

[morningGrid, morningBins] = computeTimeCells(timestamps10sFixed, morningTimestamps10s, latitude);
[eveningGrid, eveningBins] = computeTimeCells(timestamps10sFixed, eveningTimestamps10s, latitude);

[predictedStormDensity, morningFtest, eveningFtest, morningFits, eveningFits] = computeParametrization(morningGrid, morningBins, eveningGrid, eveningBins, aeIntegrals, morningFourierGrid, eveningFourierGrid, timestamps10sFixed,...
    morningLatitudeForFit, eveningLatitudeForFit, morningTimestamps10sAll, eveningTimestamps10sAll, doy, msisDensity270kmNoAp, densityIndex);

save('stormfits.mat', 'morningFtest', '-v7.3');
save('stormfits.mat', 'eveningFtest', '-append');
save('stormfits.mat', 'morningFits', '-append');
save('stormfits.mat', 'eveningFits', '-append');
save('stormfits.mat', 'predictedStormDensity', '-append');

end


function [morningTimestamps10s, morningMagneticLatitude, morningDensityNoBg, morningMsisDensity, eveningTimestamps10s, ...
    morningAeProxy, eveningAeProxy, eveningMagneticLatitude, morningLatitude, eveningLatitude, eveningDensityNoBg, eveningMsisDensity, morningIndex, eveningIndex] = ...
    splitBySolarTime(timestamps10s, magneticLatitude, latitude, densityNoBg, msisDensity, aeProxy, densityIndex, solarTime)
%

morningIndices = find(solarTime <= 12);
eveningIndices = find(solarTime > 12);

morningTimestamps10s = timestamps10s(morningIndices);
morningMagneticLatitude = magneticLatitude(morningIndices);
morningLatitude = latitude(morningIndices);
morningDensityNoBg = densityNoBg(morningIndices);
morningMsisDensity = msisDensity(morningIndices);
morningAeProxy = aeProxy(morningIndices);
morningIndex = densityIndex(morningIndices);

eveningTimestamps10s = timestamps10s(eveningIndices);
eveningMagneticLatitude = magneticLatitude(eveningIndices);
eveningLatitude = latitude(eveningIndices);
eveningDensityNoBg = densityNoBg(eveningIndices);
eveningMsisDensity = msisDensity(eveningIndices);
eveningAeProxy = aeProxy(eveningIndices);
eveningIndex = densityIndex(eveningIndices);

end


function intervalsOfInterest = findInterestingIntervals(ae, timestampsAeDatenum, timestamps1minFixed, averagedDensityNoBg, epsilonQualityFlag, timestampsEpsilonDatenum, timestampsDensityDatenum, threshold)
%

[intervalsOfInterest, calmDays] = findAllStormsAboveThreshold(ae, timestampsAeDatenum, threshold);

intervalsOfInterest = connectStormsWithMuchOverlap(intervalsOfInterest);

intervalsOfInterest = removeStormsWithDatagaps(intervalsOfInterest, ae, timestampsAeDatenum, ...
    timestampsEpsilonDatenum, timestampsDensityDatenum, epsilonQualityFlag, calmDays);

printSelectedStorms(intervalsOfInterest, timestampsAeDatenum);

end

function [intervalsOfInterest, calmDays] = findAllStormsAboveThreshold(ae, timestampsAeDatenum, threshold)
%

intervalsOfInterest = zeros(1, 2);
medianCrossings = findCrossings(ae, 'smooth', 'mean');
calmDays = 2;
if length(medianCrossings) > 1
    for i = 1:length(medianCrossings) - 1
        if ~isempty(find(ae(medianCrossings(i):medianCrossings(i + 1)) >= threshold, 1))
            [intervalBegin, intervalEnd] = addCalmDaysToBothSidesOfPeak(medianCrossings(i), medianCrossings(i + 1), timestampsAeDatenum, calmDays);
            intervalsOfInterest = vertcat(intervalsOfInterest, [intervalBegin intervalEnd]);
        end
    end
end

intervalsOfInterest(1,:) = [];
if isempty(intervalsOfInterest)
    fprintf(2, '%s %d %s %d\n', 'There were no AE values above threshold: ', threshold, '. The maximum of given AE data is ', max(ae))
    error('No AE values above threshold found!')
end

end

function printSelectedStorms(intervalsOfInterest, timestampsAeDatenum)
%

fprintf('\n%s\n', 'Following storms will be analyzed:')
for i = 1:length(intervalsOfInterest(:,1))
    beginIndex = intervalsOfInterest(i,1);
    endIndex = intervalsOfInterest(i,2);
    fprintf('%d%s\n', i, ['. ', datestr(timestampsAeDatenum(beginIndex), 'yyyy-mm-dd'), ' to ', datestr(timestampsAeDatenum(endIndex), 'yyyy-mm-dd')])
end

end

function intervalsOfInterest = removeStormsWithDatagaps(intervalsOfInterest, ae, timestampsAeDatenum, ...
    timestampsEpsilonDatenum, timestampsDensityDatenum, epsilonQualityFlag, calmDays)
%

fprintf('\n%s\n', '')

indicesInDay = 24 * 60;
intervalsToRemove = [];
for i = 1:length(intervalsOfInterest(:,1))
    aeIndices = intervalsOfInterest(i,1):intervalsOfInterest(i,2);
    [~, densityIntervalBeginIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(aeIndices(1))));
    [~, densityIntervalEndIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(aeIndices(end))));
    [~, epsilonIntervalBeginIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(aeIndices(1))));
    [~, epsilonIntervalEndIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(aeIndices(end))));
    epsilonIntervalIndices = epsilonIntervalBeginIndex:epsilonIntervalEndIndex;
    densityIntervalIndices = densityIntervalBeginIndex:densityIntervalEndIndex;
    
    stormBeginIndex = aeIndices(1) + calmDays * indicesInDay;
    stormEndIndex = aeIndices(end) - calmDays * indicesInDay;
    [~, densityStormBeginIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(stormBeginIndex)));
    [~, densityStormEndIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(stormEndIndex)));
    [~, epsilonStormBeginIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(stormBeginIndex)));
    [~, epsilonStormEndIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(stormEndIndex)));
    epsilonStormIndices = length(find(epsilonQualityFlag(epsilonStormBeginIndex:epsilonStormEndIndex)));
    densityStormIndices = length(find(ismember(timestampsDensityDatenum(densityStormBeginIndex:densityStormEndIndex), timestampsAeDatenum(stormBeginIndex:stormEndIndex))));
    epsilonBeforeIndices = length(find(epsilonQualityFlag(epsilonIntervalBeginIndex:epsilonStormBeginIndex)));
    densityBeforeIndices = length(find(ismember(timestampsDensityDatenum(densityIntervalBeginIndex:densityStormBeginIndex), timestampsAeDatenum(aeIndices(1):stormBeginIndex))));
    epsilonAfterIndices = length(find(epsilonQualityFlag(epsilonStormEndIndex:epsilonIntervalEndIndex)));
    densityAfterIndices = length(find(ismember(timestampsDensityDatenum(densityStormEndIndex:densityIntervalEndIndex), timestampsAeDatenum(stormEndIndex:aeIndices(end)))));
    
    valuesShouldHaveAfter = length(timestampsAeDatenum(stormEndIndex:aeIndices(end)));
    valuesShouldHaveBefore = length(timestampsAeDatenum(aeIndices(1):stormBeginIndex));
    valuesShouldHaveStorm = length(timestampsAeDatenum(stormBeginIndex:stormEndIndex));
    
    if densityStormIndices / valuesShouldHaveStorm < 0.75 ||...
       (densityAfterIndices / valuesShouldHaveAfter < 0.75 && densityBeforeIndices / valuesShouldHaveBefore < 0.75)
        intervalsToRemove = [intervalsToRemove i];
        fprintf(2, '%s%d%s\n', ['Warning: Storm between dates ', datestr(timestampsAeDatenum(aeIndices(1)), 'yyyy-mm-dd'), ...
            ' and ', datestr(timestampsAeDatenum(aeIndices(end)), 'yyyy-mm-dd'), ' with max AE: '], max(ae(aeIndices)), ...
            ' has too large density data gaps. It will be omitted.');
    end
    
    if epsilonStormIndices / valuesShouldHaveStorm < 0.5 ||...
       (epsilonAfterIndices / valuesShouldHaveAfter < 0.3 && epsilonBeforeIndices / valuesShouldHaveBefore < 0.3)
        intervalsToRemove = [intervalsToRemove i];
        fprintf(2, '%s%d%s\n', ['Warning: Storm between dates ', datestr(timestampsAeDatenum(aeIndices(1)), 'yyyy-mm-dd'), ...
            ' and ', datestr(timestampsAeDatenum(aeIndices(end)), 'yyyy-mm-dd'), ' with max AE: '], max(ae(aeIndices)), ...
            ' has too large solar wind data gaps. It will be omitted.');
    end
      
end

intervalsToConserve = setdiff(1:length(intervalsOfInterest(:,1)), intervalsToRemove);
intervalsOfInterest = intervalsOfInterest(intervalsToConserve,:);

if isempty(intervalsOfInterest)
    fprintf(2, '\n%s \n', 'All storms over threshold contained too large data gaps. Decrease the threshold!')
    error('Too large data gaps during all requested storms')
end

end

function intervalsOfInterest = connectStormsWithMuchOverlap(intervalsOfInterest)
%

if length(intervalsOfInterest(:,1)) > 1
    previousOverlappingStorms = 0;
    newIntervals = intervalsOfInterest;
    indicesToRemove = [];
    for i = 1:length(intervalsOfInterest(:,1)) - 1
        nextStormOverlap = intervalsOfInterest(i, 2) - intervalsOfInterest(i + 1,1);
        nextStormLength = intervalsOfInterest(i + 1, 2) - intervalsOfInterest(i + 1,1);
        if nextStormOverlap / nextStormLength > 0.4
            if previousOverlappingStorms > 0
                newIntervals(i - previousOverlappingStorms,:) = [newIntervals(i - previousOverlappingStorms, 1) intervalsOfInterest(i + 1, 2)];
            else
                newIntervals(i,:) = [intervalsOfInterest(i, 1) intervalsOfInterest(i + 1, 2)];
            end
            indicesToRemove = [indicesToRemove; (i + 1)];
            previousOverlappingStorms = previousOverlappingStorms + 1;
        else
            previousOverlappingStorms = 0;
        end
    end
    indicesToConserve = setdiff(1:length(newIntervals(:,1)), indicesToRemove);
    intervalsOfInterest = newIntervals(indicesToConserve,:);
end

end

function [intervalBegin, intervalEnd] = addCalmDaysToBothSidesOfPeak(peakBeginIndex, peakEndIndex, timestampsAeDatenum, calmDays)
%

firstAeDay = timestampsAeDatenum(1);
lastAeDay = timestampsAeDatenum(end);
peakBegin = timestampsAeDatenum(peakBeginIndex);
peakEnd = timestampsAeDatenum(peakEndIndex);

if floor(peakBegin - calmDays) < firstAeDay
    intervalBegin = firstAeDay;
else
    intervalBegin = floor(peakBegin - calmDays);
end

if ceil(peakEnd + calmDays) > lastAeDay
    intervalEnd = floor(lastAeDay);
else
    intervalEnd = ceil(peakEnd + calmDays);
end

intervalBegin = find(timestampsAeDatenum == intervalBegin);
intervalEnd = find(timestampsAeDatenum == intervalEnd);

end


function printCorrResults(corrResults)
%

celldisp(corrResults);

end

function [predictedDensity, morningFtest, eveningFtest, morningFits, eveningFits] = computeParametrization(morningGrid, morningBins, eveningGrid, eveningBins, aeIntegrals, morningFourier, eveningFourier, timestamps10sFixed,...
    morningLatitude, eveningLatitude, morningTimestamps10s, eveningTimestamps10s, doy, msis270kmNoAp, densityIndex)
%

fprintf('%s\n', 'Computing parametrizations. This may take up to two to three hours.')

aeIntFixed = aeIntegrals ./ mean(aeIntegrals(:));

targetCount = length(morningBins) + length(eveningBins);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running Parametrizations, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

morningProxy = zeros(length(timestamps10sFixed), length(morningBins));
morningFits = cell(length(morningBins), 1);

parfor i = 1:length(morningBins)
    timeThisLat = morningGrid{i};
    thisLatIndices = ismember(timestamps10sFixed, timeThisLat);
    densityThisLat = densityIndex(thisLatIndices);
    fourierFit = morningFourier{i};
    
    fullMatrix = [fourierFit(doy) aeIntFixed];
    x5x7 = fullMatrix(:,5) .* fullMatrix(:,7);
    x5x8 = fullMatrix(:,5) .* fullMatrix(:,8);
    x9x10 = fullMatrix(:,9) .* fullMatrix(:,10);
    x7square = fullMatrix(:,7) .^2;
    fullMatrix = [fullMatrix(:,[1 2 3 4 5 9]) x5x7 x5x8 x9x10 x7square];
    
    proxyMatrix = fullMatrix(thisLatIndices,:);
    linearModel = fitlm(proxyMatrix, densityThisLat);
    
    finalProxyDensity = feval(linearModel, fullMatrix);
    morningFits{i} = linearModel;
    
    resultArray = table2array(anova(linearModel, 'component'));
    morningFtest(:,i) = resultArray(:,4);     
    
    morningProxy(:,i) = finalProxyDensity; %interp1(timeThisLat, finalProxyDensity, timestamps10sFixed, 'linear', 0);
    
    p.progress;
end

eveningProxy = zeros(length(timestamps10sFixed), length(eveningBins));
eveningFits = cell(length(eveningBins), 1);
parfor i = 1:length(eveningBins)
    timeThisLat = eveningGrid{i};
    thisLatIndices = ismember(timestamps10sFixed, timeThisLat);
    densityThisLat = densityIndex(thisLatIndices);
    fourierFit = eveningFourier{i};
    
    fullMatrix = [fourierFit(doy) aeIntFixed];
    x5x7 = fullMatrix(:,5) .* fullMatrix(:,7);
    x5x8 = fullMatrix(:,5) .* fullMatrix(:,8);
    x9x10 = fullMatrix(:,9) .* fullMatrix(:,10);
    x7square = fullMatrix(:,7) .^2;
    fullMatrix = [fullMatrix(:,[1 2 3 4 5 9]) x5x7 x5x8 x9x10 x7square];

    proxyMatrix = fullMatrix(thisLatIndices,:);
    linearModel = fitlm(proxyMatrix, densityThisLat);
    
    finalProxyDensity = feval(linearModel, fullMatrix);
    eveningFits{i} = linearModel;
    
    resultArray = table2array(anova(linearModel, 'component'));
    eveningFtest(:,i) = resultArray(:,4);    
    
    eveningProxy(:,i) = finalProxyDensity;%interp1(timeThisLat, finalProxyDensity, timestamps10sFixed, 'linear', 0);
    
    p.progress;
end
p.stop;

morningFtest = morningFtest';
eveningFtest = eveningFtest';

[morningTimestampsGrid, morningLatitudeGrid] = meshgrid(timestamps10sFixed, morningBins);
morningProxy = morningProxy';
[eveningTimestampsGrid, eveningLatitudeGrid] = meshgrid(timestamps10sFixed, eveningBins);
eveningProxy = eveningProxy';

% morningDensity = [];
% eveningDensity = [];
% morningTimes = [];
% eveningTimes = [];
% intervals = 3;
% intervalLength = floor(length(timestamps10sFixed) / intervals);
% for i = 1:intervals
%     if i < intervals
%         k = (i - 1) * intervalLength + 1 : i * intervalLength;
%     else
%         k = (i - 1) * intervalLength + 1 : length(timestamps10sFixed);
%     end
%     
%     beginTime = timestamps10sFixed(k(1));
%     endTime = timestamps10sFixed(k(end));
%     
%     morningInterp = interp2(morningTimestampsGrid(:,k), morningLatitudeGrid(:,k), morningProxy(:,k), morningTimestamps10s, morningLatitude, 'spline');
%     eveningInterp = interp2(eveningTimestampsGrid(:,k), eveningLatitudeGrid(:,k), eveningProxy(:,k), eveningTimestamps10s, eveningLatitude, 'spline');
%     
%     morningIndices = morningTimestamps10s >= beginTime & morningTimestamps10s <= endTime;
%     eveningIndices = eveningTimestamps10s >= beginTime & eveningTimestamps10s <= endTime;
%     morningInterp(~morningIndices) = [];
%     eveningInterp(~eveningIndices) = [];
%     morningTimes = [morningTimes; morningTimestamps10s(morningIndices)];
%     eveningTimes = [eveningTimes; eveningTimestamps10s(eveningIndices)];
%     
%     morningDensity = [morningDensity; morningInterp];
%     eveningDensity = [eveningDensity; eveningInterp];
% end
% 
% t = [morningTimes; eveningTimes];

morningDensity = interp2(morningTimestampsGrid, morningLatitudeGrid, morningProxy, morningTimestamps10s, morningLatitude, 'linear');
eveningDensity = interp2(eveningTimestampsGrid, eveningLatitudeGrid, eveningProxy, eveningTimestamps10s, eveningLatitude, 'linear');
t = [morningTimestamps10s; eveningTimestamps10s];
proxy = [morningDensity; eveningDensity];
[t, order, ~] = unique(t);
proxy = proxy(order);

tNoNans = t(~isnan(proxy));
proxyNoNans = proxy(~isnan(proxy));
tToInterpolate = t;

geomagneticPrediction = interp1(tNoNans, proxyNoNans, tToInterpolate, 'linear', 0);
predictedDensity = geomagneticPrediction + msis270kmNoAp;


end

function plotParametrizationResults(morningFourierGrid, eveningFourierGrid, morningBins, eveningBins)
%

[morningDoyGrid, morningBinGrid] = meshgrid(1:365, morningBins);
[eveningDoyGrid, eveningBinGrid] = meshgrid(1:365, eveningBins);

figure;
subplot(2,1,1);
surf(morningDoyGrid, morningBinGrid, morningFourierGrid);
title('Morning Fourier coeffs');
ylabel('Geographic Latitude')
xlabel('doy')
xlim([1 365]);
ylim([min(morningBinGrid(:,1)) max(morningBinGrid(:,1))]);
colorbar;
view(2);
shading interp

subplot(2,1,2);
surf(eveningDoyGrid, eveningBinGrid, eveningFourierGrid);
title('Evening Fourier coeffs');
ylabel('Geographic Latitude')
xlabel('doy')
xlim([1 365]);
ylim([min(eveningBinGrid(:,1)) max(eveningBinGrid(:,1))]);
colorbar;
view(2);
shading interp

end


function [timeByLatitude, latitudeBins] = computeTimeCells(allTimestamps10s, timestamps10s, latitude)
%

timestamps10s = unique(cell2vector(timestamps10s));
latitude = latitude(ismember(allTimestamps10s, timestamps10s));

[minLatitude, maxLatitude] = findInterpolationLimits(latitude);

latitudeBins = minLatitude + 1 : 3 : maxLatitude - 1;
if maxLatitude - 1 > max(latitudeBins)
    latitudeBins(end + 1) = max(latitudeBins) + 3;
end

timeByLatitude = cell(length(latitudeBins), 1);
indicesToRemove = false(length(latitudeBins),1);
parfor i = 1:length(latitudeBins)
    indices = (latitudeBins(i) - 1 < latitude & latitude <= latitudeBins(i) + 1);
    if isempty(find(indices, 1))
        indicesToRemove(i) = 1;
        continue;
    end
    timeByLatitude{i} = timestamps10s(indices);
end

latitudeBins = latitudeBins(~indicesToRemove);
timeByLatitude = timeByLatitude(~indicesToRemove);

end


function [vector] = cell2vector(cellArray)
%

vector = [];

for i = 1:numel(cellArray)
    thisCellValues = cellArray{i};
    
    if isrow(thisCellValues)
        thisCellValues = thisCellValues';
    end
    
    vector = [vector; thisCellValues];
end

end


function compareGoceDensityToModel(goceDensity, modelDensity, ae, timestampsAeDatenum, timestampsDensityDatenum, results)
%

ratio = goceDensity ./ modelDensity;

plotTimeseriesOfRatio(ratio, ae, timestampsAeDatenum, timestampsDensityDatenum);

plotCorrelation(modelDensity, goceDensity, 'NRLMSISE00 density', 'GOCE measured density', 1, results);

plotHistogramFits(ratio);

end


function plotTimeseriesOfRatio(ratio, ae, timestampsAeDatenum, timestampsDensityDatenum)
%

ratioTrend = removePeriodicBackground(ratio, 125, 6, 0);

figure;
[hAx,hLine1,hLine2] = plotyy(timestampsDensityDatenum, ratio, timestampsAeDatenum, ae);
set(hLine1, 'LineStyle', 'none', 'Marker', '.')
% For some reason Matlab can't handle axis limits properly, so they need to
% be set manually.
ratioYmin = 0.1 * floor(min(ratio) / 0.1);
ratioYmax = 0.1 * ceil(max(ratio) / 0.1);
set(hAx(1), 'YLim', [ratioYmin ratioYmax], 'YTick', ratioYmin:0.1:ratioYmax);
set(hAx, 'XLim', [min(timestampsDensityDatenum) max(timestampsDensityDatenum)], 'XMinorTick', 'on'); 
aeYmax = 500 * ceil(max(ae) / 500);
set(hAx(2),'YLim', [0 aeYmax], 'YTick', 0:250:aeYmax);
datetick(hAx(1), 'x', 'yyyy-mm', 'keepticks', 'keeplimits')
datetick(hAx(2), 'x', 'yyyy-mm', 'keepticks', 'keeplimits')
rotateticklabel(hAx(1), 90);
rotateticklabel(hAx(2), 90);
hold on;
plot(timestampsDensityDatenum, ratioTrend, 'r-');
hold off;
ylabel(hAx(1), 'GOCE / NRLMSISE00 density')
ylabel(hAx(2), 'AE')
title('Ratio of GOCE NON-normalized densities to NRLMSISE00 prediction')

meanRatio = mean(ratio);
errRatio = std(ratio);
meanReciprocal = mean(1 ./ ratio);
errReciprocal = std(1 ./ ratio);
textString1 = ['Mean GOCE / NRLMSISE00 density: ', num2str(meanRatio, '%.2f'), ' ± ', num2str(errRatio, '%.2f')];
textString2 = ['Mean NRLMSISE00 / GOCE density: ', num2str(meanReciprocal, '%.2f'), ' ± ', num2str(errReciprocal, '%.2f')];
     
ylimits = get(gca, 'ylim');
xlimits = get(gca, 'xlim');
textYLocation = ylimits(2) - 0.05 * (ylimits(2) - ylimits(1));
textXLocation = xlimits(2) - 0.05 * (xlimits(2) - xlimits(1));
text(textXLocation, textYLocation, textString1, 'FontSize', 9, ...
    'VerticalAlignment','top', 'HorizontalAlignment','right');

textYLocation = ylimits(2) - 0.1 * (ylimits(2) - ylimits(1));
textXLocation = xlimits(2) - 0.05 * (xlimits(2) - xlimits(1));
text(textXLocation, textYLocation, textString2, 'FontSize', 9, ...
    'VerticalAlignment','top', 'HorizontalAlignment','right');

end

function plotHistogramFits(ratio)
%

figure;
ratio = sort(ratio);
minX = floor(min(ratio) * 100) / 100;
maxX = ceil(max(ratio) * 100) / 100;

[numOfElementsInBin, x] = hist(ratio, minX:0.01:maxX);

pdNormal = fitdist(ratio, 'Normal');
pdLognormal = fitdist(ratio, 'Lognormal');
pdTlocationScale = fitdist(ratio, 'tlocationscale');

pdfNormal = pdf(pdNormal, ratio);
pdfLognormal = pdf(pdLognormal, ratio);
pdfTlocationScale = pdf(pdTlocationScale, ratio);

dx = diff(x(1:2));
bar(x,numOfElementsInBin / sum(numOfElementsInBin * dx))
title('Goce/Msis histogram with pdf fits')
xlabel('Goce/NRLMSISE00 density ratio')
ylabel('Percentage of data');
xlim([0.25 1.25])

hold all
normalHandle = plot(ratio, pdfNormal, 'r', 'Linewidth', 2);
lognHandle = plot(ratio, pdfLognormal, 'y', 'Linewidth', 2);
tlocHandle = plot(ratio, pdfTlocationScale, 'g', 'Linewidth', 2);
hold off

h = [normalHandle, lognHandle, tlocHandle];
legend(h, 'Normal', 'Lognormal', 't Location-Scale'); 

figure;
xValuesLognormal = random(pdLognormal, length(ratio), 1);
xValuesTlocScale = random(pdTlocationScale, length(ratio), 1);
xValuesTlocScale(xValuesTlocScale <= 0) = [];

tlocHandle = normplot(xValuesTlocScale);
hold all;
lognHandle = normplot(xValuesLognormal);
ratioHandle = normplot(ratio);

set(tlocHandle(1), 'Color', 'g')
set(lognHandle(1), 'Color', 'y')
set(tlocHandle(2), 'LineStyle', 'none')
set(lognHandle(2), 'LineStyle', 'none')
set(tlocHandle(3), 'LineStyle', 'none')
set(lognHandle(3), 'LineStyle', 'none')

h = [ratioHandle(3), lognHandle(1), tlocHandle(1), ratioHandle(1)];
legend(h, 'Normal', 'Lognormal', 't Location-Scale', 'Ratio', 'Location', 'NorthWest'); 

end

function th=rotateticklabel(h,rot,demo)
%ROTATETICKLABEL rotates tick labels
%   TH=ROTATETICKLABEL(H,ROT) is the calling form where H is a handle to
%   the axis that contains the XTickLabels that are to be rotated. ROT is
%   an optional parameter that specifies the angle of rotation. The default
%   angle is 90. TH is a handle to the text objects created. For long
%   strings such as those produced by datetick, you may have to adjust the
%   position of the axes so the labels don't get cut off.
%
%   Of course, GCA can be substituted for H if desired.
%
%   TH=ROTATETICKLABEL([],[],'demo') shows a demo figure.
%
%   Known deficiencies: if tick labels are raised to a power, the power
%   will be lost after rotation.
%
%   See also datetick.

% Copyright (c) 2005, Andrew Bliss
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

%   Written Oct 14, 2005 by Andy Bliss
%   Copyright 2005 by Andy Bliss

%DEMO:
if nargin==3
    x=[now-.7 now-.3 now];
    y=[20 35 15];
    figure
    plot(x,y,'.-')
    datetick('x',0,'keepticks')
    h=gca;
    set(h,'position',[0.13 0.35 0.775 0.55])
    rot=90;
end

%set the default rotation if user doesn't specify
if nargin==1
    rot=90;
end
%make sure the rotation is in the range 0:360 (brute force method)
while rot>360
    rot=rot-360;
end
while rot<0
    rot=rot+360;
end
%get current tick labels
a=get(h,'XTickLabel');
%erase current tick labels from figure
set(h,'XTickLabel',[]);
%get tick label positions
b=get(h,'XTick');
c=get(h,'YTick');
%make new tick labels
if rot<180
    th=text(b,repmat(c(1)-.1*(c(2)-c(1)),length(b),1),a,'HorizontalAlignment','right','rotation',rot);
else
    th=text(b,repmat(c(1)-.1*(c(2)-c(1)),length(b),1),a,'HorizontalAlignment','left','rotation',rot);
end

end


function [ae, ap, absB, vBz, akasofuEpsilon, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningAeProxy, eveningAeProxy, morningTimestamps10s, eveningTimestamps10s, timestamps1minOut, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = sliceToInterestingIntervals(ae, ap, absB, vBz, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningAeProxy, eveningAeProxy, morningTimestamps10s, eveningTimestamps10s, timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed,...
 density3h, morningMagneticLatitude, eveningMagneticLatitude, timestamps10sFixed, timestamps1min, timestampsAbsB, akasofuEpsilon, timestampsDatenum, intervalsOfInterest)
% 

aeTemp = ae; apTemp = ap;
averagedDensityNoBgTemp = averagedDensityNoBg; morningDensityNoBgTemp = morningDensityNoBg; 
eveningDensityNoBgTemp = eveningDensityNoBg; morningMsisDensityTemp = morningMsisDensity; eveningAeProxyTemp = eveningAeProxy;
eveningMsisDensityTemp = eveningMsisDensity; morningAeProxyTemp = morningAeProxy; morningTimestamps10sTemp = morningTimestamps10s; 
eveningTimestamps10sTemp = eveningTimestamps10s; timestamps1minFixedTemp = timestamps1minFixed; 
timestamps3hTemp = timestamps3h; timestamps3hFixedTemp = timestamps3hFixed; density3hTemp = density3h; morningMagneticLatitudeTemp = morningMagneticLatitude;
eveningMagneticLatitudeTemp = eveningMagneticLatitude; timestampsEpsilonTemp = timestampsEpsilon; akasofuEpsilonTemp = akasofuEpsilon;
vBzTemp = vBz; timestampsAbsBTemp = timestampsAbsB; absBTemp = absB; timestampsDatenumTemp = timestampsDatenum; 

ae = {}; ap = {};  averagedDensityNoBg = {}; morningDensityNoBg = {};
eveningDensityNoBg = {}; morningMsisDensity = {}; eveningMsisDensity = {}; morningTimestamps10s = {};
eveningTimestamps10s = {}; timestamps1minFixed = {}; timestamps3h = {}; density3h = {}; morningMagneticLatitude = {};
eveningMagneticLatitude = {}; timestampsAbsB = {}; akasofuEpsilon = {}; vBz = {}; timestampsEpsilon = {}; absB = {}; timestampsDatenum = {};
timestamps3hFixed = {}; morningAeProxy = {}; eveningAeProxy = {};

cellArrayLength = length(intervalsOfInterest(:,1));

for i = 1:cellArrayLength
    beginIndex = intervalsOfInterest(i,1);
    endIndex = intervalsOfInterest(i,2);
    timestamps1minOut{i} = timestamps1min(beginIndex:endIndex);
    ae{i} = aeTemp(beginIndex:endIndex);
    
    threeHinSec = 3 * 60 * 60;
    [~, beginIndex3h] = min(abs(timestamps3hTemp - threeHinSec * round(timestamps1min(beginIndex) / threeHinSec)));
    [~, endIndex3h] = min(abs(timestamps3hTemp - threeHinSec * round(timestamps1min(endIndex) / threeHinSec)));
    ap{i} = apTemp(beginIndex3h:endIndex3h);
    timestamps3h{i} = timestamps3hTemp(beginIndex3h:endIndex3h);
    
    [~, beginIndex3hFixed] = min(abs(timestamps3hFixedTemp - threeHinSec * round(timestamps1min(beginIndex) / threeHinSec)));
    [~, endIndex3hFixed] = min(abs(timestamps3hFixedTemp - threeHinSec * round(timestamps1min(endIndex) / threeHinSec)));
    density3h{i} = density3hTemp(beginIndex3hFixed:endIndex3hFixed);
    timestamps3hFixed{i} = timestamps3hFixedTemp(beginIndex3hFixed:endIndex3hFixed);
    
    [~, beginIndexAbsB] = min(abs(timestampsAbsBTemp - timestamps1min(beginIndex)));
    [~, endIndexAbsB] = min(abs(timestampsAbsBTemp - timestamps1min(endIndex)));
    absB{i} = absBTemp(beginIndexAbsB:endIndexAbsB);
    timestampsAbsB{i} = timestampsAbsBTemp(beginIndexAbsB:endIndexAbsB);
    
    [~, beginIndexEpsilon] = min(abs(timestampsEpsilonTemp - timestamps1min(beginIndex)));
    [~, endIndexEpsilon] = min(abs(timestampsEpsilonTemp - timestamps1min(endIndex)));
    akasofuEpsilon{i} = akasofuEpsilonTemp(beginIndexEpsilon:endIndexEpsilon);
    vBz{i} = vBzTemp(beginIndexEpsilon:endIndexEpsilon);
    timestampsEpsilon{i} = timestampsEpsilonTemp(beginIndexEpsilon:endIndexEpsilon);
    
    [~, beginIndex1min] = min(abs(timestamps1minFixedTemp - timestamps1min(beginIndex)));
    [~, endIndex1min] = min(abs(timestamps1minFixedTemp - timestamps1min(endIndex)));
    averagedDensityNoBg{i} = averagedDensityNoBgTemp(beginIndex1min:endIndex1min);
    timestamps1minFixed{i} = timestamps1minFixedTemp(beginIndex1min:endIndex1min);
    
    [~, beginIndexMorning10s] = min(abs(morningTimestamps10sTemp - timestamps1min(beginIndex)));
    [~, endIndexMorning10s] = min(abs(morningTimestamps10sTemp - timestamps1min(endIndex)));   
    [~, beginIndexEvening10s] = min(abs(eveningTimestamps10sTemp - timestamps1min(beginIndex)));
    [~, endIndexEvening10s] = min(abs(eveningTimestamps10sTemp - timestamps1min(endIndex)));
    morningMsisDensity{i} = morningMsisDensityTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningMsisDensity{i} = eveningMsisDensityTemp(beginIndexEvening10s:endIndexEvening10s);
    morningDensityNoBg{i} = morningDensityNoBgTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningDensityNoBg{i} = eveningDensityNoBgTemp(beginIndexEvening10s:endIndexEvening10s);
    morningTimestamps10s{i} = morningTimestamps10sTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningTimestamps10s{i} = eveningTimestamps10sTemp(beginIndexEvening10s:endIndexEvening10s);
    morningMagneticLatitude{i} = morningMagneticLatitudeTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningMagneticLatitude{i} = eveningMagneticLatitudeTemp(beginIndexEvening10s:endIndexEvening10s);
    morningAeProxy{i} = morningAeProxyTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningAeProxy{i} = eveningAeProxyTemp(beginIndexEvening10s:endIndexEvening10s);
    
    [~, beginIndex10s] = min(abs(timestamps10sFixed - timestamps1min(beginIndex)));
    [~, endIndex10s] = min(abs(timestamps10sFixed - timestamps1min(endIndex)));
    timestampsDatenum{i} = timestampsDatenumTemp(beginIndex10s:endIndex10s);
end

end
