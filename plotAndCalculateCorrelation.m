function [results, averageIntegral, timestampsAverInt] = plotAndCalculateCorrelation(firstDatenum, timestampsGeom, timestampsMorning, timestampsEvening, geomIndex, morningDensity, eveningDensity, latitude, latitudeTimestamps, indexName, plotFigures, results, timeseriesHandle)
% [r, r2] = plotAndCalculateCorrelation(ae, averagedDensity, timelag)

density = [morningDensity; eveningDensity];
timestampsFixed = [timestampsMorning; timestampsEvening];
[timestampsFixed, order] = unique(timestampsFixed);
density = density(order);

if round(mean(diff(timestampsGeom)) / 3600) <= 2
    density = smooth(density, 541);
    latitude = latitude(ismember(latitudeTimestamps, timestampsFixed));
    timestampsNorth = timestampsFixed(latitude > 80);
    times = timestampsNorth(find(diff(timestampsNorth) > 45 * 60) + 1);
    ind = ismember(timestampsFixed, times);
    density = density(ind);
    timestampsFixed = timestampsFixed(ind);
    timestampsFixed = 60 * round(timestampsFixed / 60);
end

[timelag, timelagInHours, bestIntegral, averageGoodLag, averageIntegral] = compareDensityToGeomIndexIntegral(density, geomIndex, timestampsGeom, timestampsFixed, indexName, plotFigures);

[rowNum, ~] = size(results);
emptyCells = cellfun(@isempty,results);
[~, emptyColPositions] = find(emptyCells);
colNum = min(emptyColPositions);
if rowNum == 2
    colNum = length(results(rowNum,:)) + 1;
    if strcmpi(indexName, 'ae')
        results{1, colNum} = 'Max Aver. AE Int.';
        results{1, colNum + 1} = ['Int window: ',indexName];
    else
        results{1, colNum} = ['Int window: ',indexName];
    end
end

if strcmpi(indexName, 'ae')
    results{rowNum, colNum} = max(averageIntegral);
    results{rowNum, colNum + 1} = timelagInHours;
else
    results{rowNum, colNum} = timelagInHours;
end

timestampsFixed = timestampsFixed(ismember(timestampsFixed, timestampsGeom));
densityAverIntIndices = ismember(timestampsFixed, timestampsGeom(averageGoodLag + 1:end));

if strcmpi(indexName, 'ae') && plotFigures ~=0
    figure(timeseriesHandle);
    subplot(2,2,2)
    secondsInDay = 24 * 60 * 60;
    timestampsInDays = timestampsGeom(averageGoodLag + 1:end) / secondsInDay + firstDatenum;
    timestampsInDaysFixed = timestampsFixed(densityAverIntIndices) / secondsInDay + firstDatenum;
    [hAx,~,~] = plotyy(timestampsInDays, averageIntegral, timestampsInDaysFixed, density(densityAverIntIndices));
    ylabel(hAx(1), 'AE Integral')
    ylabel(hAx(2), 'Density')
    title('Raw AE integral vs FFT smoothed density');
    set(hAx, 'XLim', [min(timestampsInDaysFixed) max(timestampsInDaysFixed)]);
    datetick(hAx(1), 'x', 'dd', 'keeplimits')
    datetick(hAx(2), 'x', 'dd', 'keeplimits')
    grid on;
    
    figure('Color', 'white');
    subplot(3,1,1)
    plot(timestampsInDaysFixed, density(densityAverIntIndices) * 1e-12, 'k')
    ylabel('Density (kg/m^3)', 'FontSize', 18);
    axis tight
    datetick('x', 'dd', 'keeplimits')
    
    subplot(3,1,2)
    plot(timestampsInDays, averageIntegral, 'k')
    ylabel('21-hour AE Integral', 'FontSize', 18);
    axis tight
    datetick('x', 'dd', 'keeplimits')
    
    subplot(3,1,3)
    plot(timestampsInDays, geomIndex(averageGoodLag + 1:end), 'k')
    xlabel('Calendar day on May 2011', 'FontSize', 18);
    ylabel('AE (nT)', 'FontSize', 14);
    axis tight
    datetick('x', 'dd', 'keeplimits')
end

densityBestIntIndices = ismember(timestampsFixed, timestampsGeom(timelag + 1:end));
bestIntegral = bestIntegral(ismember(timestampsGeom(timelag + 1:end), timestampsFixed));
results = plotCorrelation(bestIntegral, density(densityBestIntIndices), [indexName, ' Best Integral'], 'Density', plotFigures, results);

averageIntIndices = ismember(timestampsFixed, timestampsGeom(averageGoodLag + 1:end));
densityAverageInt = density(averageIntIndices);
timestampsAverInt = timestampsGeom(averageGoodLag + 1:end);
averageIntegralFixed = averageIntegral(ismember(timestampsAverInt, timestampsFixed));
results = plotCorrelation(averageIntegralFixed, densityAverageInt, [indexName, ' Average Integral'], 'Density', plotFigures, results);

if round(mean(diff(timestampsGeom)) / 3600) <= 2
    geomIndex = smooth(geomIndex, 90);
end
shortenedIndicesContinuous = ismember(timestampsGeom, timestampsAverInt);
geomIndex = geomIndex(shortenedIndicesContinuous);
timestampsGeom = timestampsGeom(shortenedIndicesContinuous);
shortenedIndicesGaps = ismember(timestampsFixed, timestampsGeom);
timestampsFixed = timestampsFixed(shortenedIndicesGaps);
density = density(shortenedIndicesGaps);

geomIndexFixed = geomIndex(ismember(timestampsGeom, timestampsFixed));
densityFixed = density(ismember(timestampsFixed, timestampsGeom));

if strcmpi(indexName, 'Akasofu Epsilon')
    % Solar indices are lagged by 6h
    % Cut 6h from the end of timestamps
    lagtime = 540; %min
elseif ~isempty(strfind(upper(indexName), '|B|'))
    lagtime = 630; %min
elseif ~isempty(strfind(upper(indexName), '|V|'))
    lagtime = 540; %min
elseif ~isempty(strfind(upper(indexName), 'AE'))
    lagtime = 450; %min
else
    lagtime = 360;
end

timestampsLimitedFromEnd = timestampsGeom(ismember(timestampsGeom, timestampsGeom - lagtime * 60));
% Lag these limited timestamps by 6h
timestampsLagged = timestampsLimitedFromEnd + lagtime * 60;
% Pick geomIndex as geomIndex([tFixed(1), tfixed(2), ...])
geomIndex6hAgo = geomIndex(ismember(timestampsGeom, timestampsLimitedFromEnd));
geomIndex6hAgo = geomIndex6hAgo(ismember(timestampsLagged, timestampsFixed));
% Pick density values as density([tFixed(1)+6h, tFixed(2)+6h, ...])
densityShorter = density(ismember(timestampsFixed, timestampsLagged));
results = plotCorrelation(geomIndex6hAgo, densityShorter, indexName, 'Density at 270 km', plotFigures, results);

end

function [bestLag, bestLagInHours, geomIndexBestInt, averageGoodLag, averageLagInt] = compareDensityToGeomIndexIntegral(density, ...
    geomIndex, timestamps, timestampsFixed, indexName, plotFigures)
%
maxDays = 3;
if strcmpi(indexName, 'ap'); maxLag = maxDays * 8;
else maxLag = maxDays * 24 * 60; end

indicesInHour = 60;
if strcmpi(indexName, 'ae')
    averageGoodLag = 21 * indicesInHour;
elseif strcmpi(indexName, 'Akasofu Epsilon')
    averageGoodLag = 30 * indicesInHour;
elseif ~isempty(strfind(upper(indexName), '|B|'))
    averageGoodLag = 41 * indicesInHour;
elseif ~isempty(strfind(upper(indexName), '|V|'))
    averageGoodLag = 43 * indicesInHour;
else
    averageGoodLag = 8;
end

cumulativeGeomIndex = cumsum(geomIndex);
correlations = zeros(1, maxLag + 1);
lags = zeros(1, maxLag + 1);
correlationType = 'Spearman';
parfor lag = 1:maxLag
    geomIndexInt = cumulativeGeomIndex(lag + 1 : end) - cumulativeGeomIndex(1 : end - lag);
    timestampsInt = timestamps(lag + 1 : end);
    densityIndices = ismember(timestampsFixed, timestampsInt);
    integralIndices = ismember(timestampsInt, timestampsFixed);
    correlations(lag + 1) = corr(density(densityIndices), geomIndexInt(integralIndices), 'type', correlationType);
    lags(lag + 1) = lag;
end
correlations(1) = corr(density(ismember(timestampsFixed, timestamps)), geomIndex(ismember(timestamps, timestampsFixed)),...
    'type', correlationType);
lags(1) = 0;
if strcmpi(indexName, 'ap'); lagsInHours = lags * 3;
else lagsInHours = lags / 60; end

maxCorr = abs(correlations) == max(abs(correlations));
bestLag = lags(find(maxCorr, 1, 'last'));
bestLagInHours = lagsInHours(find(maxCorr, 1, 'last'));
if bestLag > 0
    geomIndexBestInt = cumulativeGeomIndex(bestLag + 1 : end) - cumulativeGeomIndex(1 : end - bestLag);
else
    geomIndexBestInt = geomIndex;
end
averageLagInt = cumulativeGeomIndex(averageGoodLag + 1 : end) - cumulativeGeomIndex(1 : end - averageGoodLag);
fprintf('\n');

if plotFigures ~= 0
    figure;
    plot(lagsInHours, correlations);
    title([indexName, ' integral optimal window length'])
    ylabel([correlationType, ' correlation']);
    xlabel('lags / h')
    correlationMaximum = max(abs(correlations));
    correlationMaximum = correlationMaximum(1);
    integralWindowSize = lagsInHours(abs(correlations) == correlationMaximum);
    integralWindowSize = integralWindowSize(1);
    ylimits = get(gca, 'ylim');
    line([integralWindowSize integralWindowSize], [ylimits(1) correlationMaximum], 'LineStyle', '--');
    textYLocation = mean([ylimits(1) correlationMaximum]);
    text(integralWindowSize, textYLocation, ['\leftarrow Lag = ', num2str(integralWindowSize), ' h'], 'FontSize', 14)
end
% if strcmpi(indexName, 'ae'); integralWindowSize = integralWindowSize * 60;
% elseif ~isempty(strfind(upper(indexName), '|B|')); integralWindowSize = integralWindowSize * 60 / 4;
% else integralWindowSize = integralWindowSize / 3;end

end