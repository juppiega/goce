function [results, averageIntegral, timestampsAverInt] = plotAndCalculateCorrelation(firstDatenum, timestamps, timestampsFixed, geomIndex, density, indexName, plotFigures, results, timeseriesHandle)
% [r, r2] = plotAndCalculateCorrelation(ae, averagedDensity, timelag)


[timelag, timelagInHours, bestIntegral, averageGoodLag, averageIntegral] = compareDensityToGeomIndexIntegral(density, geomIndex, timestamps, timestampsFixed, indexName, plotFigures);

%figure;
%plotyy(timestamps, geomIndex, timestampsFixed, density);

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

timestampsFixed = timestampsFixed(ismember(timestampsFixed, timestamps));
densityBestIntIndices = ismember(timestampsFixed, timestamps(timelag + 1:end));
    
if strcmpi(indexName, 'ae') && plotFigures ~=0
    figure(timeseriesHandle);
    subplot(2,2,2)
    secondsInDay = 24 * 60 * 60;
    timestampsInDays = timestamps(timelag + 1:end) / secondsInDay + firstDatenum;
    timestampsInDaysFixed = timestampsFixed(densityBestIntIndices) / secondsInDay + firstDatenum;
    [hAx,~,~] = plotyy(timestampsInDays, bestIntegral, timestampsInDaysFixed, density(densityBestIntIndices));
    ylabel(hAx(1), 'AE Integral')
    ylabel(hAx(2), 'Density')
    title('Raw AE integral vs FFT smoothed density');
    set(hAx, 'XLim', [min(timestampsInDaysFixed) max(timestampsInDaysFixed)]);
    datetick(hAx(1), 'x', 'dd', 'keeplimits')
    datetick(hAx(2), 'x', 'dd', 'keeplimits')
    grid on;
end

bestIntegral = bestIntegral(ismember(timestamps(timelag + 1:end), timestampsFixed));
results = plotCorrelation(bestIntegral, density(densityBestIntIndices), [indexName, ' Best Integral'], 'Density', plotFigures, results);

averageIntIndices = ismember(timestampsFixed, timestamps(averageGoodLag + 1:end));
densityAverageInt = density(averageIntIndices);
timestampsAverInt = timestamps(averageGoodLag + 1:end);
averageIntegralFixed = averageIntegral(ismember(timestampsAverInt, timestampsFixed));
results = plotCorrelation(averageIntegralFixed, densityAverageInt, [indexName, ' Average Integral'], 'Density', plotFigures, results);

shortenedIndicesContinuous = ismember(timestamps, timestampsAverInt);
geomIndex = geomIndex(shortenedIndicesContinuous);
timestamps = timestamps(shortenedIndicesContinuous);
shortenedIndicesGaps = ismember(timestampsFixed, timestamps);
timestampsFixed = timestampsFixed(shortenedIndicesGaps);
density = density(shortenedIndicesGaps);

geomIndexFixed = geomIndex(ismember(timestamps, timestampsFixed));
densityFixed = density(ismember(timestampsFixed, timestamps));
if strcmpi(indexName, 'ae')
    geomIndexNoBg = removePeriodicBackground(geomIndexFixed, 125, 1, 0);
    geomIndexNoBg = normalize(geomIndexNoBg, geomIndexFixed);
    results = plotCorrelation(geomIndexNoBg, densityFixed, indexName, 'Density at 270 km', plotFigures, results);
elseif strcmpi(indexName, 'Akasofu Epsilon') || ~isempty(strfind(upper(indexName), '|B|')) ...
        || ~isempty(strfind(upper(indexName), '|V|'))
    % Solar indices are lagged by 6h
    % Cut 6h from the end of timestamps
    timestampsLimitedFromEnd = timestamps(ismember(timestamps, timestamps - 6 * 60 * 60));
    % Lag these limited timestamps by 6h
    timestampsLagged = timestampsLimitedFromEnd + 6 * 60 * 60;
    % Pick geomIndex as geomIndex([tFixed(1), tfixed(2), ...])
    geomIndex6hAgo = geomIndex(ismember(timestamps, timestampsLimitedFromEnd));
    geomIndex6hAgo = geomIndex6hAgo(ismember(timestampsLagged, timestampsFixed));
    % Pick density values as density([tFixed(1)+6h, tFixed(2)+6h, ...])
    densityShorter = density(ismember(timestampsFixed, timestampsLagged));
    results = plotCorrelation(geomIndex6hAgo, densityShorter, indexName, 'Density at 270 km', plotFigures, results);
else % ap
    results = plotCorrelation(geomIndexFixed, densityFixed, indexName, 'Density at 270 km', plotFigures, results);
end

end

function [bestLag, bestLagInHours, geomIndexBestInt, averageGoodLag, averageLagInt] = compareDensityToGeomIndexIntegral(density, ...
    geomIndex, timestamps, timestampsFixed, indexName, plotFigures)
%
maxDays = 3;
if strcmpi(indexName, 'ap'); maxLag = maxDays * 8;
else maxLag = maxDays * 24 * 60; end

indicesInHour = 60;
if strcmpi(indexName, 'ae')
    averageGoodLag = 22 * indicesInHour;
elseif strcmpi(indexName, 'Akasofu Epsilon')
    averageGoodLag = 38 * indicesInHour;
elseif ~isempty(strfind(upper(indexName), '|B|'))
    averageGoodLag = 34 * indicesInHour;
elseif ~isempty(strfind(upper(indexName), '|V|'))
    averageGoodLag = 34 * indicesInHour;
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

bestLag = lags(find(correlations == max(correlations), 1, 'last'));
bestLagInHours = lagsInHours(find(correlations == max(correlations), 1, 'last'));
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
    integralWindowSize = lagsInHours(correlations == max(correlations));
    ylimits = get(gca, 'ylim');
    line([integralWindowSize integralWindowSize], [ylimits(1) max(correlations)], 'LineStyle', '--');
    textYLocation = mean([ylimits(1) max(correlations)]);
    text(integralWindowSize, textYLocation, ['\leftarrow Lag = ', num2str(integralWindowSize), ' h'], 'FontSize', 14)
end
% if strcmpi(indexName, 'ae'); integralWindowSize = integralWindowSize * 60;
% elseif ~isempty(strfind(upper(indexName), '|B|')); integralWindowSize = integralWindowSize * 60 / 4;
% else integralWindowSize = integralWindowSize / 3;end

end