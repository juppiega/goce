function plotPeakAnalysis(limitedTimestamps, ae, densityByLatitude, crossingTimes, minAllowedLatitude, maxAllowedLatitude, timeOfDay)
% writeAndPlotPeakAnalysis(timestamps, ae, splinedDensity, computeLatitutes)

persistent aeDensityLagByLatHandle
%global densityPeakMaxTimesByLatHandle

ae = removePeriodicBackground(ae, 432, 1, 0); % 432 min =^ 0.3 days
% ae = smooth(ae, 431);
[ae, ~, ~] = limitToNearPeak(ae, 'noSmooth', 'mean');
for i = 1:length(densityByLatitude(1,:))
   samplingFreq = 86400 / mean(diff(crossingTimes(:,i)));
   densityByLatitude(:,i) = removePeriodicBackground(densityByLatitude(:,i), 0.3, samplingFreq, 0);
end

oneDegreeStep = minAllowedLatitude:maxAllowedLatitude;
parfor i = 1:length(oneDegreeStep)
    splinedDensity(:,i) = spline(crossingTimes(:,i), densityByLatitude(:,i), limitedTimestamps);
end

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
    
    maxlag = 24 * 60;
    lagsAllIntervalLats = zeros(length(densityInterval), 1);
    maxIndicesAllIntervalLats = zeros(length(densityInterval), 1);
    for k = densityInterval
       [correlations, lags] = xcorr(splinedDensity(:,k), ae, maxlag, 'coeff');
       correlations = correlations(lags >= 0);
       lags = lags(lags >= 0); 
       allLatitudesAeLag(k) = lags(correlations == max(correlations));
       lagsAllIntervalLats(k - min(densityInterval) + 1) = lags(correlations == max(correlations));
    end
    
    if catTimelag > 0
        timelag = [timelag mean(lagsAllIntervalLats)];
        errInLag = [errInLag (std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats)))];  
    else
        timelag(length(i:-step:minAllowedLatitude) + 1) = mean(lagsAllIntervalLats);
        errInLag(length(i:-step:minAllowedLatitude) + 1) = (std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats)));
    end
end

southernInterval = -step/2:-step:minAllowedLatitude;
southernInterval(ismember(southernInterval, minAllowedLatitude)) = [];
for i = southernInterval
    densityInterval = i - minAllowedLatitude + 1 - step : i - minAllowedLatitude + 1;
    if ~isempty(find(densityInterval < 1, 1))
        densityInterval = 1 : i - minAllowedLatitude + 1;
    end
    
    maxlag = 24 * 60;
    lagsAllIntervalLats = zeros(length(densityInterval), 1);
    maxIndicesAllIntervalLats = zeros(length(densityInterval), 1);
    for k = densityInterval
       [correlations, lags] = xcorr(splinedDensity(:,k), ae, maxlag, 'coeff');
       correlations = correlations(lags >= 0);
       lags = lags(lags >= 0); 
       allLatitudesAeLag(k) = lags(correlations == max(correlations));
       lagsAllIntervalLats(k - min(densityInterval) + 1) = lags(correlations == max(correlations));    
    end
    
    timelag(length(i:-step:minAllowedLatitude)) = mean(lagsAllIntervalLats);
    errInLag(length(i:-step:minAllowedLatitude)) = std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats));
end

intervalLatitudes = (southernInterval(end) : step : northernInterval(end)) + step/2;
intervalLatitudes(end) = [];
intervalLatitudes = [ round(mean([southernInterval(end) minAllowedLatitude]))...
                      intervalLatitudes                                      ...
                      round(mean([northernInterval(end) maxAllowedLatitude])) ]; 
                  
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