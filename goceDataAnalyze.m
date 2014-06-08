function goceDataAnalyze( varargin )
% goceDensityAExcorr( AEFilename, densityFilename [,beginDay, endDay])
%   Detailed explanation goes here

numOfInputArgs = nargin;
inputArgs = varargin;
[ AEFilename, densityFilename, beginDay, endDay ] = processInputArguments(inputArgs, numOfInputArgs);

[ae, averagedDensity, timestamps1min, averagedLatitude, averagedAltitude, averagedSolarTime, ...
    density, latitude, altitude, solarTime, timestamps10s, timestampsDatenum, longtitude] ...
    = readAEDensityTimestampsLatitudeAltitudeAndSolarTime( AEFilename, densityFilename, beginDay, endDay );

plotTimeseries(timestamps1min, ae, averagedDensity, 1)

[averagedDensityNoBg] = removePeriodicDensityBackground(averagedDensity);

plotTimeseries(timestamps1min, ae, averagedDensityNoBg, 4)

%figure(7);


[correlations, lags] = xcorr(averagedDensityNoBg, ae, 1440, 'coeff'); 
[timelag] = plotAndGiveMaxCrossCorrelation(correlations, lags);

[movedAe, movedAveragedDensity] = moveDataseries(ae, averagedDensityNoBg, timelag); % Viepa tama tuonne seuraavan funktion sisaan!
[r, r2] = plotAndCalculateCorrelation(movedAe, movedAveragedDensity);
sprintf('%s %f', 'Pearson correlation: ', r)
%sprintf('%s %d %s', 'Thus, ', round(r2*100), '% of variation in Density can be explained by changes in AE')

%plot3(morningTimestamps, morningLatitude, morningDensity, 10)
% TEE TAMA LOPUKSI TARKOILLA ARVOILLA!!!
%[correctedDensity] = subtractSimpleModelFromDensity(morningAe, morningDensity, morningAltitude, morningSolarTime);
%[correctedDensity] = subtractSimpleModelFromDensity(ae, density, altitude, solarTime);

mi = find(solarTime < 12);
mdensity = density(mi);
malt = altitude(mi);
mtime = timestampsDatenum(mi);
mlat = latitude(mi);
mlon = longtitude(mi);

[correctedDensity] = subtractMsisFromDensity(mdensity, malt, mtime,  mlat, mlon, 100, 100, 10, 10000);

%plotDensityLatitudeTimeSurf(density, latitude, timestamps10s)
%plotDensityLatitudeTimeSurf(correctedDensity, latitude, timestamps10s)
%plotDensityLatitudeTimeSurf(morningDensity, morningLatitude, morningTimestamps)
%plotDensityLatitudeTimeSurf(correctedDensity, mlat, mtime)

%simpleDataPlot(correctedDensity, latitude, timestamps10s);
simpleDataPlot(correctedDensity, mlat, mtime);

end

function [ AEFilename, densityFilename, beginDay, endDay ] = processInputArguments(inputArgs, numOfInputArgs)
% [ AEFilename, densityFilename, beginDay, endDay ] = processInputArguments(varargin, nargin)

if numOfInputArgs == 2
    AEFilename = inputArgs{1};
    densityFilename = inputArgs{2};
    beginDay = 1;
    endDay = -1;
    
elseif numOfInputArgs == 4
    AEFilename = inputArgs{1};
    densityFilename = inputArgs{2};
    beginDay = inputArgs{3};
    endDay = inputArgs{4};
    
else
    sprintf('%s %d %s', 'Wrong number of input arguments: ', numOfInputArgs, '. See >>help goceDensityAExcorr')
    throw(MException('VerifyInput:InputArgumentError'))
end

end

function [ae, averagedDensity, timestamps1min, averagedLatitude, averagedAltitude, averagedSolarTime, ...
    density, latitude, altitude, solarTime, timestamps10s, timestampsDatenum, longtitude] ...
    = readAEDensityTimestampsLatitudeAltitudeAndSolarTime( AEFilename, densityFilename, beginDay, endDay )
% [ae, averagedDensity, timestamps] = readAEDensityAndTimestamps( AEFilename, densityFilename )

aeFile = fopen(AEFilename);
densityFile = fopen(densityFilename);

if aeFile == -1
    disp('aeFile open unsuccesful')
    return;
end

if densityFile == -1
    disp('densityFile open unsuccesful')
    return;
end

numOfHeaderLines = 15;
numOfAEValsInDay = 24 * 60;
firstAEline = numOfHeaderLines + (beginDay - 1) * numOfAEValsInDay;
if endDay > 0
    blockSizeAE = (endDay - beginDay + 1) * numOfAEValsInDay;
    aeData = textscan(aeFile, '%s %s %f %f %f %f %f', blockSizeAE, 'MultipleDelimsAsOne',1, 'HeaderLines',firstAEline);
else
    aeData = textscan(aeFile, '%s %s %f %f %f %f %f', 'MultipleDelimsAsOne',1, 'HeaderLines',firstAEline);
end

numOfDensityValsInDay = 24 * 60 * 6;
firstDensityLine = (beginDay - 1) * numOfDensityValsInDay;
if endDay > 0
    blockSizeDensity = (endDay - beginDay + 1) * numOfDensityValsInDay;
    densityData = textscan(densityFile, '%s %s %s %f %f %f %f %f %f %f %f %f', blockSizeDensity, ...
        'commentStyle','#', 'HeaderLines',firstDensityLine);
else
    densityData = textscan(densityFile, '%s %s %s %f %f %f %f %f %f %f %f %f', ...
        'commentStyle','#', 'HeaderLines',firstDensityLine);
end

ae = aeData{4};
density = densityData{9};
longtitude = densityData{5};
latitude = densityData{6};
altitude = densityData{4};
solarTime = densityData{7};

averagedDensity = smooth(density, 7);
averagedDensity = averagedDensity(1:6:end);
ae = ae(1:length(averagedDensity));

averagedLatitude = smooth(latitude, 7);
averagedLatitude = averagedLatitude(1:6:end);

averagedAltitude = smooth(altitude, 7);
averagedAltitude = averagedAltitude(1:6:end);

averagedSolarTime = smooth(solarTime, 7);
averagedSolarTime = averagedSolarTime(1:6:end);

timestampString = strcat(densityData{1}, densityData{2});
timestampsDatenum = datenum(timestampString, 'yyyy-mm-ddHH:MM:SS');
timestamps10s = timestampsDatenum - timestampsDatenum(1);
timestamps1min = timestamps10s(1:length(averagedDensity));

if fclose('all') ~= 0
    disp('File close unsuccesful')
end

end

function [averagedDensity] = removePeriodicDensityBackground( averagedDensity )
% [averagedDensity] = removePeriodicDensityBackground( averagedDensity )

numOfSamples = length(averagedDensity);
densityFFT = fft(averagedDensity);
plotFreq = (0:numOfSamples - 1) / numOfSamples; % Unit is 1 / min
period = 1 ./ plotFreq;
plotIndices = find(period < 1500);

h(1) = subplot(2,2,2);
plot(period(plotIndices), abs(densityFFT(plotIndices)))
title('Density FFT')
xlabel('T / min')

indicesToRemove = period < 100;
densityFFT(indicesToRemove) = 0;

h(2) = subplot(2,2,3);
plot(period(plotIndices), abs(densityFFT(plotIndices)))
linkaxes(h)
title('Stripped density FFT')
xlabel('T / min')

averagedDensity = abs(ifft(densityFFT));

end

function plotTimeseries(timestamps, ae, averagedDensity, plotNumber)
% plotTimeseries(timestamps, ae, averagedDensity)

subplot(2,2,plotNumber)
secondsInDay = 60 * 60 * 24;
timestampsInDays = timestamps / secondsInDay;
[hAx,hLine1,hLine2] = plotyy(timestampsInDays, ae, timestampsInDays, averagedDensity);
title('Timeseries of AE and density')
xlabel('t / days')
ylabel(hAx(1), 'AE')
ylabel(hAx(2), 'Density')
grid on

end

function [densityAETimelag] = plotAndGiveMaxCrossCorrelation(correlations, lags)
% plotCrossCorrelation(correlations, lags)

lagsInHours = lags/60;
figure(2);
plot(lagsInHours, correlations);
title('AE and Density xcorr')
xlabel('lags / h')
densityAETimelag = lagsInHours(correlations == max(correlations));
display(densityAETimelag)
maxCorrel = max(correlations);
ylimits = get(gca, 'ylim');
line([densityAETimelag densityAETimelag], [ylimits(1) maxCorrel], 'LineStyle', '--');

densityAETimelag = densityAETimelag * 60;

end

function [ae, averagedDensity] = moveDataseries(ae, averagedDensity, timelag)
% [ae, averagedDensity] = moveDataseries(ae, averagedDensity, timelag)

if timelag > 0
   averagedDensity = averagedDensity((timelag + 1) : end);
   ae = ae(1 : (length(ae) - timelag));
elseif timelag < 0
   averagedDensity = averagedDensity(1 : (length(averagedDensity) - timelag));
   ae = ae((timelag + 1) : end);
end

end

function [r, r2] = plotAndCalculateCorrelation(ae, averagedDensity)
% [r, r2] = plotAndCalculateCorrelation(ae, averagedDensity)

figure(3);
plot(ae, averagedDensity, '.')
p = polyfit(ae, averagedDensity, 1);
m = p(1);
b = p(2);
ylimits = get(gca, 'ylim');
xlimits = get(gca, 'xlim');
regLineYAxisCross = m * xlimits(1) + b;
regLineXAxisCross = (ylimits(2) - b) / m;
line([xlimits(1) regLineXAxisCross], [regLineYAxisCross ylimits(2)]);
title('Density vs AE (at most likely timelag)')
xlabel('AE')
ylabel('Density')

r = corr(ae, averagedDensity);
r2 = r * r;
end

function plotDensityLatitudeTimeSurf(correctedDensity, latitude, timestamps)
% plotDensityLatitudeTimeSurf(averagedDensity, averagedLatitude, timestamps)

figure(4);
secondsInDay = 60 * 60 * 24;
timestampsInDays = timestamps / secondsInDay;
[tInterp, latitudeInterp] = meshgrid(0:120/secondsInDay:timestampsInDays(length(timestampsInDays)), ...
    min(latitude):1:max(latitude));
densityInterp = griddata(timestampsInDays, latitude, correctedDensity, tInterp, latitudeInterp, 'cubic');
surf(tInterp, latitudeInterp, densityInterp, 'EdgeColor', 'None')
caxis([min(correctedDensity) max(correctedDensity)])
colormap jet(500)
colorbar
view(2);
xlabel('t / days')
ylabel('Lat (°)')
zlabel('Density')

% hold on;
% maxZ = ones(size(timestamps)) * max(correctedDensity);
% plot3(timestampsInDays, averagedLatitude, maxZ, 'r--');
% hold off;
end

function [correctedDensity] = subtractSimpleModelFromDensity(ae, density, altitude, solarTime)
% [correctedDensity] = subtractSimpleModelFromDensity(ae, averagedDensity, averagedAltitude)

[calmIndices] = findCalmPeriod(ae);
modelDensity = density(calmIndices);
modelAltitude = altitude(calmIndices);
modelSolarTime = solarTime(calmIndices);

F = scatteredInterpolant(modelAltitude, modelSolarTime, modelDensity, 'linear', 'nearest'); % Convex hull interpolation
densityModel = zeros(length(density), 1);

%pointsOutsideFAltitude = find((averagedAltitude > max(modelAltitude)) | (averagedAltitude < min(modelAltitude)));
%pointsOutsideFSolarTime = find(averagedSolarTime > max(modelSolarTime) | averagedSolarTime < min(modelSolarTime));
%pointsOutsideF = vertcat(pointsOutsideFAltitude, pointsOutsideFSolarTime);

%pointsInsideF = setdiff(1:length(averagedDensity), pointsOutsideF);

%densityModel(pointsInsideF) = griddata(modelAltitude, modelSolarTime, modelDensity, ...
 %   averagedAltitude(pointsInsideF), averagedSolarTime(pointsInsideF), 'v4');
%densityModel(pointsOutsideF) = F(averagedAltitude(pointsOutsideF), averagedSolarTime(pointsOutsideF));
densityModel = F(altitude, solarTime);
densityModel(densityModel < 0) = 0;

figure(5);
plot(altitude, density, '.', altitude, densityModel, 'o')
xlabel('alt / m')
ylabel('density')
title('Model prediction (o) vs real values (.)')

figure(6);
plot(solarTime, density, '.', solarTime, densityModel, 'o')
xlabel('t / h')
ylabel('density')
title('Model prediction (o) vs real values (.)')

correctedDensity = density - densityModel;

end

function [calmIndices] = findCalmPeriod(ae)
% [calmIndexes] = findCalmPeriod(ae)

rotationPeriod = 90;
calmTimeLength = 10 * rotationPeriod;
stdVector = zeros(length(ae) - calmTimeLength);
for k = 1:(length(ae) - calmTimeLength)
    stdVector(k) = std(ae(k : (k + calmTimeLength)));
end

[minVals, minIndex] = min(stdVector);
beginIndex = minIndex(1);
calmIndices = beginIndex : beginIndex + calmTimeLength;

end

function simpleDataPlot(correctedDensity, latitude, timestamps10s)
% simpleDataPlot(correctedDensity, latitude, timestamps10s)

figure(7);
timestampsInDays = timestamps10s / 60 / 60 / 24;
scatter(timestampsInDays, latitude, 60, correctedDensity, '.')
%view(2)
xlabel('t / d')
ylabel('latitude')

end

function [correctedDensity] = subtractMsisFromDensity(density, altitude, timestampsDatenum, latitude, longtitude, F107A, F107, ap, numOfModelsForOrbit)
% [correctedDensity] = subtractMsisFromDensity(density, altitude, timestampsDatenum, latitude, longtitude, F107A, F107, ap, numOfModelsForOrbit)

correctedDensity = zeros(size(density));

modelingIndices = 1:length(density);

altitudeKm = altitude / 1000;
gOverCM3ToKGOverM3 = 1000;

parfor i = modelingIndices
    model = run_nrlmsise00(altitudeKm(i), timestampsDatenum(i), latitude(i), longtitude(i), F107A, F107, ap);
    correctedDensity(i) = density(i) - gOverCM3ToKGOverM3*model(:,7);
end

end
