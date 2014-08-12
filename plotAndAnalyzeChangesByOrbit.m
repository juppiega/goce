function results = plotAndAnalyzeChangesByOrbit(firstDatenum, densityNoBg, magneticLatitude, averagedDensityNoBg, timestamps1minFixed, ...
    timestamps10s, timeOfDay, plotFigures, results)
% 

[limitedTimestamps, limitedLatitude, limitedDensity, orbits, minAllowedLatitude, maxAllowedLatitude] = ...
    giveInterpolationQualifyingOrbits(timestamps10s, magneticLatitude, densityNoBg);

results = loopAroundOrbits(orbits, timestamps1minFixed, limitedTimestamps, limitedLatitude, averagedDensityNoBg, limitedDensity, ...
    plotFigures, timeOfDay, firstDatenum, minAllowedLatitude, maxAllowedLatitude, results);

end


function [beginOrbit, beginStormOrbit, endOrbit, endColorOrbit, endColorMapOrbit, calmOrbits] = findBeginAndEndOrbits(orbits, timestamps1minFixed, limitedTimestamps, averagedDensityNoBg)
%

calmOrbits = 5;

averagedDensityNoBg = averagedDensityNoBg(ismember(timestamps1minFixed, limitedTimestamps));
timestamps1minFixed = timestamps1minFixed(ismember(timestamps1minFixed, limitedTimestamps));
[~, peakBeginIndex, peakEndIndex] = limitToNearPeak(averagedDensityNoBg, 'noSmooth', 'mean');
peakBeginIndex = find(limitedTimestamps == timestamps1minFixed(peakBeginIndex));
peakEndIndex = find(limitedTimestamps == timestamps1minFixed(peakEndIndex));

beginStormOrbit = find(orbits(:,1) <= peakBeginIndex, 1, 'last');
if beginStormOrbit > calmOrbits
    beginOrbit = beginStormOrbit - calmOrbits;
else
    beginOrbit = 1;
end

endOrbit = find(orbits(:,2) >= peakEndIndex, 1, 'first');
if endOrbit < length(orbits(:,2)) - calmOrbits
    endOrbit = endOrbit + calmOrbits;
else
    endOrbit = length(orbits(:,2));
end

maxNumOfColorOrbits = 10;
extraColorMapOrbits = 5;

if endOrbit - beginOrbit - 2 * calmOrbits > maxNumOfColorOrbits + extraColorMapOrbits
    endColorOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits - 1;
    endColorMapOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits + extraColorMapOrbits - 1;
elseif endOrbit - beginOrbit - 2 * calmOrbits > maxNumOfColorOrbits
    endColorOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits - 1;
    endColorMapOrbit = endOrbit - calmOrbits;
else
    endColorOrbit = endOrbit;
    endColorMapOrbit = endOrbit;
end

end

function [limitedTimestamps, limitedLatitude, limitedDensity, orbits, minAllowedLatitude, maxAllowedLatitude] = ...
    giveInterpolationQualifyingOrbits(timestamps10s, magneticLatitude, densityNoBg)
%

[latBeginIndex, latEndIndex] = limitLatitudeToIntegerMultipleOfOrbitalPeriod(magneticLatitude);
limitedLatitude = magneticLatitude(latBeginIndex:latEndIndex);
limitedTimestamps = timestamps10s(latBeginIndex:latEndIndex);
limitedDensity = densityNoBg(latBeginIndex:latEndIndex);
[~, exactOrbitIndices] = splitIntoOrbits(limitedLatitude, limitedTimestamps);
limitedLatitude = limitedLatitude(exactOrbitIndices);
limitedTimestamps = limitedTimestamps(exactOrbitIndices);
limitedDensity = limitedDensity(exactOrbitIndices);
[orbits, ~] = splitIntoOrbits(limitedLatitude, limitedTimestamps);

[minAllowedLatitude, maxAllowedLatitude] = findInterpolationLimits(limitedLatitude);
orbitsToDelete = find(abs(limitedLatitude(orbits(:,2)) - limitedLatitude(orbits(:,1))) < (maxAllowedLatitude - minAllowedLatitude));
newIndices = 1:length(limitedLatitude);
for i = 1:length(orbitsToDelete)
    newIndices = setdiff(newIndices, (orbits(orbitsToDelete(i), 1) : orbits(orbitsToDelete(i), 2)));
end
limitedTimestamps = limitedTimestamps(newIndices);
limitedLatitude = limitedLatitude(newIndices);
limitedDensity = limitedDensity(newIndices);
[orbits, ~] = splitIntoOrbits(limitedLatitude, limitedTimestamps);


end

function [beginOrbit, endColorMapOrbit, TADPlotIndices, TADplot, relativeResidues, calmOrbits, results] = ...
    loopAroundOrbits(orbits, timestamps1minFixed, limitedTimestamps, limitedLatitude, averagedDensityNoBg, ...
    limitedDensity, plotFigures, timeOfDay, firstDatenum, minAllowedLatitude, maxAllowedLatitude, results)
%

persistent densityByOrbitFigHandle
persistent densityByOrbitAxesHandle

if plotFigures ~= 0
    if ~isempty(strfind(lower(timeOfDay), 'morning')); densityByOrbitFigHandle = figure; subplotNum = 1; else subplotNum = 2; end
    figure(densityByOrbitFigHandle);
    densityByOrbitAxesHandle(subplotNum) = subplot(2,1,subplotNum);
    hold all;
end

linehandles = [];
relativeResidues = nan(size(limitedLatitude));
TADplot = nan(size(limitedLatitude));

[beginOrbit, beginStormOrbit, endOrbit, endColorOrbit, endColorMapOrbit, calmOrbits] = findBeginAndEndOrbits(orbits, timestamps1minFixed, limitedTimestamps, averagedDensityNoBg);

loopOrbits = beginOrbit:endOrbit;
TADPlotIndices = nan(size(limitedLatitude));
for i = 1:length(loopOrbits)
    indices = orbits(loopOrbits(i),1) : orbits(loopOrbits(i),2);
    smoothedDensity150s = smooth(limitedDensity(indices), 15);
    relativeResidues(indices) = (limitedDensity(indices) - smoothedDensity150s) ./ smoothedDensity150s;
    if (loopOrbits(i) <= endColorOrbit || loopOrbits(i) > endOrbit - calmOrbits) && plotFigures ~= 0
        h = plot(limitedLatitude(indices), smoothedDensity150s, 'LineWidth', 2);
        linehandles = [linehandles h];
    end
    
    if loopOrbits(i) <= endColorMapOrbit
        TADPlotIndices(indices) = indices;            
        smoothedDensity11800km = smooth(limitedDensity(indices), 151);
        smoothedDensity2600km = smooth(limitedDensity(indices), 33); 
        TADplot(indices) = (smoothedDensity2600km - smoothedDensity11800km) ./ smoothedDensity11800km;
    end    
end

[timeMatrix, latitudeMatrix, densityMatrix] = advancedScatteredInterpolate(limitedTimestamps, limitedLatitude, TADplot, TADPlotIndices,...
    minAllowedLatitude, maxAllowedLatitude);

estimateSpectralPowerDensity(timeMatrix, densityMatrix)

[results, magneticLatitudeForResidue, relativeResidues] = writeSmallWaveDataToResults(limitedLatitude, relativeResidues, timeOfDay, results);

if plotFigures ~= 0

    plotResiduePlot(magneticLatitudeForResidue, relativeResidues, timeOfDay)
    plotTADcolormap(limitedTimestamps, limitedLatitude, timeMatrix, latitudeMatrix, densityMatrix, minAllowedLatitude, maxAllowedLatitude, ...
        orbits, beginOrbit, calmOrbits, endColorMapOrbit, firstDatenum, timeOfDay)
    plotDensityByOrbit(linehandles, calmOrbits, timeOfDay, densityByOrbitAxesHandle)
end

end

function plotDensityByOrbit(linehandles, calmOrbits, timeOfDay, densityByOrbitAxesHandle)
%

stormOrbits = (calmOrbits + 1 : length(linehandles) - calmOrbits);
[~,cmap] = cmapline('colormap',jet,'filled', 'lines', linehandles(stormOrbits)');
colormap(cmap)
colorbar

set(linehandles(1:calmOrbits), 'LineStyle', '--')
set(linehandles(1:calmOrbits), 'Color', 'k')
endLineHandles = length(linehandles) - calmOrbits + 1 : length(linehandles);
set(linehandles(endLineHandles), 'LineStyle', '-')
set(linehandles(endLineHandles), 'Color', 'k')
hold off;
xlabel('Geomagnetic latitude')
ylabel('Density at 270 km')
xlim([-90 90])
title(timeOfDay)
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '150-s smoothed density along orbit', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
if ~isempty(strfind(lower(timeOfDay), 'evening'))
   linkaxes(densityByOrbitAxesHandle) 
end

end

function plotResiduePlot(magneticLatitudeForResidue, relativeResidues, timeOfDay)
%

persistent residueFigHandle
persistent residueAxisHandle

confIntervalStep = 10;
confIntervalX = [];
confIntervalMean = [];
confIntervalLower = [];
confIntervalUpper = [];
for i = -90:confIntervalStep:90 - confIntervalStep
    indicesInInterval = find(magneticLatitudeForResidue < i + 5 & magneticLatitudeForResidue >= i);
    if ~isempty(indicesInInterval)
        confIntervalX = [confIntervalX mean([i (i + confIntervalStep)])];
        confIntervalMean = [confIntervalMean mean(relativeResidues(indicesInInterval))];
        [~, ~, ci, ~] = ttest(relativeResidues(indicesInInterval));        
        confIntervalLower = [confIntervalLower ci(1)];
        confIntervalUpper = [confIntervalUpper ci(2)];
    end
end

if ~isempty(strfind(lower(timeOfDay), 'morning')); 
    residueFigHandle = figure('Color','w');
    residueAxisHandle = [];
    confIntervalColor = 'b'; 
    axisNum = 1;
    %mlock;
else
    confIntervalColor = 'r';
    axisNum = 2;
    %munlock;
end
figure(residueFigHandle);
residueAxisHandle = [residueAxisHandle ciplot(confIntervalLower, confIntervalUpper, confIntervalX, confIntervalColor)];
xlim([min(confIntervalX) max(confIntervalX)]);
hold all;
plot(confIntervalX, confIntervalMean, 'k--');
hold off;
xlabel('Geomagnetic latitude')
ylabel('Relative residues 95 % confidence')
title('100-1000 km variations')
if axisNum == 1
    hold on
else
    hold off
    legend(residueAxisHandle, '<-Morning', 'Evening->');
end

end


function plotTADcolormap(limitedTimestamps, limitedLatitude, timeMatrix, latitudeMatrix, densityMatrix,minAllowedLatitude, maxAllowedLatitude, orbits, ...
    beginOrbit, calmOrbits, endColorMapOrbit, firstDatenum, timeOfDay)
%

timeMatrix = timeMatrix / secondsInDay + firstDatenum;
referenceDay = datestr(min(timeMatrix(:)), 'mmmm dd, yyyy');
timeMatrix = timeMatrix - datenum(referenceDay, 'mmmm dd, yyyy');

limitedTimestamps = limitedTimestamps / secondsInDay + firstDatenum;
limitedTimestamps = limitedTimestamps - datenum(referenceDay, 'mmmm dd, yyyy');

firstOrbitTrackInd = orbits(beginOrbit + calmOrbits, 1);
lastOrbitTrackInd = orbits(endColorMapOrbit,2);
orbitTrackLatitude = limitedLatitude(firstOrbitTrackInd:lastOrbitTrackInd);
orbitTrackTime = limitedTimestamps(firstOrbitTrackInd:lastOrbitTrackInd);
orbitTrackPlotHeight = ones(size(orbitTrackLatitude)) * max(densityMatrix(:));

figure;
surf(timeMatrix, latitudeMatrix, densityMatrix,'EdgeColor', 'none');
view(2);
hold on;
h = plot3(orbitTrackTime, orbitTrackLatitude, orbitTrackPlotHeight, '.k');
set(h, 'MarkerSize', 3)
ylabel('Magnetic Latitude')
title(['1300-5200km changes (TADs) [(2600 km smooth - 10400 km smooth) / 10400 km smooth] ', timeOfDay])
colorbar
colormap jet(500)
ylim([minAllowedLatitude maxAllowedLatitude]);
xlim([min(timeMatrix(:)) max(timeMatrix(:))]);
xlabel(['Days since the UTC beginning of ', referenceDay])

end

function [timeMatrix, latitudeMatrix, densityMatrix] = advancedScatteredInterpolate(limitedTimestamps, limitedLatitude, TADplot, TADplotIndices,...
    minAllowedLatitude, maxAllowedLatitude)
%

TADplot = TADplot(~isnan(TADplot));
indices = TADplotIndices(~isnan(TADplotIndices));
TADlatitude = limitedLatitude(indices);
TADtime = limitedTimestamps(indices);

oneTenthDegreeStep = minAllowedLatitude:0.2:maxAllowedLatitude;

for i = 1:length(oneTenthDegreeStep)
    regriddedTime(:,i) = latitudeCrossingTimes(TADlatitude, TADtime, oneTenthDegreeStep(i)); 
end

regriddedDensity = interp1(TADtime, TADplot, regriddedTime, 'spline');

numOfOrbits = length(regriddedTime(:,1));
numOfValuesInOrbit = length(regriddedTime(1,:));
for i = 1:numOfValuesInOrbit
    timeThisLatitude = regriddedTime(:,i);
    densityThisLatitude = regriddedDensity(:,i);

    tInterp = interp1(1:numOfOrbits, timeThisLatitude, 1:1/60:numOfOrbits);
    interpolatedDensity = interp1(timeThisLatitude, densityThisLatitude, tInterp, 'spline');
    thisLatitude = ones(length(tInterp), 1) * oneTenthDegreeStep(i);
    latitudeMatrix(:,i) = thisLatitude;
    densityMatrix(:,i) = interpolatedDensity;
    timeMatrix(:,i) = tInterp;
end

indicesToRemove = findMatrixIndicesInDatagap(timeMatrix, TADtime);
timeMatrix(indicesToRemove) = nan(1);
latitudeMatrix(indicesToRemove) = nan(1);
densityMatrix(indicesToRemove) = 0;

end

function [results, magneticLatitudeForResidue, relativeResidues] = writeSmallWaveDataToResults(limitedLatitude, relativeResidues, timeOfDay, results)
%

magneticLatitudeForResidue = limitedLatitude(~isnan(relativeResidues));
relativeResidues = relativeResidues(~isnan(relativeResidues));

stdSouth = std(relativeResidues(magneticLatitudeForResidue < 0));
stdNorth = std(relativeResidues(magneticLatitudeForResidue > 0));

[rowNum, ~] = size(results);
emptyCells = cellfun(@isempty,results);
[~, emptyColPositions] = find(emptyCells);
colNum = min(emptyColPositions);
if rowNum == 2
    colNum = length(results(rowNum,:)) + 1;
    results{1, colNum}     = ['Std SH ', timeOfDay];
    results{1, colNum + 1} = ['Std NH ', timeOfDay];
    results{1, colNum + 2} = ['Std SH/NH ', timeOfDay];
end
results{rowNum, colNum}     = stdSouth;
results{rowNum, colNum + 1} = stdNorth;
results{rowNum, colNum + 2} = stdSouth / stdNorth;

end

function estimateSpectralPowerDensity(timeMatrix, densityMatrix)
%

samplingFreq = 6;

[rowsToRemove, ~] = find(isnan(timeMatrix));
rowsToConserve = setdiff(1:length(timeMatrix(:,1)), rowsToRemove);
densityMatrixNoNans = densityMatrix(rowsToConserve,:);
test = length(densityMatrixNoNans(:,1));
rowsToRemove = mod(test, 120);
densityMatrixNoNans(1:rowsToRemove,:) = [];

loopIteration = 1;
for i = 1:120:length(densityMatrixNoNans(:,1))
    densityMatrixForFFT(loopIteration,:) = reshape(densityMatrixNoNans(i:i + 119,:), 1, []);
    loopIteration = loopIteration + 1;
end

densityMatrixForFFT = densityMatrixForFFT - repmat(mean(densityMatrixForFFT, 2), 1, length(densityMatrixForFFT(1,:)));
numOfSamples = length(densityMatrixForFFT(1,:)) * 1;
densityFFT = abs(fft(densityMatrixForFFT(1,:)', numOfSamples));% .^ 2 ./ (numOfSamples * samplingFreq);
densityFFT = mean(densityFFT, 2);

freq = (0:numOfSamples - 1) * samplingFreq / numOfSamples;
periodInMinutes = 1 ./ freq;

earthCircumferenceAtGoceOrbit = 2 * pi * (6371 + 270);
goceOrbitPeriod = 89;
speed = earthCircumferenceAtGoceOrbit / goceOrbitPeriod;
periodInKm = speed * periodInMinutes;

figure;
plot(1:length(densityMatrixForFFT(1,:)), densityMatrixForFFT(1,:));

minimumResolution = 2 * min(periodInKm);
maximumResolution = 16000;

%densityFFT(periodInKm < minimumResolution | periodInKm > maximumResolution) = 0;

%densityFFT = smooth(densityFFT,5);
%fullPower = trapz(densityFFT);
relativeSpectralDensity = densityFFT;% / fullPower;

figure;
plot(periodInKm, relativeSpectralDensity,'.');
xlim([minimumResolution maximumResolution])

end