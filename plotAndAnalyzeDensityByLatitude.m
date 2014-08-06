function [results] = plotAndAnalyzeDensityByLatitude(firstDatenum, ae, timestamps1min, timestamps1minFixed, correctedDensity, msisDensity,...
    timestamps10s, magneticLatitude, timeOfDay, plotFigures, results)
% plotCorrectedDensityLatitudes(ae, timestamps1min, correctedDensity, timestamps10s, latitude, timestampsDatenum, computeLatitudes);
    
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

    plotDensityLatitudeTimeSurf(firstDatenum, magneticLatitude, timestamps10s, latitudeMatrix, timeMatrix, goceDensityMatrix, msisDensityMatrix, timeOfDay);
end

limitedTimestamps = limitedTimestamps(ismember(limitedTimestamps, timestamps1minFixed));

aeShort = ae(ismember(timestamps1min, limitedTimestamps));
if plotFigures ~= 0
    plotPeakAnalysis(limitedTimestamps, aeShort, densityByLatitude, crossingTimes, minAllowedLatitude, maxAllowedLatitude, timeOfDay);
end

results = plotAndAnalyzeByHemisphere(firstDatenum, densityByLatitude, ae, crossingTimes, timestamps1min, timeOfDay, oneDegreeStep, plotFigures, results);

end
