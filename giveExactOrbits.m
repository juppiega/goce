function [limitedTimestamps, limitedLatitude, minAllowedLatitude, maxAllowedLatitude] = giveExactOrbits(timestamps10s, magneticLatitude, removeOrbits)
%

[latBeginIndex, latEndIndex] = limitLatitudeToIntegerMultipleOfOrbitalPeriod(magneticLatitude);
limitedLatitude = magneticLatitude(latBeginIndex:latEndIndex);
limitedTimestamps = timestamps10s(latBeginIndex:latEndIndex);
[~, exactOrbitIndices] = splitIntoOrbits(limitedLatitude, limitedTimestamps);
limitedLatitude = limitedLatitude(exactOrbitIndices);
limitedTimestamps = limitedTimestamps(exactOrbitIndices);
[orbits, ~] = splitIntoOrbits(limitedLatitude, limitedTimestamps);

[minAllowedLatitude, maxAllowedLatitude] = findInterpolationLimits(limitedLatitude);

%if nargin < 3
    if satelliteIsGoingSouth(limitedLatitude)
        orbitsEndingTooNorth = find(limitedLatitude(orbits(:,2)) > minAllowedLatitude);
        orbitsBeginningTooSouth = find(limitedLatitude(orbits(:,1)) < maxAllowedLatitude);
        orbitsToDelete = unique(vertcat(orbitsBeginningTooSouth, orbitsEndingTooNorth));
    else
        orbitsEndingTooSouth = find(limitedLatitude(orbits(:,2)) < maxAllowedLatitude);
        orbitsBeginningTooNorth = find(limitedLatitude(orbits(:,1)) > minAllowedLatitude);
        orbitsToDelete = unique(vertcat(orbitsBeginningTooNorth, orbitsEndingTooSouth));
    end
    newIndices = 1:length(limitedLatitude);
    for i = 1:length(orbitsToDelete)
        newIndices = setdiff(newIndices, (orbits(orbitsToDelete(i), 1) : orbits(orbitsToDelete(i), 2)));
    end
    limitedTimestamps = limitedTimestamps(newIndices);
    limitedLatitude = limitedLatitude(newIndices);
%end

end