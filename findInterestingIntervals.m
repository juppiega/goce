function intervalsOfInterest = findInterestingIntervals(ae, timestampsAeDatenum, epsilonQualityFlag, timestampsEpsilonDatenum, densityNoBg, timestampsDensityDatenum, threshold)
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

indicesInDay = 24 * 60;
if length(intervalsOfInterest(:,1)) > 1
    lengthOfOverlap = 0;
    newIntervals = intervalsOfInterest;
    indicesToRemove = [];
    for i = 1:length(intervalsOfInterest(:,1)) - 1
        nextStormOverlap = intervalsOfInterest(i, 2) - intervalsOfInterest(i + 1,1);
        nextStormLength = intervalsOfInterest(i + 1, 2) - intervalsOfInterest(i + 1,1);
        if nextStormOverlap / nextStormLength > 0.4
            if lengthOfOverlap > 0
                newIntervals(i - lengthOfOverlap,:) = [newIntervals(i - lengthOfOverlap, 1) intervalsOfInterest(i + 1, 2)];
            else
                newIntervals(i,:) = [intervalsOfInterest(i, 1) intervalsOfInterest(i + 1, 2)];
            end
            indicesToRemove = [indicesToRemove; (i + 1)];
            lengthOfOverlap = lengthOfOverlap + 1;
        else
            lengthOfOverlap = 0;
        end
    end
    indicesToConserve = setdiff(1:length(newIntervals(:,1)), indicesToRemove);
    intervalsOfInterest = newIntervals(indicesToConserve,:);
end

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
    
    valuesShouldHaveAfter = length(ae(stormEndIndex:aeIndices(end)));
    valuesShouldHaveBefore = length(ae(aeIndices(1):stormBeginIndex));
    valuesShouldHaveStorm = length(ae(stormBeginIndex:stormEndIndex));
    
    if densityStormIndices / valuesShouldHaveStorm < 0.5 ||...
       (densityAfterIndices / valuesShouldHaveAfter < 0.75 && densityBeforeIndices / valuesShouldHaveBefore < 0.75)
        intervalsToRemove = [intervalsToRemove i];
        fprintf(2, '%s\n', ['Warning: Storm between dates ', datestr(timestampsAeDatenum(aeIndices(1)), 'yyyy-mm-dd'), ...
            ' and ', datestr(timestampsAeDatenum(aeIndices(end)), 'yyyy-mm-dd'), ' has too large density data gaps. It will be omitted.']);
    end
    
    if epsilonStormIndices / valuesShouldHaveStorm < 0.5 ||...
       (epsilonAfterIndices / valuesShouldHaveAfter < 0.3 && epsilonBeforeIndices / valuesShouldHaveBefore < 0.3)
        intervalsToRemove = [intervalsToRemove i];
        fprintf(2, '%s\n', ['Warning: Storm between dates ', datestr(timestampsAeDatenum(aeIndices(1)), 'yyyy-mm-dd'), ...
            ' and ', datestr(timestampsAeDatenum(aeIndices(end)), 'yyyy-mm-dd'), ' has too large Akasofu epsilon data gaps. It will be omitted.']);
    end
      
end

intervalsToConserve = setdiff(1:length(intervalsOfInterest(:,1)), intervalsToRemove);
intervalsOfInterest = intervalsOfInterest(intervalsToConserve,:);

fprintf('%s\n', 'Following storms will be analyzed:')
for i = 1:length(intervalsOfInterest(:,1))
    beginIndex = intervalsOfInterest(i,1);
    endIndex = intervalsOfInterest(i,2);
    fprintf('%s\n', [datestr(timestampsAeDatenum(beginIndex), 'yyyy-mm-dd'), ' to ', datestr(timestampsAeDatenum(endIndex), 'yyyy-mm-dd')])
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