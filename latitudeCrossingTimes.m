function [crossingTimes] = latitudeCrossingTimes(limitedLatitude, limitedTimestamps, computeLatitude)
% latitudeCrossingTimes(limitedLatitude, limitedTimestamps, computeLatitude)

subtractedLatitude = limitedLatitude - computeLatitude;
crossingIndices = find(subtractedLatitude(1:end-1) .* subtractedLatitude(2:end) < 0);
crossingIndices = crossingIndices(abs(limitedLatitude(crossingIndices) - limitedLatitude(crossingIndices + 1)) < 20);
[crossingIndices] = checkConsistency(crossingIndices);
tBeforeCrossing = limitedTimestamps(crossingIndices);
tAfterCrossing = limitedTimestamps(crossingIndices + 1);
subLatBeforeCrossing = abs(subtractedLatitude(crossingIndices));
subLatAfterCrossing = abs(subtractedLatitude(crossingIndices + 1));
crossingTimes = (subLatAfterCrossing .* tBeforeCrossing + subLatBeforeCrossing .* tAfterCrossing) ...
                                    ./ (subLatAfterCrossing + subLatBeforeCrossing);
end

function [crossingIndices] = checkConsistency(crossingIndices)
% [crossingIndices] = checkConsistency(crossingIndices)

diffOfIndices = diff(crossingIndices);
inconsistentIndices = find(diffOfIndices > 1.5 * mean(diffOfIndices));
incIndicesApproxVals = round((crossingIndices(inconsistentIndices) + crossingIndices(inconsistentIndices + 1)) ./ 2);
for i = 1:length(inconsistentIndices)
   crossingIndices = [crossingIndices(1:inconsistentIndices(i)); incIndicesApproxVals(i); ...
                        crossingIndices(inconsistentIndices(i) + 1:end)]; 
   inconsistentIndices = inconsistentIndices + 1;
end
end