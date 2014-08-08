function indicesToRemove = findMatrixIndicesInDatagap(timeMatrix, timestamps)
%

hourInSeconds = 60 * 60;
gapBeginIndices = find(diff(timestamps) > hourInSeconds);
gapBeginTimes = timestamps(gapBeginIndices);
gapEndTimes = timestamps(gapBeginIndices + 1);

indicesToRemove = zeros(size(timeMatrix));
for i = 1:length(gapBeginTimes)
    indicesInGap = timeMatrix > gapBeginTimes(i) & timeMatrix < gapEndTimes(i);
    indicesToRemove(indicesInGap) = 1;
end

indicesToRemove = logical(indicesToRemove);

end