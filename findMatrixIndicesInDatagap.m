function indicesToRemove = findMatrixIndicesInDatagap(timeMatrix, timestamps)
%

hourInSeconds = 60 * 60;
gapBeginIndices = find(diff(timestamps) > hourInSeconds);
gapBeginTimes = timestamps(gapBeginIndices);
gapEndTimes = timestamps(gapBeginIndices + 1);

indicesToRemove = [];
for i = 1:length(gapBeginTimes)
    indicesInGap = find(timeMatrix > gapBeginTimes(i) & timeMatrix < gapEndTimes(i));
    indicesToRemove = [indicesToRemove; indicesInGap];
end

end