function [stormBeginInd, stormEndInd, combinedInd] = findStormsForSat(rhoStruct, dstOrAE, threshold, offset, afterDays, uniqueInd)

satInd = zeros(size(rhoStruct.data));

if isfield(rhoStruct,'goce') satInd(rhoStruct.goce) = 1; end
if isfield(rhoStruct,'champ') satInd(rhoStruct.champ) = 2; end
if isfield(rhoStruct,'grace') satInd(rhoStruct.grace) = 3; end
if isfield(rhoStruct,'swarm') satInd(rhoStruct.swarm) = 4; end
if isfield(rhoStruct,'de2') satInd(rhoStruct.de2) = 5; end
if isfield(rhoStruct,'aeC') satInd(rhoStruct.aeC) = 6; end
if isfield(rhoStruct,'aeE') satInd(rhoStruct.aeE) = 7; end
if isfield(rhoStruct,'aeENace') satInd(rhoStruct.aeENace) = 8; end
if isfield(rhoStruct,'aeEOss') satInd(rhoStruct.aeEOss) = 9; end
if isfield(rhoStruct,'guvi') satInd(rhoStruct.guvi) = 10; end
if isfield(rhoStruct,'aeros') satInd(rhoStruct.aeros) = 11; end
if isfield(rhoStruct,'aeCOss') satInd(rhoStruct.aeCOss) = 12; end

rhoStruct.satInd = satInd;
if ~all(satInd == 0) && any(satInd == 0)
    error('Something went wrong!')
end

if strcmpi(dstOrAE, 'Dst')
    thresholdInd = find(rhoStruct.dst <= threshold);
elseif strcmpi(dstOrAE, 'AE')
    thresholdInd = find(any(rhoStruct.aeInt(:,4:18) >= threshold, 2));
elseif strcmpi(dstOrAE, 'AE2h')
    thresholdInd = find(rhoStruct.aeInt(:,1) >= threshold);
elseif strcmpi(dstOrAE, 'AE4h')
    thresholdInd = find(rhoStruct.aeInt(:,2) >= threshold);
elseif strcmpi(dstOrAE, 'AE8h')
    thresholdInd = find(rhoStruct.aeInt(:,3) >= threshold);
elseif strcmpi(dstOrAE, 'AE16h')
    thresholdInd = find(rhoStruct.aeInt(:,4) >= threshold);
elseif strcmpi(dstOrAE, 'AE30h')
    thresholdInd = find(rhoStruct.aeInt(:,6) >= threshold);
elseif strcmpi(dstOrAE, 'AE40h')
    thresholdInd = find(rhoStruct.aeInt(:,7) >= threshold);
end
stormTimes = rhoStruct.timestamps(thresholdInd);
stormChangeInd = find(diff((stormTimes)) > 2);
stormBegin = [stormTimes(1); rhoStruct.timestamps(thresholdInd(stormChangeInd+1))];
stormEnd = [rhoStruct.timestamps(thresholdInd(stormChangeInd)); stormTimes(end)];

stormSat = rhoStruct.satInd(thresholdInd);
stormBeginSat = [stormSat(1); rhoStruct.satInd(thresholdInd(stormChangeInd+1))];
stormEndSat = [rhoStruct.satInd(thresholdInd(stormChangeInd)); stormSat(end)];

rmInd = stormBeginSat ~= stormEndSat;
stormBegin(rmInd) = [];
stormEnd(rmInd) = [];
stormBeginSat(rmInd) = [];
stormSat = stormBeginSat;

if strcmpi(dstOrAE, 'Dst')    
    stormBegin = stormBegin - 2 - 2;
else
    stormBegin = stormBegin - 1;
end
if nargin >= 5
    stormEnd = stormEnd + afterDays;
else
    stormEnd = stormEnd + 1;
end

stormBeginInd = findNearestIndices(stormBegin, stormSat, rhoStruct);
stormEndInd = findNearestIndices(stormEnd, stormSat, rhoStruct);

intervalLengths = stormEndInd-stormBeginInd+1;
rmInd = intervalLengths <= 1;
stormBeginInd(rmInd) = [];
stormEndInd(rmInd) = [];
intervalLengths(rmInd) = [];
if nargout == 3
    k = 0;
    combinedInd = zeros(sum(intervalLengths), 1);
    for i = 1:length(stormBeginInd)
        combinedInd(k+1:k+intervalLengths(i)) = stormBeginInd(i):stormEndInd(i);
        k = k + intervalLengths(i);
    end
end

if nargin >= 4
    stormBeginInd = stormBeginInd + offset;
    stormEndInd = stormEndInd + offset;
    combinedInd = combinedInd + offset;
end

if nargin >= 6 && uniqueInd
    combinedInd = unique(combinedInd);
    ind = false(size(rhoStruct.data));
    ind(combinedInd) = true;
    combinedInd = ind;
end

end

function ind = findNearestIndices(singleTimes, stormSat, rhoStruct)

ind = zeros(size(singleTimes));
for i = 1:length(singleTimes)
    tempInd = find(rhoStruct.timestamps >= singleTimes(i) & rhoStruct.satInd == stormSat(i), 1);
    if isempty(tempInd); tempInd = length(rhoStruct.timestamps); end
    ind(i) = tempInd;
end

end