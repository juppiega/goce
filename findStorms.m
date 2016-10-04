function [stormBeginInd, stormEndInd, combinedInd, satInfo] = findStorms(rhoStruct, dstOrAE, threshold)

goceStruct = removeDataPoints(rhoStruct, [rhoStruct.champ'; rhoStruct.grace'; rhoStruct.swarm']);
champStruct = removeDataPoints(rhoStruct, [rhoStruct.goce'; rhoStruct.grace'; rhoStruct.swarm']);
graceStruct = removeDataPoints(rhoStruct, [rhoStruct.champ'; rhoStruct.goce'; rhoStruct.swarm']);
swarmStruct = removeDataPoints(rhoStruct, [rhoStruct.champ'; rhoStruct.goce'; rhoStruct.grace']);

[goceBI, goceEI, goceCI] = findStormsForSat(goceStruct, dstOrAE, threshold);
[champBI, champEI, champCI] = findStormsForSat(champStruct, dstOrAE, threshold);
[graceBI, graceEI, graceCI] = findStormsForSat(graceStruct, dstOrAE, threshold);
[swarmBI, swarmEI, swarmCI] = findStormsForSat(swarmStruct, dstOrAE, threshold);

satInfo(1:length(goceBI)) = 0; k = length(goceBI);
satInfo(k+1:k+length(champBI)) = 1; k = k + length(champBI);
satInfo(k+1:k+length(graceBI)) = 2; k = k + length(graceBI);
satInfo(k+1:k+length(swarmBI)) = 3;

champBegin = rhoStruct.champ(1)-1;
graceBegin = rhoStruct.grace(1)-1;
goceBegin = rhoStruct.goce(1)-1;
swarmBegin = rhoStruct.swarm(1)-1;

stormBeginInd(1:length(goceBI)) = goceBI+goceBegin; k = length(goceBI);
stormBeginInd(k+1:k+length(champBI)) = champBI+champBegin; k = k + length(champBI);
stormBeginInd(k+1:k+length(graceBI)) = graceBI+graceBegin; k = k + length(graceBI);
stormBeginInd(k+1:k+length(swarmBI)) = swarmBI+swarmBegin;

stormEndInd(1:length(goceBI)) = goceEI+goceBegin; k = length(goceBI);
stormEndInd(k+1:k+length(champBI)) = champEI+champBegin; k = k + length(champBI);
stormEndInd(k+1:k+length(graceBI)) = graceEI+graceBegin; k = k + length(graceBI);
stormEndInd(k+1:k+length(swarmBI)) = swarmEI+swarmBegin;

combinedInd = [goceCI+goceBegin; champCI+champBegin; graceCI+graceBegin; swarmCI+swarmBegin];

end

function [stormBeginInd, stormEndInd, combinedInd] = findStormsForSat(rhoStruct, dstOrAE, threshold)

if strcmpi(dstOrAE, 'Dst')
    thresholdInd = find(rhoStruct.dst <= threshold);
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

stormBegin = stormBegin - 1;
stormEnd = stormEnd + 1;

stormBeginInd = findNearestIndices(stormBegin, rhoStruct.timestamps);
stormEndInd = findNearestIndices(stormEnd, rhoStruct.timestamps);

if nargout == 3
    k = 0;
    intervalLengths = stormEndInd-stormBeginInd+1;
    combinedInd = zeros(sum(intervalLengths), 1);
    for i = 1:length(stormBeginInd)
        combinedInd(k+1:k+intervalLengths(i)) = stormBeginInd(i):stormEndInd(i);
        k = k + intervalLengths(i);
    end
end

end

function ind = findNearestIndices(singleTimes, timeVector)

ind = zeros(size(singleTimes));
for i = 1:length(singleTimes)
    tempInd = find(timeVector >= singleTimes(i), 1);
    if isempty(tempInd); tempInd = length(timeVector); end
    if tempInd == 1
        a = 1;
    end
    ind(i) = tempInd;
end

end