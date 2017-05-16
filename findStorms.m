function [stormBeginInd, stormEndInd, combinedInd, satInfo] = findStorms(rhoStruct, dstOrAE, threshold)

goceStruct = removeDataPoints(rhoStruct, [rhoStruct.champ'; rhoStruct.grace'; rhoStruct.swarm']);
champStruct = removeDataPoints(rhoStruct, [rhoStruct.goce'; rhoStruct.grace'; rhoStruct.swarm']);
graceStruct = removeDataPoints(rhoStruct, [rhoStruct.champ'; rhoStruct.goce'; rhoStruct.swarm']);
swarmStruct = removeDataPoints(rhoStruct, [rhoStruct.champ'; rhoStruct.goce'; rhoStruct.grace']);

if length(goceStruct.data) > 0; [goceBI, goceEI, goceCI] = findStormsForSat(goceStruct, dstOrAE, threshold); else goceBI = []; end
if length(champStruct.data) > 0; [champBI, champEI, champCI] = findStormsForSat(champStruct, dstOrAE, threshold); else champBI = []; end
if length(graceStruct.data) > 0; [graceBI, graceEI, graceCI] = findStormsForSat(graceStruct, dstOrAE, threshold); else graceBI = []; end
if length(swarmStruct.data) > 0; [swarmBI, swarmEI, swarmCI] = findStormsForSat(swarmStruct, dstOrAE, threshold); else swarmBI = []; end

satInfo(1:length(goceBI)) = 0; k = length(goceBI);
satInfo(k+1:k+length(champBI)) = 1; k = k + length(champBI);
satInfo(k+1:k+length(graceBI)) = 2; k = k + length(graceBI);
satInfo(k+1:k+length(swarmBI)) = 3;

if length(champStruct.data) > 0; champBegin = rhoStruct.champ(1)-1; end
if length(graceStruct.data) > 0; graceBegin = rhoStruct.grace(1)-1; end
if length(goceStruct.data) > 0; goceBegin = rhoStruct.goce(1)-1; end
if length(swarmStruct.data) > 0; swarmBegin = rhoStruct.swarm(1)-1; end

if length(goceStruct.data) > 0; stormBeginInd(1:length(goceBI)) = goceBI+goceBegin; k = length(goceBI); end
if length(champStruct.data) > 0; stormBeginInd(k+1:k+length(champBI)) = champBI+champBegin; k = k + length(champBI); end
if length(graceStruct.data) > 0; stormBeginInd(k+1:k+length(graceBI)) = graceBI+graceBegin; k = k + length(graceBI); end
if length(swarmStruct.data) > 0; stormBeginInd(k+1:k+length(swarmBI)) = swarmBI+swarmBegin; end

if length(goceStruct.data) > 0;stormEndInd(1:length(goceBI)) = goceEI+goceBegin; k = length(goceBI);end
if length(champStruct.data) > 0;stormEndInd(k+1:k+length(champBI)) = champEI+champBegin; k = k + length(champBI);end
if length(graceStruct.data) > 0;stormEndInd(k+1:k+length(graceBI)) = graceEI+graceBegin; k = k + length(graceBI);end
if length(swarmStruct.data) > 0;stormEndInd(k+1:k+length(swarmBI)) = swarmEI+swarmBegin;end

combinedInd = [];
if length(goceStruct.data) > 0; combinedInd = [combinedInd; goceCI+goceBegin];end
if length(champStruct.data) > 0; combinedInd = [combinedInd; champCI+champBegin];end
if length(graceStruct.data) > 0; combinedInd = [combinedInd; graceCI+graceBegin];end
if length(swarmStruct.data) > 0; combinedInd = [combinedInd; swarmCI+swarmBegin];end

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