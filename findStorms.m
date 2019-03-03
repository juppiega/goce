function [stormBeginInd, stormEndInd, combinedInd, satInfo] = findStorms(rhoStruct, dstOrAE, threshold)

goceStruct = removeDataPoints(rhoStruct, [rhoStruct.champ'; rhoStruct.grace'; rhoStruct.swarm'], true, true, true, true);
champStruct = removeDataPoints(rhoStruct, [rhoStruct.goce'; rhoStruct.grace'; rhoStruct.swarm'], true, true, true, true);
graceStruct = removeDataPoints(rhoStruct, [rhoStruct.champ'; rhoStruct.goce'; rhoStruct.swarm'], true, true, true, true);
swarmStruct = removeDataPoints(rhoStruct, [rhoStruct.champ'; rhoStruct.goce'; rhoStruct.grace'], true, true, true, true);

if length(goceStruct.data) > 0; [goceBI, goceEI, goceCI] = findStormsForSat(goceStruct, dstOrAE, threshold,0,2); else goceBI = []; end
if length(champStruct.data) > 0; [champBI, champEI, champCI] = findStormsForSat(champStruct, dstOrAE, threshold,0,2); else champBI = []; end
if length(graceStruct.data) > 0; [graceBI, graceEI, graceCI] = findStormsForSat(graceStruct, dstOrAE, threshold,0,2); else graceBI = []; end
if length(swarmStruct.data) > 0; [swarmBI, swarmEI, swarmCI] = findStormsForSat(swarmStruct, dstOrAE, threshold,0,2); else swarmBI = []; end

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

