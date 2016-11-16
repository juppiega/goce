load ilData originalRhoStruct

% % Remove unnecessary observations.
t0 = datenum('2003-10-20');
t1 = datenum('2003-11-05');
windowLen = 6;

removeInd = ~ismember(1:length(originalRhoStruct.data), originalRhoStruct.champ);
removeTimes = originalRhoStruct.timestamps < t0 - windowLen/24/2 | ...
                originalRhoStruct.timestamps > t1;
removeInd(removeTimes) = true;
assimiStruct = removeDataPoints(originalRhoStruct, removeInd, false,true,true,true);

aver = averageRho(assimiStruct,true);