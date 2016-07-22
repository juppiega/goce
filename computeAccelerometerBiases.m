function [goceScale, graceScale] = computeAccelerometerBiases(rhoStruct)

reducedLst = rhoStruct.solarTime;
reducedLst(reducedLst >= 12) = reducedLst(reducedLst >= 12) - 12;
goceLst = reducedLst(rhoStruct.goce);
champLst = reducedLst(rhoStruct.champ);
graceLst = reducedLst(rhoStruct.grace);
t0 = rhoStruct.timestamps(rhoStruct.goce(1));
goceTimes = round((rhoStruct.timestamps(rhoStruct.goce) - t0)*86400);
champTimes = round((rhoStruct.timestamps(rhoStruct.champ) - t0)*86400);
graceTimes = round((rhoStruct.timestamps(rhoStruct.grace) - t0)*86400);

goceInd = find(ismember(goceTimes, champTimes));
graceInd = find(ismember(graceTimes, champTimes));
champGoce = find(ismember(champTimes, goceTimes));
champGrace = find(ismember(champTimes, graceTimes));

goceChampNonCoplanar = find(abs(goceLst(goceInd) - champLst(champGoce)) > 1);
graceChampNonCoplanar = find(abs(graceLst(graceInd) - champLst(champGrace)) > 1);

goceInd(goceChampNonCoplanar) = [];
graceInd(graceChampNonCoplanar) = [];
champGoce(goceChampNonCoplanar) = [];
champGrace(graceChampNonCoplanar) = [];

goceStruct = removeDataPoints(rhoStruct, goceInd + rhoStruct.goce(1) - 1);
graceStruct = removeDataPoints(rhoStruct, graceInd + rhoStruct.grace(1) - 1);
champGoceStruct = removeDataPoints(rhoStruct, champGoce + rhoStruct.champ(1) - 1);
champGraceStruct = removeDataPoints(rhoStruct, champGrace + rhoStruct.champ(1) - 1);

[~,msisGoce] = computeMsis(goceStruct);
[~,msisGrace] = computeMsis(graceStruct);
[~,msisChampGoce] = computeMsis(champGoceStruct);
[~,msisChampGrace] = computeMsis(champGraceStruct);

goceMsisBias = mean(goceStruct.data ./ msisGoce);
graceMsisBias = mean(graceStruct.data ./ msisGrace);
champGoceMsisBias= mean(champGoceStruct.data ./ msisChampGoce);
champGraceMsisBias = mean(champGraceStruct.data ./ msisChampGrace);

goceScale = champGoceMsisBias ./ goceMsisBias;
graceScale = champGraceMsisBias ./ graceMsisBias;

save('satScales.mat', 'goceScale')
save('satScales.mat', 'graceScale', '-append')

end