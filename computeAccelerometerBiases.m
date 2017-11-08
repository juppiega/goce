function [goceScale, graceScale] = computeAccelerometerBiases(rhoStruct)

% rhoStruct = removeDataPoints(rhoStruct, rhoStruct.timestamps >= datenum('2010-01-01') & ...
%     ismember((1:length(rhoStruct.data))', rhoStruct.grace), ...
%     false,true,false,false);

rhoStruct = removeDataPoints(rhoStruct, abs(rhoStruct.latitude) > 60 | rhoStruct.timestamps > datenum('2006-01-01'),...
    true, true, true,true);

reducedLst = rhoStruct.solarTime;
reducedLst(reducedLst >= 12) = reducedLst(reducedLst >= 12) - 12;
goceLst = reducedLst(rhoStruct.goce);
champLst = reducedLst(rhoStruct.champ);
graceLst = reducedLst(rhoStruct.grace);
t0 = rhoStruct.timestamps(rhoStruct.champ(1));
goceTimes = round((rhoStruct.timestamps(rhoStruct.goce) - t0)*86400);
champTimes = round((rhoStruct.timestamps(rhoStruct.champ) - t0)*86400);
graceTimes = round((rhoStruct.timestamps(rhoStruct.grace) - t0)*86400);

% goceInd = find(ismember(goceTimes, champTimes));
% graceInd = find(ismember(graceTimes, champTimes));
% champGoce = find(ismember(champTimes, goceTimes));
% champGrace = find(ismember(champTimes, graceTimes));
% 
% goceChampNonCoplanar = find(abs(goceLst(goceInd) - champLst(champGoce)) > 1);
% graceChampNonCoplanar = find(abs(graceLst(graceInd) - champLst(champGrace)) > 1);
% 
% goceInd(goceChampNonCoplanar) = [];
% graceInd(graceChampNonCoplanar) = [];
% champGoce(goceChampNonCoplanar) = [];
% champGrace(graceChampNonCoplanar) = [];
% 
% goceInd = goceInd(1:10:end);
% graceInd = graceInd(1:10:end);
% champGoce = champGoce(1:10:end);
% champGrace = champGrace(1:10:end);

if isempty(goceTimes)
    goceInd = []; champGoce = [];
else
    [goceInd, champGoce] = computeCoplanarInd(goceTimes, goceLst, champTimes, champLst);
end
[graceInd, champGrace] = computeCoplanarInd(graceTimes, graceLst, champTimes, champLst);

N = length(rhoStruct.data);
%goceRmInd = setdiff(1:N, goceInd + rhoStruct.goce(1) - 1);
graceRmInd = setdiff(1:N, graceInd + rhoStruct.grace(1) - 1);
%champGoceRmInd = setdiff(1:N, champGoce + rhoStruct.champ(1) - 1);
champGraceRmInd = setdiff(1:N, champGrace + rhoStruct.champ(1) - 1);

%goceStruct = removeDataPoints(rhoStruct, goceRmInd);
graceStruct = removeDataPoints(rhoStruct, graceRmInd);
%champGoceStruct = removeDataPoints(rhoStruct, champGoceRmInd);
champGraceStruct = removeDataPoints(rhoStruct, champGraceRmInd);

model = 'dtm';
if strcmpi(model,'msis')
    %[~,msisGoce] = computeMsis(goceStruct);
    [~,msisGrace] = computeMsis(graceStruct);
    %[~,msisChampGoce] = computeMsis(champGoceStruct);
    [~,msisChampGrace] = computeMsis(champGraceStruct);
else
    %[~,msisGoce] = computeDtm(goceStruct);
    [~,msisGrace] = computeDtm(graceStruct);
    %[~,msisChampGoce] = computeDtm(champGoceStruct);
    [~,msisChampGrace] = computeDtm(champGraceStruct);

end

%goceMsisBias = mean(goceStruct.data ./ msisGoce);
graceMsisBias = mean(graceStruct.data ./ msisGrace);
%champGoceMsisBias= mean(champGoceStruct.data ./ msisChampGoce);
champGraceMsisBias = mean(champGraceStruct.data ./ msisChampGrace);

%goceScale = champGoceMsisBias ./ goceMsisBias;
goceScale = 1.0;
graceScale = champGraceMsisBias ./ graceMsisBias;

save('satScales.mat', 'goceScale')
save('satScales.mat', 'graceScale', '-append')

end

function [satInd, otherSatInd] = computeCoplanarInd(satTimes, satLst, otherTimes, otherLst)

maxLstDiff = 0.5;

otherLstAtSatTimes = interp1(otherTimes, otherLst, satTimes);
satLstAtOtherTimes = interp1(satTimes, satLst, otherTimes);

satInd = find(abs(satLst - otherLstAtSatTimes) <= maxLstDiff);
otherSatInd = find(abs(otherLst - satLstAtOtherTimes) <= maxLstDiff);

end