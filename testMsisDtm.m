clear S
load('ilData.mat','HeStruct');
rmInd = ~ismember(1:length(HeStruct.data),HeStruct.de2);
S = removeDataPoints(HeStruct,rmInd,true,true,true,true);

N = length(S.data);
msisVals = zeros(N,5);
[~, ~, msisVals(:,1), msisVals(:,2), msisVals(:,3), msisVals(:,4), msisVals(:,5)] = computeMsis(S);
msis_OM = mean(S.data ./ msisVals(:,3))
hist(S.data ./ msisVals(:,3),20)

dtmVals = zeros(N,4);
[~, ~, dtmVals(:,1), dtmVals(:,2), dtmVals(:,3), dtmVals(:,4)] = computeDtm(S);
dtm_OM = mean(S.data ./ dtmVals(:,3))