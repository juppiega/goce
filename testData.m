function [] = testData (doy1, doy2, satellite)


load('ilData.mat','rhoStruct')
rhoStruct = removeAndFixData(rhoStruct,0);
N = length(rhoStruct.data);
if strcmpi(satellite,'goce')
    removeInd = ~ismember(1:N, rhoStruct.goce);
elseif strcmpi(satellite,'champ')
    removeInd = ~ismember(1:N, rhoStruct.champ);
elseif strcmpi(satellite,'grace')
    removeInd = ~ismember(1:N, rhoStruct.grace);
elseif strcmpi(satellite,'swarm')
    removeInd = ~ismember(1:N, rhoStruct.swarm);
elseif strcmpi(satellite,'all')
    removeInd = false(N,1);
elseif strcmpi(satellite,'allNoSwarm')
    removeInd = ismember(1:N, rhoStruct.swarm);
end
rhoStruct = removeDataPoints(rhoStruct, removeInd, true, true, true, true);
[yr,~,~,~,~,~] = datevec(rhoStruct.timestamps);
yearVec = [yr, repmat([1,1,0,0,0], length(yr), 1)];
rhoStruct.doy = rhoStruct.timestamps - datenum(yearVec) + 1;

%ind = doy1 < rhoStruct.doy & rhoStruct.doy < doy2;
%rhoStruct = removeDataPoints(rhoStruct, ind,true,true,true,true);

load('ilComparison.mat','ilRho')
load('msisDtmComparison.mat')
ilRho(removeInd) = [];
msisRho(removeInd) = [];
dtmRho(removeInd) = [];

dt = 2;
k = 1;
for doy = doy1:dt:doy2
    ind = doy < rhoStruct.doy & rhoStruct.doy < doy + dt;
    RMS_il(k)  = rms(rhoStruct.data(ind) ./ ilRho(ind) - 1);
    RMS_dtm(k) = rms(rhoStruct.data(ind) ./ dtmRho(ind) - 1);
    numObs(k) = sum(ind);
    k = k + 1;
end

figure;
[hAx, hIL, hCount] = plotyy(doy1:dt:doy2, RMS_il, doy1:dt:doy2, numObs);
set(hCount, 'color', 'k', 'linestyle', '-');
ylim(hAx(1), [0 1]);
hold on;
plot(hAx(1), doy1:dt:doy2, RMS_dtm);
legend(hAx(1),'il','dtm')

end
