function [] = testData (doy1, doy2, satellite)


load('ilData.mat','rhoStruct')
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

load('ilComparison.met','ilRho')
%ilRho(ind) = [];

dt = 2;
k = 1;
for doy = doy1:dt:doy2
    ind = doy < rhoStruct.doy & rhoStruct.doy < doy + dt;
    RMS(k) = rms(rhoStruct.data(ind) ./ ilRho(ind) - 1);
    k = k + 1;
end

figure;
plot(doy1:dt:doy2, RMS);

end
