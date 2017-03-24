function [] = plotRadarData(radarNum)

if exist('ilData.mat', 'file')
    load ilData.mat
else
    error('File ilData.mat not found!')
end


[rhoStruct, TempStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct] = ...
    removeAndFixData(rhoStruct, 0, TempStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct, true);

removeInd = lbT0Struct.index ~= radarNum;
lbT0Struct = removeDataPoints(lbT0Struct, removeInd);
lbT0Struct = computeVariablesForFit(lbT0Struct);

dt = 2;
meanData = [];
stdData = [];
solarTime = dt/2 : dt : 24-dt/2;
for t = 0:dt:24-dt
    ind = t <= lbT0Struct.solarTime & lbT0Struct.solarTime < t+dt;
    meanData = [meanData; mean(lbT0Struct.data(ind))];
    stdData = [stdData; std(lbT0Struct.data(ind))];
end

figure;
errorbar(solarTime, meanData, stdData, 's');
xlabel('Solar Time [h]', 'fontsize', 15)
ylabel('T_0 [K]', 'fontsize', 15)
set(gca,'fontsize',15)
xlim([0, 24])

dt = 1;
meanData = [];
stdData = [];
plotT = dt/2 : dt : 12-dt/2;
timeVec = datevec(lbT0Struct.timestamps);
month = timeVec(:,2);
for t = 0:dt:12-dt
    ind = t <= month & month < t+dt;
    meanData = [meanData; mean(lbT0Struct.data(ind))];
    stdData = [stdData; std(lbT0Struct.data(ind))];
end

figure;
errorbar(plotT, meanData, stdData, 's');
xlabel('Month', 'fontsize', 15)
ylabel('T_0 [K]', 'fontsize', 15)
set(gca,'fontsize',15)
xlim([0, 12])

end