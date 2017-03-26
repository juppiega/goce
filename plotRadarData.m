function [] = plotRadarData(radarNum)

if exist('ilData.mat', 'file')
    load ilData.mat
else
    error('File ilData.mat not found!')
end


[rhoStruct, TempStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct] = ...
    removeAndFixData(rhoStruct, 0, TempStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct, true);

removeInd = lbT0Struct.index ~= radarNum;
%removeInd(lbT0Struct.FA > 120) = true;
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
title(num2str(radarNum),'fontsize',15);

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
title(num2str(radarNum),'fontsize',15);

dFA = 20;
meanData = [];
stdData = [];
FA = lbT0Struct.FA;
plotFA = 60 : dFA : 220;
for f = min(plotFA)-dFA/2 : dFA : max(plotFA)-dFA/2
    ind = f <= FA & FA < f+dFA;
    meanData = [meanData; mean(lbT0Struct.data(ind))];
    stdData = [stdData; std(lbT0Struct.data(ind))];
end

figure;
subplot(1,2,1)
hist(lbT0Struct.FA, 20)
set(gca,'fontsize',15)
title(num2str(radarNum),'fontsize',15);

subplot(1,2,2)
errorbar(plotFA, meanData, stdData, 's');
xlabel('FA', 'fontsize', 15)
ylabel('T_0 [K]', 'fontsize', 15)
set(gca,'fontsize',15)
title(num2str(radarNum),'fontsize',15);
axis tight

end