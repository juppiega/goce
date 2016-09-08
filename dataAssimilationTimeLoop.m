function dataAssimilationTimeLoop(modelString, assimilationWindowLength, ensembleSize)
% INPUT:
%     modelString: 'dummy' for test.
%     assimilationWindowLength: in hours.


if strcmpi(modelString, 'dummy')
    loopDummy('2003-10-28', '2003-11-03', 'CHAMP', 'CHAMP', assimilationWindowLength, ensembleSize);
end

end

function loopDummy(beginDay, endDay, assSatellite, plotSatellite, windowLen, ensembleSize)

load ilData originalRhoStruct

% Remove unnecessary observations.
t0 = datenum(beginDay);
t1 = datenum(endDay);
if strcmpi(assSatellite, 'CHAMP')
    removeInd = ~ismember(1:length(originalRhoStruct.data), originalRhoStruct.champ);
    removeTimes = originalRhoStruct.timestamps < t0 - windowLen/24/2 | ...
                    originalRhoStruct.timestamps > t1;
    removeInd(removeTimes) = true;
    dataStruct = removeDataPoints(originalRhoStruct, removeInd, false,true,true,true);
end

if strcmpi(plotSatellite, 'CHAMP')
    removeInd = ~ismember(1:length(originalRhoStruct.data), originalRhoStruct.champ);
    removeTimes = originalRhoStruct.timestamps < t0 | ...
                    originalRhoStruct.timestamps > t1;
    removeInd(removeTimes) = true;
    plotStruct = removeDataPoints(originalRhoStruct, removeInd, false,true,true,true);
end

N = length(plotStruct.data);
numFields = 4;
plotStruct.rho = zeros(N,numFields);
plotStruct.O = zeros(N,numFields);
plotStruct.He = zeros(N,numFields);
plotStruct.N2 = zeros(N,numFields);
plotStruct.Ar = zeros(N,numFields);
plotStruct.O2 = zeros(N,numFields);
plotStruct.Tex = zeros(N,numFields);
plotStruct.T0 = zeros(N,numFields);
plotStruct.dT = zeros(N,numFields);
plotStruct.T = zeros(N,numFields);

ensemble = createInitialEnsemble('dummy', ensembleSize);

dataStruct.sigma = 0.1*dataStruct.data;

assBegin = t0 - windowLen/24/2;
assEnd = t0 + windowLen/24/2;
while assBegin < t1
    removeInd = dataStruct.timestamps < assBegin | dataStruct.timestamps >= assEnd;
    S = removeDataPoints(dataStruct, removeInd, false,true,true,true);
    S.sigma(removeInd) = [];
    ensemble = assimilateDataAndUpdateEnsemble(ensemble, @dummyThermosphere, S);
    
    assTime = mean([assBegin; assEnd]);
    windowTimes = [assTime, assTime + windowLen/24];
    plotStruct = updatePlotStruct(ensemble, windowTimes, plotStruct, @dummyThermosphere);
    
    assBegin = assBegin  + windowLen/24;
    assEnd = assEnd + windowLen/24;
end

figure;
plot(plotStruct.timestamps, plotStruct.data, plotStruct.timestamps, plotStruct.rho(:,1));
datetick('x')
set(gca,'fontsize', 15)
ylabel('rho', 'fontsize',15)

end

function plotStruct = updatePlotStruct(ensemble, windowTimes, plotStruct, obsOper)

removeInd = plotStruct.timestamps < windowTimes(1) | plotStruct.timestamps >= windowTimes(2);
S = removeDataPoints(plotStruct, removeInd, false,true,true,true);

bestModel = mean(ensemble, 2);
plotStruct.rho(~removeInd,1) = obsOper(bestModel, S);


end