function dataAssimilationTimeLoop(modelString, assimilationWindowLength, ensembleSize)
% INPUT:
%     modelString: 'dummy' for test.
%     assimilationWindowLength: in hours.


if strcmpi(modelString, 'dummy')
    loopDummy('2003-10-27', '2003-11-10', 'CHAMP', 'CHAMP', assimilationWindowLength, ensembleSize);
end

end

function loopDummy(beginDay, endDay, assSatellite, plotSatellite, windowLen, ensembleSize)

load ilData rhoStruct

% % Remove unnecessary observations.
t0 = datenum(beginDay);
t1 = datenum(endDay);
% if strcmpi(assSatellite, 'CHAMP')
%     removeInd = ~ismember(1:length(rhoStruct.data), rhoStruct.champ);
%     removeTimes = rhoStruct.timestamps < t0 - windowLen/24/2 | ...
%                     rhoStruct.timestamps > t1;
%     removeInd(removeTimes) = true;
%     dataStruct = removeDataPoints(rhoStruct, removeInd, false,true,true,true);
% end
% 
% if strcmpi(plotSatellite, 'CHAMP')
%     removeInd = ~ismember(1:length(rhoStruct.data), rhoStruct.champ);
%     removeTimes = rhoStruct.timestamps < t0 | ...
%                     rhoStruct.timestamps > t1;
%     removeInd(removeTimes) = true;
%     plotStruct = removeDataPoints(rhoStruct, removeInd, false,true,true,true);
% end

dataStruct = createSyntheticStorm(t0, t1, 120, 'sine');
plotStruct = dataStruct;

N = length(plotStruct.data);
numFields = 3;
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

dataStruct.sigma = 0.05*dataStruct.data;
d = zeros(N,1);
r = zeros(N,1);
assTimes = [];
covsum = [];

assBegin = t0 - windowLen/24/2;
assEnd = t0 + windowLen/24/2;
while assBegin < t1
    removeInd = dataStruct.timestamps < assBegin | dataStruct.timestamps >= assEnd;
    S = removeDataPoints(dataStruct, removeInd);
    S.sigma(removeInd) = [];
    [ensemble,d(~removeInd),r(~removeInd),c] = ...
        assimilateDataAndUpdateEnsemble(ensemble, @dummyThermosphere, S, true);
    covsum = [covsum; sum(diag(c))];
    
    assTime = mean([assBegin; assEnd]);
    assTimes = [assTimes; assTime];
    windowTimes = [assTime, assTime + windowLen/24];
    plotStruct = updatePlotStruct(ensemble, windowTimes, plotStruct, @dummyThermosphere);
    
    assBegin = assBegin  + windowLen/24;
    assEnd = assEnd + windowLen/24;
end

% ensRho = [];
% for i = 1:numFields
%     ensRho = [ensRho, computeOrbitAverage(plotStruct.rho(:,i), ...
%         plotStruct.latitude, plotStruct.timestamps)];
% end
% 
% [dataRho, plotTime] = computeOrbitAverage(plotStruct.data, plotStruct.latitude, plotStruct.timestamps);

figure;
plot(repmat(plotStruct.timestamps,1,numFields+1), [plotStruct.data, plotStruct.rho]);
datetick('x')
set(gca,'fontsize', 15)
ylabel('rho', 'fontsize',15)
ylim([0.8*min(plotStruct.data), 1.25*max(plotStruct.data)])
legend(plotSatellite,'Ens. mean', 'Upper 95%', 'Lower 95%')

figure;
subplot(2,1,1)
plot(repmat(plotStruct.timestamps,1,2), [d, r])
legend('Innovations', 'Residuals');
datetick('x')
set(gca,'fontsize', 15)

subplot(2,1,2)
plot(assTimes, covsum);
legend('sum(diag(B))');
datetick('x')
set(gca,'fontsize', 15)


end

function dataStruct = createSyntheticStorm(t0, t1, dt, shape)
% t0, t1 in days, dt in seconds

timestamps = (t0:dt/86400:t1)';
z = zeros(length(timestamps),1);
alt = 400*ones(length(timestamps),1);
dataStruct = struct('data', z, 'timestamps', timestamps, 'longitude', z, 'latitude', z, 'altitude', alt, 'solarTime', z, 'aeInt', zeros(length(timestamps),7),...
                    'F', z, 'FA', z, 'apNow', z, 'ap3h', z, 'ap6h', z,'ap9h', z, 'ap12To33h', z, 'ap36To57h', z, 'Ap', z);

rho = 6.7E-12*ones(size(timestamps));
N = length(rho);
m = floor(N/2);
if strcmpi(shape,'triangle')
    ascEnd = floor(6*3600/dt);
    descEnd = floor(48*3600/dt);
    base = rho(1);
    maxVal = 2*rho(1);
    ascGrad = (maxVal-base)./(ascEnd);
    descGrad = -(maxVal-base)./(descEnd);
    
    q = floor(m/2); 
    for i = q:q+descEnd
        if i >= N; break; end
        rho(i) = maxVal + (i-q)*descGrad;
    end
    for i = m:m+ascEnd
        if i >= N; break; end
        rho(i) = base + (i-m)*ascGrad;
    end
    l = m + ascEnd+1;
    for i = l:l+descEnd
        if i >= N; break; end
        rho(i) = maxVal + (i-l)*descGrad;
    end
elseif strcmpi(shape,'sine')
    for i = 1:N
        rho(i) = rho(1) + 3E-12*sin(2*pi*i/(N/2));
    end
end

rho = rho + rho.*0.05.*randn(size(rho));

dataStruct.data = rho;
               
end

function plotStruct = updatePlotStruct(ensemble, windowTimes, plotStruct, obsOper)

removeInd = plotStruct.timestamps < windowTimes(1) | plotStruct.timestamps >= windowTimes(2);
S = removeDataPoints(plotStruct, removeInd);
% 
% bestModel = mean(ensemble, 2);
% plotStruct.rho(~removeInd,1) = obsOper(bestModel, S);

ensPredictions = zeros(length(S.data), size(ensemble,2));

for i = 1:size(ensemble,2)
    ensPredictions(:,i) = obsOper(ensemble(:,i), S);
end

plotStruct.rho(~removeInd,1) = mean(ensPredictions,2);
plotStruct.rho(~removeInd,2) = quantile(ensPredictions,0.95,2);
plotStruct.rho(~removeInd,3) = quantile(ensPredictions,0.05,2);
if any(plotStruct.rho < 0)
    a = 1;
end


end