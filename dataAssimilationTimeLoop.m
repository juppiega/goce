function dataAssimilationTimeLoop(modelString, assimilationWindowLength, ensembleSize, TexStd)
% INPUT:
%     modelString: 'dummy' for test, 'full' for complete model
%     assimilationWindowLength: in hours.


loopModelAssimilation('2002-08-10', '2002-08-24', 'CHAMP', 'GRACE', modelString, assimilationWindowLength, ensembleSize, TexStd);


end

function loopModelAssimilation(beginDay, endDay, assSatellite, plotSatellite, modelString, windowLen, ensembleSize, TexStd)

load('ilData.mat', 'rhoStruct', 'OStruct', 'HeStruct', 'N2Struct', 'ArStruct', 'O2Struct')

% % Remove unnecessary observations.
t0 = datenum(beginDay);
t1 = datenum(endDay);
if strcmpi(assSatellite, 'CHAMP')
    removeInd = ~ismember(1:length(rhoStruct.data), rhoStruct.champ);
    removeTimes = rhoStruct.timestamps < t0 - windowLen/24/2 | ...
                    rhoStruct.timestamps > t1;
    removeInd(removeTimes) = true;
    assimiStruct = removeDataPoints(rhoStruct, removeInd, false,true,true,true);
elseif strcmpi(assSatellite, 'GRACE')
    removeInd = ~ismember(1:length(rhoStruct.data), rhoStruct.grace);
    removeTimes = rhoStruct.timestamps < t0 - windowLen/24/2 | ...
                    rhoStruct.timestamps > t1;
    removeInd(removeTimes) = true;
    assimiStruct = removeDataPoints(rhoStruct, removeInd, false,true,true,true);
end

if strcmpi(plotSatellite, 'CHAMP')
    removeInd = ~ismember(1:length(rhoStruct.data), rhoStruct.champ);
    removeTimes = rhoStruct.timestamps < t0 | ...
                    rhoStruct.timestamps > t1;
    removeInd(removeTimes) = true;
    plotStruct = removeDataPoints(rhoStruct, removeInd, false,true,true,true);
    plotStruct = computeVariablesForFit(plotStruct);
elseif strcmpi(plotSatellite, 'GRACE')
    removeInd = ~ismember(1:length(rhoStruct.data), rhoStruct.grace);
    removeTimes = rhoStruct.timestamps < t0 | ...
                    rhoStruct.timestamps > t1;
    removeInd(removeTimes) = true;
    plotStruct = removeDataPoints(rhoStruct, removeInd, false,true,true,true);
    plotStruct = computeVariablesForFit(plotStruct);
end

%dataStruct = createSyntheticStorm(t0, t1, 120, 'triangle');
%plotStruct = dataStruct;

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

ensemble = createInitialEnsemble(modelString, ensembleSize);

if strcmpi(modelString,'dummy')
    refModel = dummyThermosphere(zeros(size(ensemble,1),1), plotStruct, 1);
    modelOperator = @dummyThermosphere;
else
    assimiStruct = addCoeffsToStruct(assimiStruct, OStruct, HeStruct, N2Struct, ArStruct, O2Struct);
    plotStruct = addCoeffsToStruct(plotStruct, OStruct, HeStruct, N2Struct, ArStruct, O2Struct);
    refModel = il_model_operator(zeros(size(ensemble,1),1), plotStruct, 1);
    refModel = 1.05*mean(plotStruct.data./refModel) * refModel; % TESTAUS
    modelOperator = @il_model_operator;
end

%assimiStruct.sigma = 0.05*assimiStruct.data;
a = assimiStruct.sigma ./ assimiStruct.data;
d = zeros(N,1);
r = zeros(N,1);
assTimes = [];
covdiag = zeros(0,size(ensemble,1));

assBegin = t0 - windowLen/24/2;
assEnd = t0 + windowLen/24/2;
prevRmInd = ones(size(assimiStruct.timestamps));
step = 1;
while assBegin < t1
    removeInd = assimiStruct.timestamps < assBegin | assimiStruct.timestamps >= assEnd | ~prevRmInd;
    S = removeDataPoints(assimiStruct, removeInd);
    if isempty(S.data)
        continue
    end
    S = computeVariablesForFit(S);
    
    ensMean = mean(ensemble, 2);
    ensStd = std(ensemble, 0, 2);
    %ensemble = bsxfun(@plus, (1 + inflFac)*bsxfun(@minus, ensemble, ensMean), ensMean);
    if step > 1
        ensemble(3,:) = (TexStd)/ensStd(3) * (ensemble(3,:)-ensMean(3)) + ensMean(3);
        %ensemble(2,:) = (0.2)/ensStd(2) * (ensemble(2,:)-ensMean(2)) + ensMean(2);
    end
    [ensemble,d(~removeInd),r(~removeInd),c,P] = ...
        assimilateDataAndUpdateEnsemble(ensemble, modelOperator, S, true);
    covdiag = [covdiag; diag(c)'];
    
    assTime = mean([assBegin; assEnd]);
    assTimes = [assTimes; assTime];
    windowTimes = [assTime, assTime + windowLen/24];
    plotStruct = updatePlotStruct(ensemble, windowTimes, plotStruct, modelOperator);
    
    assBegin = assBegin  + windowLen/24;
    assEnd = assEnd + windowLen/24;
    prevRmInd = removeInd;
    step = step + 1;
end

d = d(1:N);
r = r(1:N);

ensRho = [];
for i = 1:numFields
    ensRho = [ensRho, computeOrbitAverage(plotStruct.rho(:,i), ...
        plotStruct.latitude, plotStruct.timestamps)];
end

[dataRho, plotTime] = computeOrbitAverage(plotStruct.data, plotStruct.latitude, plotStruct.timestamps);
[refOrbAver, plotTime] = computeOrbitAverage(refModel, plotStruct.latitude, plotStruct.timestamps);
t = plotTime;
i = t(end)-10 <= t & t <= t(end); % Last 10 days
ensRho(:,1) = 1.05*mean(dataRho(i)./ensRho(i,1)) * ensRho(:,1); % TESTAUS

figure;
% subplot(2,1,1)
% if strcmpi(modelString,'dummy')
    plot(repmat(plotTime,1,numFields+1), [dataRho, ensRho], 'linewidth', 2.0);
    legend(plotSatellite,'Ens. mean', '+ 1 \sigma', '- 1 \sigma')
% else
%     plot(repmat(plotTime,1,numFields+2), [dataRho, refOrbAver, ensRho], 'linewidth', 2.0);
%     legend(plotSatellite, 'IL', 'Ens. mean', '+ 1 \sigma', '- 1 \sigma')
% end
datetick('x')
title('Rho','fontsize',15)
set(gca,'fontsize', 15)
ylabel('rho', 'fontsize',15)
axis tight
ylim([0.8*min(plotStruct.data), 1.25*max(plotStruct.data)])

% subplot(2,1,2)
% plot(repmat(plotStruct.timestamps,1,numFields), plotStruct.T, 'linewidth', 2.0);
% datetick('x')
% title('Tex','fontsize',15)
% set(gca,'fontsize', 15)
% ylabel('rho', 'fontsize',15)
% ylim([0.8*min(plotStruct.data), 1.25*max(plotStruct.data)])
% legend(plotSatellite,'Ens. mean', 'Upper 95%', 'Lower 95%')
% axis tight

[refOrbAver] = computeOrbitAverage(refModel, plotStruct.latitude, plotStruct.timestamps);

refRMS = rms(dataRho(i)./refOrbAver(i)-1);
ensMeanRMS = rms(dataRho(i)./ensRho(i,1)-1);

fprintf('Unadjusted RMS: %f\n', refRMS);
fprintf('EnKF FGAT RMS: %f\n', ensMeanRMS);
fprintf('Improvement: %f\n', (1 - (ensMeanRMS./refRMS))*100);

figure;
subplot(3,1,1)
plot(repmat(plotStruct.timestamps,1,2), [d, r])
legend('Innovations', 'Residuals');
datetick('x')
set(gca,'fontsize', 15)

subplot(3,1,2)
plot(assTimes, sum(covdiag,2));
legend('sum(diag(B))');
datetick('x')
set(gca,'fontsize', 15)

subplot(3,1,3)
tVar = sqrt(covdiag(:,1:3));
plot(repmat(assTimes,1,3), tVar);
legend('\sigma T0', '\sigma dT','\sigma Tex');
datetick('x')
set(gca,'fontsize', 15)

figure;
plot(plotStruct.timestamps, plotStruct.Tex)
title('Tex', 'fontsize', 15)
datetick('x')
set(gca,'fontsize', 15)

figure('renderer', 'zbuffer');
[x,y] = meshgrid(0:0.5:24, -90:5:90);
T.latitude = y(:);
T.solarTime = x(:);
T.altitude = 400*ones(size(x(:)));
T.doy = 180*ones(size(x(:)));
T.timestamps = now*ones(size(x(:)));
T.longitude = 180*ones(size(x(:)));
T = computeVariablesForFit(T);
z = bulge(ensMean(3:end),T);
z = reshape(z, size(x,1), size(x,2));
surf(x,y,z,'edgecolor','none')
shading interp
colorbar
view(2);
xlabel('lst','fontsize',15)
ylabel('lat','fontsize',15)
title('Tex correction','fontsize',15)
axis tight


figure;
surf(log10(abs(c)),'edgecolor','none');
title('Model covariance', 'fontsize', 15)
view(2);
colorbar;
set(gca,'fontsize', 15)
set(gca,'ydir','reverse')
axis tight;


writeCorrectedModelToFiles(plotStruct, plotSatellite);

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
if isempty(S.data)
    return
end
S = computeVariablesForFit(S);
% 
% bestModel = mean(ensemble, 2);
% plotStruct.rho(~removeInd,1) = obsOper(bestModel, S);

ensPredictions = zeros(length(S.data), size(ensemble,2));
s = struct('O', [], 'N2', [], 'He', [], 'Ar', [], 'O2', [], 'Tex', [],...
                            'dT', [], 'T', [], 'T0', []);
ensVals = repmat(s,size(ensemble,2),1);

for i = 1:size(ensemble,2)
    [ensPredictions(:,i), ensVals(i)] = obsOper(ensemble(:,i), S);
end

Tarray = [ensVals.T];
plotStruct.T(~removeInd,1) = mean(Tarray,2);
plotStruct.T(~removeInd,2) = mean(Tarray,2) - std(Tarray,[],2);
plotStruct.T(~removeInd,3) = mean(Tarray,2) + std(Tarray,[],2);

TexArray = [ensVals.Tex];
plotStruct.Tex(~removeInd,1) = mean(TexArray,2);
plotStruct.Tex(~removeInd,2) = mean(TexArray,2) - std(TexArray,[],2);
plotStruct.Tex(~removeInd,3) = mean(TexArray,2) + std(TexArray,[],2);

plotStruct.rho(~removeInd,1) = mean(ensPredictions,2);
plotStruct.rho(~removeInd,2) = mean(ensPredictions,2) - std(ensPredictions,[],2);
plotStruct.rho(~removeInd,3) = mean(ensPredictions,2) + std(ensPredictions,[],2);
if any(plotStruct.rho < 0)
    a = 1;
end


end

function writeCorrectedModelToFiles(plotStruct, plotSatellite)

date = floor(plotStruct.timestamps(1));
[year,~,~,~,~,~] = datevec(date);
doy = round(plotStruct.doy(1));
foldername = ['tim/',plotSatellite, '_il_', sprintf('%.2d',year-2000),'_',sprintf('%.3d',doy)];
if exist(foldername,'dir')
    rmdir(foldername,'s');
end
mkdir(foldername);

while plotStruct.timestamps(end) > date
    ind = find(date <= plotStruct.timestamps & plotStruct.timestamps < date + 1);
    doy = round(plotStruct.doy(ind(1)));
    sec = (plotStruct.timestamps(ind) - date) * 86400;
    
    sigma_rho = plotStruct.rho(ind,3) - plotStruct.rho(ind,1);
    
    Dil = struct('data', plotStruct.rho(ind,1));
    Doy = struct('data', doy);
    Height = struct('data', plotStruct.altitude(ind));
    Lat = struct('data', plotStruct.latitude(ind));
    Lon = struct('data', plotStruct.longitude(ind));
    LocTim = struct('data', plotStruct.solarTime(ind));
    Sec = struct('data', sec);
    U_rho = struct('data', sigma_rho);
    Year = struct('data', year - 2000);
    
    filename = [foldername,'/',plotSatellite,'_il_',sprintf('%.2d',year-2000),...
        '_',sprintf('%.3d',doy),'.mat'];
    save(filename, 'Dil', 'Doy', 'Height', 'Lat', 'Lon', 'LocTim', 'Sec', 'U_rho', 'Year');
    
    date = date + 1;
    
end

stat = system(['tar -zcf ',foldername,'.tar.gz ',foldername]);

end