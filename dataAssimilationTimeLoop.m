function dataAssimilationTimeLoop(modelString, assimilationWindowLength, ensembleSize)
% INPUT:
%     modelString: 'dummy' for test, 'full' for complete model
%     assimilationWindowLength: in hours.

Fstd = sqrt(6);
%Fstd = 8;
loopModelAssimilation('2002-08-01', '2003-07-31', 'CH/GR', 'CHAMP', modelString, assimilationWindowLength, ensembleSize, Fstd);


end

function loopModelAssimilation(beginDay, endDay, assSatellite, plotSatellite, modelString, windowLen, ensembleSize, Fstd)

load('ilData.mat', 'rhoStruct', 'OStruct', 'HeStruct', 'N2Struct', 'ArStruct', 'O2Struct')

% % Remove unnecessary observations.
t0 = datenum(beginDay);
t1 = datenum(endDay);
if strcmpi(assSatellite, 'CHAMP')
    removeInd = ~ismember(1:length(rhoStruct.data), rhoStruct.champ);
elseif strcmpi(assSatellite, 'GRACE')
    removeInd = ~ismember(1:length(rhoStruct.data), rhoStruct.grace);
elseif strcmpi(assSatellite, 'CH/GR')
    removeInd = ~ismember(1:length(rhoStruct.data), [rhoStruct.grace, rhoStruct.champ]);
end

removeTimes = rhoStruct.timestamps < t0 - windowLen/24/2 | ...
                rhoStruct.timestamps > t1;
removeInd(removeTimes) = true;
assimiStruct = removeDataPoints(rhoStruct, removeInd, false,true,true,true);

if strcmpi(plotSatellite, 'CHAMP')
    removeInd = ~ismember(1:length(rhoStruct.data), rhoStruct.champ);
    removeTimes = rhoStruct.timestamps < t0 | ...
                    rhoStruct.timestamps > t1;
    removeInd(removeTimes) = true;
    analStruct = removeDataPoints(rhoStruct, removeInd, false,true,true,true);
    analStruct = computeVariablesForFit(analStruct);
elseif strcmpi(plotSatellite, 'GRACE')
    removeInd = ~ismember(1:length(rhoStruct.data), rhoStruct.grace);
    removeTimes = rhoStruct.timestamps < t0 | ...
                    rhoStruct.timestamps > t1;
    removeInd(removeTimes) = true;
    analStruct = removeDataPoints(rhoStruct, removeInd, false,true,true,true);
    analStruct = computeVariablesForFit(analStruct);
end

%dataStruct = createSyntheticStorm(t0, t1, 120, 'triangle');
%plotStruct = dataStruct;

N = length(analStruct.data);
numFields = 3;
analStruct.rho = zeros(N,numFields);
analStruct.O = zeros(N,numFields);
analStruct.He = zeros(N,numFields);
analStruct.N2 = zeros(N,numFields);
analStruct.Ar = zeros(N,numFields);
analStruct.O2 = zeros(N,numFields);
analStruct.Tex = zeros(N,numFields);
analStruct.T0 = zeros(N,numFields);
analStruct.dT = zeros(N,numFields);
analStruct.T = zeros(N,numFields);

ensemble = createInitialEnsemble(modelString, ensembleSize);

if strcmpi(modelString,'dummy')
    refModel = dummyThermosphere(zeros(size(ensemble,1),1), analStruct, 1);
    modelOperator = @dummyThermosphere;
else
    assimiStruct = addCoeffsToStruct(assimiStruct, OStruct, HeStruct, N2Struct, ArStruct, O2Struct);
    analStruct = addCoeffsToStruct(analStruct, OStruct, HeStruct, N2Struct, ArStruct, O2Struct);
    
    if ~exist('rhoStruct_il.mat', 'file')
        load optCoeff.mat
        TexStruct.coeffInd = TexInd;
        coeffStruct = struct('TexCoeff' , optCoeff(TexInd),... 
        'OCoeff', optCoeff(OInd),...
        'N2Coeff' , optCoeff(N2Ind),...
        'HeCoeff' , optCoeff(HeInd),...
        'O2Coeff' , optCoeff(O2Ind),...
        'ArCoeff' , optCoeff(ArInd),...
        'dTCoeff', dTCoeffs,...
        'T0Coeff', T0Coeffs);
        
        numBiasesStruct = struct('O', 5, 'N2', 6,...
        'He', 5, 'Ar', 2, 'O2', 0); 
        
        [ilRho] = computeComparisonData(rhoStruct, coeffStruct, numBiasesStruct);
        
        save('rhoStruct_il.mat','ilRho')
    end
    
    load rhoStruct_il
    
    refModel = ilRho(~removeInd);
    %refModel = 1.05*mean(plotStruct.data./refModel) * refModel; % TESTAUS
    modelOperator = @il_model_operator;
end
bgStruct = analStruct;

%assimiStruct.sigma = 0.05*assimiStruct.data;
a = assimiStruct.sigma ./ assimiStruct.data;
d = zeros(N,1);
r = zeros(N,1);
assTimes = [];
ensMeans = zeros(size(ensemble,1),1);
ensMeanF = [];
covdiag = zeros(0,size(ensemble,1));

assBegin = t0 - windowLen/24/2;
assEnd = t0 + windowLen/24/2;
prevRmInd = ones(size(assimiStruct.timestamps));
step = 1;
obsRank = ones(size(assimiStruct.data));
while assBegin < t1
    datestr(assBegin)
    removeInd = assimiStruct.timestamps < assBegin | assimiStruct.timestamps >= assEnd | ~prevRmInd;
    S = removeDataPoints(assimiStruct, removeInd);
    S.sigma = S.sigma ./ S.data;
    S.data = log(S.data);
    if isempty(S.data)
        continue
    end
    S = computeVariablesForFit(S);
    
    previousEnsemble = ensemble;
    
    %ensMean = mean(ensemble, 2);
    %ensStd = std(ensemble, 0, 2);
    %ensemble = bsxfun(@plus, (1 + inflFac)*bsxfun(@minus, ensemble, ensMean), ensMean);
%     if step > 1
%         ensemble(1,:) = (Fstd)/ensStd(1) * (ensemble(1,:)-ensMean(1)) + ensMean(1);
%         %ensemble(2,:) = (0.2)/ensStd(2) * (ensemble(2,:)-ensMean(2)) + ensMean(2);
%     end
    [ensemble,d(~removeInd),r(~removeInd),c,P,obsRank(~removeInd)] = ...
        assimilateDataAndUpdateEnsemble(ensemble, modelOperator, S, true);
    ensMean = mean(ensemble, 2);
    ensStd = std(ensemble, 0, 2);
   % if step > 1
        ensemble(1,:) = (Fstd)/ensStd(1) * (ensemble(1,:)-ensMean(1)) + ensMean(1);
        %ensemble(2,:) = (0.2)/ensStd(2) * (ensemble(2,:)-ensMean(2)) + ensMean(2);
    %end
    covdiag = [covdiag; diag(c)'];
    ensMeans = [ensMeans, ensMean];
    ensMeanF = [ensMeanF; ensMean(1)];
    
    assTime = mean([assBegin; assEnd]);
    assTimes = [assTimes; assTime];
    windowTimes = [assTime, assTime + windowLen/24];
    analStruct = updatePlotStruct(ensemble, windowTimes, analStruct, modelOperator);
    bgStruct = updatePlotStruct(previousEnsemble, windowTimes, bgStruct, modelOperator);
    
    assBegin = assBegin  + windowLen/24;
    assEnd = assEnd + windowLen/24;
    prevRmInd = removeInd;
    step = step + 1;
end
ensMeans = ensMeans(:,2:end);

d = d(1:N);
r = r(1:N);

analRho = [];
bgRho = [];
for i = 1:numFields
    analRho = [analRho, computeOrbitAverage(analStruct.rho(:,i), ...
        analStruct.latitude, analStruct.timestamps)];
    
    bgRho = [bgRho, computeOrbitAverage(bgStruct.rho(:,i), ...
        bgStruct.latitude, bgStruct.timestamps)];
end

[dataRho, plotTime] = computeOrbitAverage(analStruct.data, analStruct.latitude, analStruct.timestamps);
[refOrbAver, plotTime] = computeOrbitAverage(refModel, analStruct.latitude, analStruct.timestamps);
t = plotTime;
i = t(1)+5 <= t & t <= t(end); 
%ensRho(:,1) = 1.05*mean(dataRho(i)./ensRho(i,1)) * ensRho(:,1); % TESTAUS

figure;
% subplot(2,1,1)
% if strcmpi(modelString,'dummy')
    %plot(repmat(plotTime,1,numFields+1), [dataRho, ensRho], 'linewidth', 2.0);
    %legend(plotSatellite,'Ens. mean', '+ 1 \sigma', '- 1 \sigma')
% else
     plot(repmat(plotTime,1,numFields+2), [dataRho, refOrbAver, analRho], 'linewidth', 2.0);
     legend(plotSatellite, 'IL', 'Keskiarvo', '+ 1 \sigma', '- 1 \sigma')
% end
datetick('x')
title('Tiheys','fontsize',15)
set(gca,'fontsize', 15)
ylabel('Tiheys [kg/m^3]', 'fontsize',15)
axis tight
ylim([0.8*min(analStruct.data), 1.25*max(analStruct.data)])

% subplot(2,1,2)
% plot(repmat(plotStruct.timestamps,1,numFields), plotStruct.T, 'linewidth', 2.0);
% datetick('x')
% title('Tex','fontsize',15)
% set(gca,'fontsize', 15)
% ylabel('rho', 'fontsize',15)
% ylim([0.8*min(plotStruct.data), 1.25*max(plotStruct.data)])
% legend(plotSatellite,'Ens. mean', 'Upper 95%', 'Lower 95%')
% axis tight

%[refOrbAver] = computeOrbitAverage(refModel, plotStruct.latitude, plotStruct.timestamps);

refRMS = rms(dataRho(i)./refOrbAver(i)-1);
ensMeanRMS = rms(dataRho(i)./analRho(i,1)-1);

fprintf('Unadjusted RMS: %f\n', refRMS);
fprintf('EnKF FGAT RMS: %f\n', ensMeanRMS);
fprintf('Improvement: %f\n', (1 - (ensMeanRMS./refRMS))*100);

figure;
subplot(3,1,1)
plot(repmat(analStruct.timestamps,1,2), [d, r])
legend('Innovaatiot', 'Residuaalit');
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
legend('\sigma F30', '\sigma T0','\sigma dT');
datetick('x')
set(gca,'fontsize', 15)

figure;
plot(analStruct.timestamps, analStruct.Tex)
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
T = addCoeffsToStruct(T, OStruct, HeStruct, N2Struct, ArStruct, O2Struct);
T = computeVariablesForFit(T);
T.aeInt = mean(S.aeInt);
z = bulge(ensMean(4:end),T);
z = reshape(z, size(x,1), size(x,2));
surf(x,y,z,'edgecolor','none')
shading interp
colorbar
view(2);
xlabel('Paikallisaika','fontsize',15)
ylabel('Leveyspiiri','fontsize',15)
title('Tex -korjaus','fontsize',15)
axis tight

refRMS = [];
bgRMS = [];
analRMS = [];
for t = floor(t0):ceil(t1)-1
    ind = t < plotTime & plotTime <= t+1;
    refRMS = [refRMS; rms(dataRho(ind)./refOrbAver(ind)-1)];
    bgRMS = [bgRMS; rms(dataRho(ind)./bgRho(ind,1)-1)];
    analRMS = [analRMS; rms(dataRho(ind)./analRho(ind,1)-1)];
end
figure;
t = (floor(t0):ceil(t1)-1)';
plot(repmat(t,1,3),[refRMS, bgRMS, analRMS])
set(gca,'fontsize', 15)
xlabel('Aika','fontsize',15)
ylabel('RMSE','fontsize',15)
legend('IL','Tausta','Analyysi')
datetick('x')

figure('renderer', 'zbuffer');
zAnalysis = reshape(modelOperator(ensMean,T), size(x,1), size(x,2));
bgMean = mean(previousEnsemble,2);
zBg = reshape(modelOperator(bgMean,T), size(x,1), size(x,2));

[~,ind] = min(abs(assimiStruct.timestamps-assTimes(end)));
controlEns = [assimiStruct.FA(ind); zeros(length(ensMean)-1,1)];
zIL = reshape(modelOperator(controlEns,T), size(x,1), size(x,2));

% subplot(4,1,1)
% scatter3(S.solarTime, S.latitude, zeros(size(S.data)), 10.0, exp(S.data), ...
%     'MarkerFaceColor', 'flat','Marker','o')
% shading interp; colorbar; view(2);
% %xlabel('lst','fontsize',15)
% ylabel('lat','fontsize',15)
% title('IL','fontsize',15)
% axis tight

subplot(3,1,1)
surf(x,y,exp(zIL),'edgecolor','none')
shading interp; colorbar; view(2);
%xlabel('lst','fontsize',15)
ylabel('Leveyspiiri','fontsize',15)
title('IL','fontsize',15)
axis tight

subplot(3,1,2)
surf(x,y,exp(zBg),'edgecolor','none')
shading interp; colorbar; view(2);
%xlabel('lst','fontsize',15)
ylabel('Leveyspiiri','fontsize',15)
title('Background','fontsize',15)
axis tight

subplot(3,1,3)
surf(x,y,exp(zAnalysis),'edgecolor','none')
shading interp; colorbar; view(2);
xlabel('Paikallisaika','fontsize',15)
ylabel('Leveyspiiri','fontsize',15)
title('Analysis','fontsize',15)
axis tight

figure;
plot(assTimes, ensMeanF, assimiStruct.timestamps, assimiStruct.F)
datetick('x')
set(gca,'fontsize', 15)
xlabel('Aika','fontsize',15)
ylabel('F30','fontsize',15)

figure;
surf(log10(abs(c)),'edgecolor','none');
title('Mallin kovarianssi', 'fontsize', 15)
view(2);
colorbar;
set(gca,'fontsize', 15)
set(gca,'ydir','reverse')
axis tight;

v = ver;
if datenum(v(1).Date) > datenum('2015-01-01')
    figure;
    numBins = 20;
    t = assimiStruct.timestamps;
    i = t(1)+60 <= t & t <= t(end);
    h = histogram(obsRank(i), numBins);
    counts = get(h,'Values');
    title('Sijoituslukujen jakauma', 'fontsize', 15)
    expected = length(assimiStruct.data(i)) / numBins;
    chiSqStat = sum((counts-expected).^2 ./ expected);
    fprintf('Chi Sq.: %f\n', chiSqStat);
end

writeCorrectedModelToFiles(analStruct, plotSatellite);

filename = [plotSatellite,'_',beginDay,'_',endDay,'.mat'];
save(filename, 'bgStruct')
save(filename, 'analStruct','assimiStruct','d','r','P','obsRank','dataRho','analRho','bgRho','refOrbAver',...
    'plotTime', 'covdiag','ensMeanF','ensMeans','assTimes','ensemble','windowLen', '-append')

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

ensPredictions = exp(ensPredictions);
plotStruct.rho(~removeInd,1) = mean(ensPredictions,2);
plotStruct.rho(~removeInd,2) = quantile(ensPredictions,0.025,2);
plotStruct.rho(~removeInd,3) = quantile(ensPredictions,0.975,2);
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
    if isempty(ind)
        date = date + 1;
        continue
    end
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
