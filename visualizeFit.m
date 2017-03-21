function [] = visualizeFit(saveFolder, includeSwarm, fullscreenFigs)
aeThreshold = 0;

if nargin == 2
    fullscreenFigs = true;
end

if ~exist(saveFolder, 'file')
    mkdir(saveFolder)
end

if exist('optCoeff.mat', 'file')
    load optCoeff.mat
else
    error('File optCoeff.mat not found!')
end

% Check the existence of the data file.
if exist('ilData.mat', 'file')
    load('ilData.mat', 'originalRhoStruct', 'dTCoeffs', 'T0Coeffs')
else
    error('File ilData.mat not found!')
end

[originalRhoStruct] = removeAndFixData(originalRhoStruct, aeThreshold);

if ~includeSwarm
    removeInd = originalRhoStruct.swarm;
    originalRhoStruct = removeDataPoints(originalRhoStruct, removeInd, false, true, false, false);
end

 numBiasesStruct = struct('O', 5, 'N2', 6,...
     'He', 5, 'Ar', 2, 'O2', 0); % TODO: paivita arvot lisattyasi datasetteja

%numBiasesStruct = struct('O', 0, 'N2', 0,...
%    'He', 0, 'Ar', 0, 'O2', 0); % TODO: paivita arvot lisattyasi datasetteja

TexStruct.coeffInd = TexInd;
coeffStruct = struct('TexCoeff' , optCoeff(TexInd),... 
'OCoeff', optCoeff(OInd),...
'N2Coeff' , optCoeff(N2Ind),...
'HeCoeff' , optCoeff(HeInd),...
'O2Coeff' , optCoeff(O2Ind),...
'ArCoeff' , optCoeff(ArInd),...
'dTCoeff', dTCoeffs,...
'T0Coeff', T0Coeffs);

z = 400;
lat = -90:5:90;
lst = 0:0.5:24;
lon = 0;
doy = 1;
F = 200;
FA = 200;
aeInt = 20*ones(1,7);
Ap = 3;
lstMean = false;
lonMean = false;
latitudeMean = false;
devFromXmean = false;
sameColorBars = false;
plotSurfs(z, lat, lst, lon, doy, F, FA, aeInt, Ap, lstMean, lonMean, latitudeMean, devFromXmean, ...
    sameColorBars, 'yx', 'dT', coeffStruct, numBiasesStruct);

if exist('msisDtmComparison.mat', 'file')
    load msisDtmComparison.mat
else
    [~, msisRho, dtmRho] = computeComparisonData(originalRhoStruct, coeffStruct, numBiasesStruct);

    save('msisDtmComparison.mat', 'msisRho')
    save('msisDtmComparison.mat', 'dtmRho', '-append')
end

if exist('ilComparison.mat', 'file')
    load ilComparison.mat
else
    [ilRho] = computeComparisonData(originalRhoStruct, coeffStruct, numBiasesStruct);

    save('ilComparison.mat', 'ilRho')
end

if ~isfield(originalRhoStruct, 'dst')
    originalRhoStruct = computeDst(originalRhoStruct);
    save('ilData.mat', 'originalRhoStruct', '-append');
end

modelStruct = struct('il', ilRho, 'msis', msisRho, 'dtm', dtmRho);

%   plot3DOM(originalRhoStruct.dst, 25, originalRhoStruct.latitude, 10, originalRhoStruct.data,...
%    modelStruct, 'O/M', 'Kp3h', 'lat', saveFolder,fullscreenFigs);
% %  plot3DOM(originalRhoStruct.aeInt(:,4), 50, originalRhoStruct.solarTime, 2, originalRhoStruct.data,...
% %   modelStruct, 'O/M', 'AE16h', 'lst', saveFolder,fullscreenFigs);
% % % % plot3DOM(originalRhoStruct.aeInt(:,4), 50, originalRhoStruct.altitude, 25, originalRhoStruct.data,...
% % % %  modelStruct, 'O/M', 'AE16h', 'alt', saveFolder);
% % % % plot3DOM(originalRhoStruct.aeInt(:,4), 50, originalRhoStruct.doy, 30, originalRhoStruct.data,...
% % % %  modelStruct, 'O/M', 'AE16h', 'doy', saveFolder);
% % % % plot3DOM(originalRhoStruct.aeInt(:,4), 50, originalRhoStruct.FA, 10, originalRhoStruct.data,...
% % % %  modelStruct, 'O/M', 'AE16h', 'FA', saveFolder);
% % % % plot3DOM(originalRhoStruct.aeInt(:,4), 50, originalRhoStruct.F - originalRhoStruct.FA, ...
% % % %     10, originalRhoStruct.data,...
% % % %  modelStruct, 'O/M', 'AE16h', 'F-FA', saveFolder);
% % % % 
% % % % plot2DOM(originalRhoStruct.aeInt(:,4), 50, originalRhoStruct.data, modelStruct, 'O/M', 'AE16h', saveFolder)
% % % 
 computeStatistics(originalRhoStruct, ilRho, msisRho, dtmRho, saveFolder);
% % % % 
%     plotStormFig(originalRhoStruct, modelStruct, '2003-10-27', '2003-11-02', 'CHAMP', coeffStruct, numBiasesStruct, saveFolder,fullscreenFigs);
%     plotStormFig(originalRhoStruct, modelStruct, '2010-04-03', '2010-04-08', 'GOCE', coeffStruct, numBiasesStruct, saveFolder,fullscreenFigs);
%     plotStormFig(originalRhoStruct, modelStruct, '2007-03-22', '2007-03-26', 'GRACE', coeffStruct, numBiasesStruct, saveFolder,fullscreenFigs);
%     plotStormFig(originalRhoStruct, modelStruct, '2006-12-13', '2006-12-17', 'GRACE', coeffStruct, numBiasesStruct, saveFolder,fullscreenFigs);
%     plotStormFig(originalRhoStruct, modelStruct, '2011-05-26', '2011-05-31', 'GOCE', coeffStruct, numBiasesStruct, saveFolder,fullscreenFigs);
%     plotStormFig(originalRhoStruct, modelStruct, '2013-06-26', '2013-07-03', 'GOCE', coeffStruct, numBiasesStruct, saveFolder,fullscreenFigs);
%plotStormFig(originalRhoStruct, modelStruct, '2015-04-09', '2015-04-14', 'SWARM', coeffStruct, numBiasesStruct, saveFolder,fullscreenFigs);
%
% analyzeStormTimes(originalRhoStruct, modelStruct, saveFolder,fullscreenFigs);

end

function [Tex, dT0, T0] = findTempsForFit(varStruct, TexCoeffs, dTCoeffs, T0Coeffs)

Tex_est = evalTex(varStruct, TexCoeffs);

T0 = clamp(200, evalT0(varStruct, T0Coeffs), 1000);
dT0 = clamp(1, evalDT(varStruct, dTCoeffs), 30);
Tex = clamp(T0+1, Tex_est, 5000);

end

function [] = plotSurfs(altitude, lat, lst, lon, doy, F, FA, aeInt, Ap, lstMean, lonMean, latitudeMean, ...
    devFromXmean, sameColorBars, paramOrder, paramName, coeffStruct, numBiasesStruct)

[X, Y, xname, yname] = findSurfXY(altitude, lat, lst, lon, doy, paramOrder);

averLsts = 0:1:23;
averLats = -85:10:85;
averLons = -180:10:170;
if lstMean
    N = length(X)*length(Y)*length(averLsts);
elseif latitudeMean
    N = length(X)*length(Y)*length(averLats);
elseif lonMean
    N = length(X)*length(Y)*length(averLons);
else
    N = length(X)*length(Y);
end

%S.aeInt = 20 * ones(N, 9);
S.aeInt = bsxfun(@times, ones(N,size(aeInt,2)), aeInt);
S.latitude = zeros(N, 1);
S.longitude = zeros(N, 1);
S.solarTime = zeros(N,1);
S.altitude = zeros(N,1);
S.F = F * ones(N,1);
S.FA = FA * ones(N,1);
S.doy = zeros(N,1);

if lstMean
    [xmat, ymat, zmat] = meshgrid(X, Y, averLsts);
    S.solarTime = zmat(:);
elseif latitudeMean
    [xmat, ymat, zmat] = meshgrid(X, Y, averLats);
    S.latitude = zmat(:);
elseif lonMean
    [xmat, ymat, zmat] = meshgrid(X, Y, averLons);
    S.longitude = zmat(:);
else
    [xmat, ymat] = meshgrid(X, Y);
    if length(lst) == 1; S.solarTime(:) = lst; end
    if length(lon) == 1; S.longitude(:) = lon; end
    if length(lat) == 1; S.latitude(:) = lat; end
end

if length(doy) == 1; S.doy(:) = doy; end
if length(altitude) == 1; S.altitude(:) = altitude; end

[S, titleAddition] = assignPlotVars(S, xmat(:), ymat(:), xname, yname, lstMean, lonMean, latitudeMean);

%TODO: ap:t parametrina
S.Ap = Ap*ones(N,1); S.apNow = Ap*ones(N,1); S.ap3h = Ap*ones(N,1); S.ap6h = Ap*ones(N,1);
S.ap9h = Ap*ones(N,1); S.ap12To33h = Ap*ones(N,1); S.ap36To57h = Ap*ones(N,1);

% [lstGrid, latGrid] = meshgrid(lst, lat);
% S = computeLatLstGrid(S, lat, lst);
% S.numBiases = 0;
S = computeVariablesForFit(S, false);
S = computeGeopotentialHeight(S);

TexCoeffs = coeffStruct.TexCoeff; dTCoeffs = coeffStruct.dTCoeff;
T0Coeffs = coeffStruct.T0Coeff;

if ~strcmpi(paramName, 'Tex') && ~strcmpi(paramName, 'T0') && ~strcmpi(paramName, 'dT')
    [Tex, dT0, T0] = findTempsForFit(S, TexCoeffs, dTCoeffs, T0Coeffs);
    OlbDens = evalMajorSpecies(S, coeffStruct.OCoeff, numBiasesStruct.O);
    N2lbDens = evalMajorSpecies(S, coeffStruct.N2Coeff, numBiasesStruct.N2);
    HelbDens = evalMajorSpecies(S, coeffStruct.HeCoeff, numBiasesStruct.He);
    ArlbDens = evalMajorSpecies(S, coeffStruct.ArCoeff, numBiasesStruct.Ar);
    O2lbDens = exp(coeffStruct.O2Coeff);
end

if strcmpi(paramName, 'Tex')
    param = evalTex(S, TexCoeffs);
    msisParam = computeMsis(S);
    dtmParam = computeDtm(S);
elseif strcmpi(paramName, 'T0')
    param = evalT0(S, T0Coeffs);
    [msisParam, ~, dtmParam, ~] = computeMsisDtmLb(S);
    msisParam = (msisParam-mean(msisParam))*2 + mean(msisParam);
    dtmParam = (dtmParam-mean(dtmParam))*2 + mean(dtmParam);
elseif strcmpi(paramName, 'dT')
    param = evalDT(S, dTCoeffs);
    [~, msisParam, ~, dtmParam] = computeMsisDtmLb(S);
elseif strcmpi(paramName, 'rho')
    param = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);
    [~,msisParam] = computeMsis(S);
    [~,dtmParam] = computeDtm(S);
elseif strcmpi(paramName, 'O')
    [~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);
    [~,~,msisParam] = computeMsis(S);
    [~,~,dtmParam] = computeDtm(S);
elseif strcmpi(paramName, 'N2')
    [~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);
    [~,~,~,msisParam] = computeMsis(S);
    [~,~,~,dtmParam] = computeDtm(S);
elseif strcmpi(paramName, 'He')
    [~,~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);   
    [~,~,~,~,msisParam] = computeMsis(S);
    [~,~,~,~,dtmParam] = computeDtm(S);
elseif strcmpi(paramName, 'Ar')
    [~,~,~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);   
    [~,~,~,~,~,msisParam] = computeMsis(S);
    dtmParam = zeros(size(msisParam));
elseif strcmpi(paramName, 'O2')
    [~,~,~,~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);   
    [~,~,~,~,~,~,msisParam] = computeMsis(S);
    dtmParam = zeros(size(msisParam));
else
    error(['Unknown variable: ', paramName])
end

if lstMean
    param = reshape(param, length(Y), length(X), length(averLsts));
    msisParam = reshape(msisParam, length(Y), length(X), length(averLsts));
    dtmParam = reshape(dtmParam, length(Y), length(X), length(averLsts));
    param = mean(param, 3); 
    msisParam = mean(msisParam, 3);
    dtmParam = mean(dtmParam, 3);
    xmat = xmat(:,:,1);
    ymat = ymat(:,:,1);
elseif lonMean
    param = reshape(param, length(Y), length(X), length(averLons));
    msisParam = reshape(msisParam, length(Y), length(X), length(averLons));
    dtmParam = reshape(dtmParam, length(Y), length(X), length(averLons));
    param = mean(param, 3); 
    msisParam = mean(msisParam, 3);
    dtmParam = mean(dtmParam, 3);
    xmat = xmat(:,:,1);
    ymat = ymat(:,:,1);
elseif latitudeMean
    param = reshape(param, length(Y), length(X), length(averLats));
    msisParam = reshape(msisParam, length(Y), length(X), length(averLats));
    dtmParam = reshape(dtmParam, length(Y), length(X), length(averLats));
    weights = ones(1,1,length(averLats)); weights(1,1,:) = cosd(averLats);
    param = bsxfun(@times, param, weights); param = sum(param, 3) / sum(weights, 3);
    msisParam = bsxfun(@times, msisParam, weights); msisParam = sum(msisParam, 3) / sum(weights, 3);
    dtmParam = bsxfun(@times, dtmParam, weights); dtmParam = sum(dtmParam, 3) / sum(weights, 3);
    xmat = xmat(:,:,1);
    ymat = ymat(:,:,1);
else
    param = reshape(param, length(Y), length(X));
    msisParam = reshape(msisParam, length(Y), length(X));
    dtmParam = reshape(dtmParam, length(Y), length(X));
end

if devFromXmean
    param = bsxfun(@rdivide, param, mean(param, 2));
    msisParam = bsxfun(@rdivide, msisParam, mean(msisParam, 2));
    dtmParam = bsxfun(@rdivide, dtmParam, mean(dtmParam, 2));
end

figure('renderer', 'zbuffer');

if sameColorBars
    subplot(3,1,1);
    %clims = plotSurfSubplot(lstGrid, latGrid, dtmParam, ['DTM ', paramName], FA, doy, heights(a), 16);
    clims = plotSurfSubplot(xmat, ymat, param, paramName, titleAddition, xname, yname, 16);

    subplot(3,1,2);
    plotSurfSubplot(xmat, ymat, msisParam, paramName, ['MSIS ',titleAddition], xname, yname, 16, clims);
    %plotSurfSubplot(xmat, ymat, msisParam, ['MSIS ', paramName], FA, doy, altitude(a), 16, clims);

    subplot(3,1,3);
    plotSurfSubplot(xmat, ymat, dtmParam, paramName, ['DTM ',titleAddition], xname, yname, 16, clims);
    %plotSurfSubplot(xmat, ymat, dtmParam, ['DTM ', paramName], FA, doy, altitude(a), 16, clims);
    % plotSurfSubplot(lstGrid, latGrid, param, paramName, FA, doy, heights(a), 16, clims);
else
    subplot(3,1,1);
    %clims = plotSurfSubplot(lstGrid, latGrid, dtmParam, ['DTM ', paramName], FA, doy, heights(a), 16);
    plotSurfSubplot(xmat, ymat, param, paramName, titleAddition, xname, yname, 16);

    subplot(3,1,2);
    plotSurfSubplot(xmat, ymat, msisParam, paramName, ['MSIS ',titleAddition], xname, yname, 16);
    %plotSurfSubplot(xmat, ymat, msisParam, ['MSIS ', paramName], FA, doy, altitude(a), 16, clims);

    subplot(3,1,3);
    plotSurfSubplot(xmat, ymat, dtmParam, paramName, ['DTM ',titleAddition], xname, yname, 16);
    %plotSurfSubplot(xmat, ymat, dtmParam, ['DTM ', paramName], FA, doy, altitude(a), 16, clims);
    % plotSurfSubplot(lstGrid, latGrid, param, paramName, FA, doy, heights(a), 16, clims);

end


end

function [] = analyzeStormTimes(rhoStruct, modelStruct, saveFolder, fullscreenFigs)

[stormBeginInd, stormEndInd, combinedInd, satInfo] = findStorms(rhoStruct, 'Dst', -50);
rawCorr = zeros(length(stormBeginInd),3);
rawOM = zeros(length(stormBeginInd),3);
rawRMS = zeros(length(stormBeginInd),3);
OACorr = zeros(length(stormBeginInd),3);
OAOM = zeros(length(stormBeginInd),3);
OARMS = zeros(length(stormBeginInd),3);

minDst = zeros(length(stormBeginInd),1);
averF81A = zeros(length(stormBeginInd),1);
t1 = cell(length(stormBeginInd),1);
t2 = cell(length(stormBeginInd),1);

for i = 1:length(stormBeginInd)
    ind = stormBeginInd(i):stormEndInd(i);
    
    lat = rhoStruct.latitude(ind);
    timestamps = rhoStruct.timestamps(ind);
    firstDay = timestamps < timestamps(1) + 1;
    
    measuredRho = rhoStruct.data(ind);
    ilRho = modelStruct.il(ind); 
    msisRho = modelStruct.msis(ind);
    dtmRho = modelStruct.dtm(ind);
    
    ilRho = ilRho * mean(measuredRho(firstDay)) / mean(ilRho(firstDay));
    msisRho = msisRho * mean(measuredRho(firstDay)) / mean(msisRho(firstDay));
    dtmRho = dtmRho * mean(measuredRho(firstDay)) / mean(dtmRho(firstDay));
    
    measuredOrbAver = computeOrbitAverage(measuredRho, lat, timestamps);
    ilOrbAver = computeOrbitAverage(ilRho, lat, timestamps);
    msisOrbAver = computeOrbitAverage(msisRho, lat, timestamps);
    dtmOrbAver = computeOrbitAverage(dtmRho, lat, timestamps);
    
    t1{i} = datestr(rhoStruct.timestamps(ind(1)));
    t2{i} = datestr(rhoStruct.timestamps(ind(end)));
    
    rawCorr(i,1) = corr(ilRho, measuredRho);
    rawCorr(i,2) = corr(msisRho, measuredRho);
    rawCorr(i,3) = corr(dtmRho, measuredRho);
    
    OACorr(i,1) = corr(ilOrbAver, measuredOrbAver);
    OACorr(i,2) = corr(msisOrbAver, measuredOrbAver);
    OACorr(i,3) = corr(dtmOrbAver, measuredOrbAver);
    
    rawOM(i,1) = mean(measuredRho ./ ilRho);
    rawOM(i,2) = mean(measuredRho ./ msisRho);
    rawOM(i,3) = mean(measuredRho ./ dtmRho);
    
    OAOM(i,1) = mean(measuredOrbAver ./ ilOrbAver);
    OAOM(i,2) = mean(measuredOrbAver ./ msisOrbAver);
    OAOM(i,3) = mean(measuredOrbAver ./ dtmOrbAver);
    
    rawRMS(i,1) = rms(measuredRho ./ ilRho-1);
    rawRMS(i,2) = rms(measuredRho ./ msisRho-1);
    rawRMS(i,3) = rms(measuredRho ./ dtmRho-1);
    
    OARMS(i,1) = rms(measuredOrbAver ./ ilOrbAver-1);
    OARMS(i,2) = rms(measuredOrbAver ./ msisOrbAver-1);
    OARMS(i,3) = rms(measuredOrbAver ./ dtmOrbAver-1);
    
    minDst(i) = min(rhoStruct.dst(ind));
    averF81A(i) = mean(rhoStruct.FA(ind));
end

outputCell = cell(length(stormBeginInd)+1, 14);
outputCell(2:end, 1) = t1;
outputCell(2:end, 2) = t2;
outputCell(2:end, 3) = num2cell(averF81A);
outputCell(2:end, 4) = num2cell(minDst);
outputCell(2:end, 5) = num2cell(satInfo);
outputCell(2:end, 6:8) = num2cell(rawCorr);
outputCell(2:end, 9:11) = num2cell(rawOM);
outputCell(2:end, 12:14) = num2cell(rawRMS);
outputCell(2:end, 15:17) = num2cell(OACorr);
outputCell(2:end, 18:20) = num2cell(OAOM);
outputCell(2:end, 21:23) = num2cell(OARMS);

outputCell(1,1:end) = {'Begin', 'End', 'Mean F10.7', 'Min Dst', 'SatInfo (0=GO:1=CH:2=GR:3=SW)',...
                        'Corr. IL','Corr. MSIS','Corr. DTM','O/M IL','O/M MSIS','O/M DTM',...
                        'RMS IL','RMS MSIS','RMS DTM', ...
                        'OA Corr. IL', 'OA Corr. MSIS', 'OA Corr. DTM',...
                        'OA O/M IL', 'OA O/M MSIS', 'OA O/M DTM',...
                        'OA RMS IL', 'OA RMS MSIS','OA RMS DTM'};

cell2csv([saveFolder, '/storms.csv'], outputCell);

fontsize = 15;
if fullscreenFigs
    figure('units','normalized','outerposition',[0 0 1 1])
else
    figure;
end
X = repmat(-minDst, 1, 2);
h = plot(X, OACorr(:,1:2), 's');
set(h(1), 'markerFaceColor', 'b');
set(h(2), 'markerFaceColor', 'g');
xlabel('-1 * min Dst', 'fontsize', fontsize);
ylabel('Orb. aver. Correlation', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
legend('AE', 'MSIS', 'location', 'southeast')
filename = [saveFolder, '/OACorr'];
saveas(gcf, filename, 'png');

if fullscreenFigs
    figure('units','normalized','outerposition',[0 0 1 1])
else
    figure;
end
X = repmat(-minDst, 1, 2);
h = plot(X, OARMS(:,1:2), 's');
set(h(1), 'markerFaceColor', 'b');
set(h(2), 'markerFaceColor', 'g');
xlabel('-1 * min Dst', 'fontsize', fontsize);
ylabel('Orb. aver. RMS', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
legend('AE', 'MSIS', 'location', 'southeast')
filename = [saveFolder, '/OARMS'];
saveas(gcf, filename, 'png');

if fullscreenFigs
    figure('units','normalized','outerposition',[0 0 1 1])
else
    figure;
end
X = repmat(-minDst, 1, 2);
h = plot(X, OAOM(:,1:2), 's');
set(h(1), 'markerFaceColor', 'b');
set(h(2), 'markerFaceColor', 'g');
xlabel('-1 * min Dst', 'fontsize', fontsize);
ylabel('Orb. aver. O/M', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
legend('AE', 'MSIS', 'location', 'southeast')
filename = [saveFolder, '/OAOM'];
saveas(gcf, filename, 'png');

end

function [S, titleAddition] = assignPlotVars(S, X, Y, xname, yname, lstMean, lonMean, latitudeMean)

if strcmpi(xname, 'altitude')
    S.altitude = X; xNum = 1;
elseif strcmpi(xname, 'lat')
    S.latitude = X; xNum = 2;
elseif strcmpi(xname, 'lst')
    S.solarTime = X; xNum = 3;
elseif strcmpi(xname, 'lon')
    S.longitude = X; xNum = 4;
elseif strcmpi(xname, 'doy')
    S.doy = X; xNum = 4;
end

if strcmpi(yname, 'altitude')
    S.altitude = Y; yNum = 1;
elseif strcmpi(yname, 'lat')
    S.latitude = Y; yNum = 2;
elseif strcmpi(yname, 'lst')
    S.solarTime = Y; yNum = 3;
elseif strcmpi(yname, 'lon')
    S.longitude = Y; yNum = 4;
elseif strcmpi(yname, 'doy')
    S.doy = Y; yNum = 4;
end

names = {'z', 'lat', 'lst', 'lon', 'doy'};
nonVecInd = setdiff(1:length(names), [xNum, yNum]);
nonVecNames = names(nonVecInd);
vals = [S.altitude(1), S.latitude(1), S.solarTime(1), S.longitude(1), S.doy(1)];
nonVecVals = vals(nonVecInd);
if lonMean
    otherInd = find(nonVecInd ~= 4);
    titleAddition = [names{otherInd}, '=', num2str(vals(otherInd)), ', lon mean'];
elseif lstMean
    otherInd = find(nonVecInd ~= 3);
    titleAddition = [names{otherInd}, '=', num2str(vals(otherInd)), ', zonal mean'];
elseif latitudeMean
    otherInd = find(nonVecInd ~= 2);
    titleAddition = [names{otherInd}, '=', num2str(vals(otherInd)), ', meridional mean'];
else
    titleAddition = [nonVecNames{1}, '=', num2str(nonVecVals(1)), ', ',...
                     nonVecNames{2}, '=', num2str(nonVecVals(2))];
end

end

function [X, Y, xname, yname] = findSurfXY(varargin)
% [X, Y] = findSurfXY(heights, lat, lst, doy, paramOrder)

lh = length(varargin{1}); llat = length(varargin{2}); llst = length(varargin{3});...
llon = length(varargin{4});    ldoy = length(varargin{5});
a = [lh, llat, llst, llon, ldoy];
names = {'altitude', 'lat', 'lst', 'lon', 'doy'};
if strcmpi(varargin{6}, 'xy')
    xind = find(a > 1, 1, 'first');
    yind = find(a > 1, 1, 'last');
else
    xind = find(a > 1, 1, 'last');
    yind = find(a > 1, 1, 'first');
end

X = varargin{xind};
Y = varargin{yind};
xname = names{xind};
yname = names{yind};

end

function [ilRho, msisRho, dtmRho] = computeComparisonData(rhoStruct, coeffStruct, numBiasesStruct)

rhoStruct = computeVariablesForFit(rhoStruct);
TexCoeffs = coeffStruct.TexCoeff; dTCoeffs = coeffStruct.dTCoeff;
T0Coeffs = coeffStruct.T0Coeff;
[Tex, dT0, T0] = findTempsForFit(rhoStruct, TexCoeffs, dTCoeffs, T0Coeffs);

OlbDens = evalMajorSpecies(rhoStruct, coeffStruct.OCoeff, numBiasesStruct.O);
N2lbDens = evalMajorSpecies(rhoStruct, coeffStruct.N2Coeff, numBiasesStruct.N2);
HelbDens = evalMajorSpecies(rhoStruct, coeffStruct.HeCoeff, numBiasesStruct.He);
ArlbDens = evalMajorSpecies(rhoStruct, coeffStruct.ArCoeff, numBiasesStruct.Ar);
O2lbDens = exp(coeffStruct.O2Coeff);

ilRho = computeRho(T0, dT0, Tex, rhoStruct.Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);

if nargout > 1
    [~, msisRho] = computeMsis(rhoStruct);
    [~, dtmRho] = computeDtm(rhoStruct);
end

end

function clims = plotSurfSubplot(xmat, ymat, param, paramName, titleAddition, xname, yname, fs, clims)

surf(xmat, ymat, param, 'edgecolor', 'none');
view(2);
colorbar;
%colormap jet
axis tight;
xlabel(xname, 'fontsize', fs)
ylabel(yname, 'fontsize', fs)
if (nargin == 9)
    caxis(clims)
end
clims = caxis;
title([paramName, ' ', titleAddition], 'fontsize', fs)
shading interp;
set(gca, 'fontsize', fs);

end

function [Tex, rho, O, N2, He] = computeDtm(S)

N = length(S.latitude);
Tex = zeros(N,1);
rho = zeros(N,1);
O = zeros(N,1);
N2 = zeros(N,1);
He = zeros(N,1);

dtm2013_mex();

targetCount = round(N / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running DTM, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

for i = 1:N
    [Tex(i),~,rho(i),~,~,d] = dtm2013_mex(S.doy(i), S.altitude(i), S.latitude(i), S.longitude(i), ...
        S.solarTime(i), S.F(i), S.FA(i), S.ap3h(i), S.Ap(i));
    O(i) = d(3); N2(i) = d(4); He(i) = d(2);
    if mod(i, 10000) == 0
        p.progress;
    end
end
p.stop;

rho = rho * 1E3;
u2g = 1.6605402e-24;
O = O / (16 * u2g);
N2 = N2 / (28 * u2g);
He = He / (4 * u2g);


end

function [] = computeStatistics(rhoStruct, ilRho, msisRho, dtmRho, saveFolder)

outputFile = fopen([saveFolder,'/','stat.out'], 'w');
fprintf(outputFile, '%s \n', 'Full model')
fprintf(outputFile, 'IL O/C: %f \n', mean(rhoStruct.data./ilRho));
fprintf(outputFile, 'MSIS O/C: %f \n', mean(rhoStruct.data./msisRho));
fprintf(outputFile, 'DTM O/C: %f \n\n', mean(rhoStruct.data./dtmRho));

fprintf(outputFile, 'IL RMS: %f \n', rms(rhoStruct.data./ilRho-1));
fprintf(outputFile, 'MSIS RMS: %f \n', rms(rhoStruct.data./msisRho-1));
fprintf(outputFile, 'DTM RMS: %f \n\n', rms(rhoStruct.data./dtmRho-1));

fprintf(outputFile, 'IL STD: %f \n', std(rhoStruct.data./ilRho));
fprintf(outputFile, 'MSIS STD: %f \n', std(rhoStruct.data./msisRho));
fprintf(outputFile, 'DTM STD: %f \n\n', std(rhoStruct.data./dtmRho));

fprintf(outputFile, 'IL CORR: %f \n', corr(rhoStruct.data,ilRho));
fprintf(outputFile, 'MSIS CORR: %f \n', corr(rhoStruct.data,msisRho));
fprintf(outputFile, 'DTM CORR: %f \n\n', corr(rhoStruct.data,dtmRho));

end

function [] = plot3DOM(x, dx, y, dy, obs, modelStruct, rmsOrOM, xName, yName, saveFolder, fullscreenFigs)

minX = dx * floor(min(x) / dx); minY = dy * floor(min(y) / dy);
maxX = dx * ceil(max(x) / dx); maxY = dy * ceil(max(y) / dy);

[Xmat, Ymat] = meshgrid(minX:dx:maxX, minY:dy:maxY);
countMat = zeros(size(Xmat)-1);
ilMat = zeros(size(Xmat)-1);
dtmMat = zeros(size(Xmat)-1);
msisMat = zeros(size(Xmat)-1);

ilVal = obs ./ modelStruct.il;
dtmVal = obs ./ modelStruct.dtm;
msisVal = obs ./ modelStruct.msis;

if strcmpi(rmsOrOM, 'rms')
    computeRms = true;
    rmsOrOM = 'RMS';
else
    computeRms = false;
    rmsOrOM = 'O/M';
end

for j = 1:size(Xmat,2)-1
    xInd = Xmat(1,j) <= x & x < Xmat(1,j+1);
    for i = 1:size(Xmat,1)-1       
        yInd = Ymat(i,1) <= y & y < Ymat(i+1,1);
        ind = xInd & yInd;
        countMat(i,j) = sum(ind);
        if computeRms
            ilMat(i,j) = rms(ilVal(ind) - 1);
            dtmMat(i,j) = rms(dtmVal(ind) - 1);
            msisMat(i,j) = rms(msisVal(ind) - 1);
        else
            ilMat(i,j) = mean(ilVal(ind));
            dtmMat(i,j) = mean(dtmVal(ind));
            msisMat(i,j) = mean(msisVal(ind));
        end
    end
end

Xmat = Xmat + dx/2; Xmat = Xmat(1:end-1,1:end-1);
Ymat = Ymat + dy/2; Ymat = Ymat(1:end-1,1:end-1);

baseTitle = [' ', rmsOrOM];

if fullscreenFigs
    figure('units','normalized','outerposition',[0 0 1 1])
else
    figure;
end
subplot(4,1,1)
surf(Xmat, Ymat, log10(countMat), 'edgecolor', 'none')
view(2);
axis tight;
h = colorbar;
axisTicks = round(10.^get(h,'YTick'));
tickLabels = strsplit(num2str(axisTicks));
set(h,'YTickLabel', tickLabels);
Title = ['Number of observations'];
title(Title);
xlabel(xName)
ylabel(yName)

subplot(4,1,2)
surf(Xmat, Ymat, ilMat, 'edgecolor', 'none')
view(2);
axis tight;
colorbar;
Title = ['IL ', baseTitle];
title(Title);
xlabel(xName)
ylabel(yName)
clims = caxis;

subplot(4,1,3)
surf(Xmat, Ymat, msisMat, 'edgecolor', 'none')
view(2);
axis tight;
colorbar;
Title = ['MSIS ', baseTitle];
title(Title);
xlabel(xName)
ylabel(yName)
caxis(clims);

subplot(4,1,4)
surf(Xmat, Ymat, dtmMat, 'edgecolor', 'none')
view(2);
axis tight;
colorbar;
Title = ['DTM ', baseTitle];
title(Title);
xlabel(xName)
ylabel(yName)
caxis(clims);

if strcmpi(rmsOrOM, 'O/M')
    rmsOrOM = 'OM';
end
filename = [saveFolder,'/','3D',xName,yName,rmsOrOM];
saveas(gcf, filename, 'png')

end

function [] = plot2DOM(x, dx, obs, modelStruct, rmsOrOM, xName, saveFolder, fullscreenFigs)

minX = dx * floor(min(x) / dx);
maxX = dx * ceil(max(x) / dx);

Xbins = minX:dx:maxX;
countVec = zeros(length(Xbins)-1,1);
ilVec = zeros(length(Xbins)-1,1);
dtmVec = zeros(length(Xbins)-1,1);
msisVec = zeros(length(Xbins)-1,1);

ilVal = obs ./ modelStruct.il;
dtmVal = obs ./ modelStruct.dtm;
msisVal = obs ./ modelStruct.msis;

if strcmpi(rmsOrOM, 'rms')
    computeRms = true;
    rmsOrOM = 'RMS';
else
    computeRms = false;
    rmsOrOM = 'O/M';
end

for i = 1:length(Xbins)-1
    ind = Xbins(i) <= x & x < Xbins(i+1);
    countVec(i) = sum(ind);
    if computeRms
        ilVec(i) = rms(ilVal(ind) - 1);
        dtmVec(i) = rms(dtmVal(ind) - 1);
        msisVec(i) = rms(msisVal(ind) - 1);
    else
        ilVec(i) = mean(ilVal(ind));
        dtmVec(i) = mean(dtmVal(ind));
        msisVec(i) = mean(msisVal(ind));
    end
end

Xbins = Xbins + dx/2; Xbins = Xbins(1:end-1);
plotX = repmat(Xbins', 1, 2);
plotY = [msisVec, dtmVec];

interpX = minX:dx/100:maxX;
interpCount = interp1(Xbins, countVec, interpX, 'nearest', 'extrap');

fontsize = 13;

if fullscreenFigs
    figure('units','normalized','outerposition',[0 0 1 1])
else
    figure;
end
[hax, hIL, hCount] = plotyy(Xbins, ilVec, interpX, interpCount);
set(hCount, 'color', 'k', 'linestyle', '-');
set(hIL, 'linewidth', 2.0);
hold all
plot(hax(1), plotX, plotY)

legend('AE', 'MSIS', 'DTM', '# obs.');
title(['Model ', rmsOrOM, ' vs. ', xName], 'fontsize', fontsize)
xlabel(xName, 'fontsize', fontsize)
ylabel(rmsOrOM, 'fontsize', fontsize)
set(hax(1), 'fontsize', fontsize, 'ycolor', 'k')
set(hax(2), 'fontsize', fontsize, 'ycolor', 'k')

if strcmpi(rmsOrOM, 'O/M')
    rmsOrOM = 'OM';
end
filename = [saveFolder,'/','2D',xName,rmsOrOM];
saveas(gcf, filename, 'png');

end

%function [] = plotLatCrossSection(paramName, AE, F, FA, lst, doy, )

function [] = plotStormFig(rhoStruct, modelStruct, date1, date2, satellite, coeffStruct, numBiasesStruct, saveFolder, fullscreenFigs)

t1 = datenum(date1);
t2 = datenum(date2);
t = rhoStruct.timestamps;
tInd = (t1 <= t & t <= t2);

if strcmpi(satellite, 'CHAMP')
    logInd = ismember(1:length(rhoStruct.data), rhoStruct.champ);
elseif strcmpi(satellite, 'GRACE')
    logInd = ismember(1:length(rhoStruct.data), rhoStruct.grace);
elseif strcmpi(satellite, 'GOCE')
    logInd = ismember(1:length(rhoStruct.data), rhoStruct.goce);
elseif strcmpi(satellite, 'SWARM')
    logInd = ismember(1:length(rhoStruct.data), rhoStruct.swarm);
else
    error(['Unrecognized satellite: ', satellite]);
end

ind = (tInd & logInd');
timestamps = t(ind);
lat = rhoStruct.latitude(ind);
lon = rhoStruct.longitude(ind);
alt = rhoStruct.altitude(ind);
measured = rhoStruct.data(ind);
ilRho = modelStruct.il(ind);
msisRho = modelStruct.msis(ind);
dtmRho = modelStruct.dtm(ind);

os = 95*60 / (mode(round(diff(timestamps*86400))));
if mod(os, 2) == 0; os = os + 1; end

if fullscreenFigs
    figure('units','normalized','outerposition',[0 0 1 1])
else
    figure;
end
plot(timestamps, smooth(measured,os), timestamps, smooth(ilRho,os), timestamps, smooth(msisRho,os), timestamps, smooth(dtmRho,os));
legend(satellite, 'OUR', 'MSIS', 'DTM')
datetick('x')
ylabel('Rho [kg/m^3]')
filename = [saveFolder,'/','2D',satellite,date1];
saveas(gcf, filename, 'png');

outputFile = fopen([saveFolder,'/','stat.out'], 'a');
fprintf('IL CORR: %f \n', corr(measured,ilRho));
fprintf('MSIS CORR: %f \n', corr(measured,msisRho));
fprintf('DTM CORR: %f \n\n', corr(measured,dtmRho));

fprintf(outputFile, '%s from %s to %s\n', satellite, date1, date2);
fprintf(outputFile, 'IL CORR: %f \n', corr(measured,ilRho));
fprintf(outputFile, 'MSIS CORR: %f \n', corr(measured,msisRho));
fprintf(outputFile, 'DTM CORR: %f \n\n', corr(measured,dtmRho));
fprintf(outputFile, 'IL O/M: %f \n', mean(measured./ilRho));
fprintf(outputFile, 'MSIS O/M: %f \n', mean(measured./msisRho));
fprintf(outputFile, 'DTM O/M: %f \n\n', mean(measured./dtmRho));
fprintf(outputFile, 'IL RMS: %f \n', rms(measured./ilRho-1));
fprintf(outputFile, 'MSIS RMS: %f \n', rms(measured./msisRho-1));
fprintf(outputFile, 'DTM RMS: %f \n\n', rms(measured./dtmRho-1));

aeIntegral = rhoStruct.aeInt(ind,5);
timestamps = (timestamps-t1)*86400;
removeInd = ~ind;
rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, true, false);

[ilRhoConstAlt, msisRhoConstAlt, dtmRhoConstAlt] = computeComparisonData(rhoStruct, coeffStruct, numBiasesStruct);

correctedDensity = measured .* (msisRho./msisRhoConstAlt);

lat = convertToMagneticCoordinates(lat, lon, alt);

i = rhoStruct.solarTime <= 12;
[limitedTimestamps, limitedLatitude, minAllowedLatitude, maxAllowedLatitude] = giveExactOrbits(timestamps(i), lat(i), false);
interpolateAndPlotByLatitude(t1, aeIntegral(i), timestamps(i), timestamps(i), lat(i), ...
    correctedDensity(i), msisRhoConstAlt(i), dtmRhoConstAlt(i), ilRhoConstAlt(i), limitedLatitude, limitedTimestamps, minAllowedLatitude, maxAllowedLatitude, 'morning', fullscreenFigs)

i = rhoStruct.solarTime > 12;
[limitedTimestamps, limitedLatitude, minAllowedLatitude, maxAllowedLatitude] = giveExactOrbits(timestamps(i), lat(i), false);
interpolateAndPlotByLatitude(t1, aeIntegral(i), timestamps(i), timestamps(i), lat(i), ...
    correctedDensity(i), msisRhoConstAlt(i), dtmRhoConstAlt(i), ilRhoConstAlt(i), limitedLatitude, limitedTimestamps, minAllowedLatitude, maxAllowedLatitude, 'evening', fullscreenFigs)

filename = [saveFolder,'/','3D',satellite,date1];
saveas(gcf, filename, 'png');

end

function [T0_msis, dT_msis, T0_dtm, dT_dtm] = computeMsisDtmLb(S)

z0 = S.altitude(1);
dz = 0.1;

N = length(S.latitude);
T0_msis = zeros(N, 1);
dT_msis = zeros(N, 1);
T0_dtm = zeros(N, 1);
dT_dtm = zeros(N, 1);

targetCount = round(N / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running MSIS, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
dtm2013_mex();

for i = 1:N
    [~, ~, ~,~,~,~,~,~,~,~,T1] = nrlmsise_mex(S.doy(i),43200,z0-dz,S.latitude(i),S.longitude(i),S.solarTime(i),...
        S.FA(i),S.F(i),S.Ap(i),S.apNow(i),S.ap3h(i),S.ap6h(i),S.ap9h(i),S.ap12To33h(i),S.ap36To57h(i));
    [~, ~, ~,~,~,~,~,~,~,~,T2] = nrlmsise_mex(S.doy(i),43200,z0,S.latitude(i),S.longitude(i),S.solarTime(i),...
        S.FA(i),S.F(i),S.Ap(i),S.apNow(i),S.ap3h(i),S.ap6h(i),S.ap9h(i),S.ap12To33h(i),S.ap36To57h(i));
    [~, ~, ~,~,~,~,~,~,~,~,T3] = nrlmsise_mex(S.doy(i),43200,z0+dz,S.latitude(i),S.longitude(i),S.solarTime(i),...
        S.FA(i),S.F(i),S.Ap(i),S.apNow(i),S.ap3h(i),S.ap6h(i),S.ap9h(i),S.ap12To33h(i),S.ap36To57h(i));
    dT_msis(i) = (T3-T1) / (2*dz);
    T0_msis(i) = T2;
    
    [~,T0,~,~,~,~,dT0] = dtm2013_mex(S.doy(i), z0, S.latitude(i), S.longitude(i), ...
    S.solarTime(i), S.F(i), S.FA(i), S.ap3h(i), S.Ap(i));

    dT_dtm(i) = dT0;
    T0_dtm(i) = T0;
    
    if mod(i, 10000) == 0
        p.progress;
    end
end
p.stop;

end


function plotDensityLatitudeTimeSurf(firstDatenum, aeIntegral, timestamps1min, magneticLatitude, timestamps10s, regriddedLatitude, regriddedTime, regriddedSatDensity, regriddedMsisDensity, regriddedDtmDensity, ...
    regriddedAeProxy, timeOfDay, fullscreenFigs)
% plotDensityLatitudeTimeSurf(averagedDensity, averagedLatitude, timestamps

persistent colormapFigHandle
persistent minDensity
persistent maxDensity

numPlotRows = 3; % 2 or 3
if ~isempty(strfind(lower(timeOfDay), 'morning'))
    %colormapFigHandle = figure('Color', 'white', 'units','normalized','outerposition',[0 0 1 1]);
    if fullscreenFigs
        colormapFigHandle = figure('units','normalized','outerposition',[0 0 1 1]);
    else
        colormapFigHandle = figure;
    end
    satSubplot = 1;
    if numPlotRows == 3
        msisDensitySubplot = 3;
        %dtmDensitySubplot = 5;
        aeProxySubplot = 5;
    end
    minDensity = min(regriddedSatDensity(:));
    maxDensity = max(regriddedSatDensity(:));
else
    satSubplot = 2;
    if numPlotRows == 3
        msisDensitySubplot = 4;
        %dtmDensitySubplot = 6;
        aeProxySubplot = 6;
    end
end
secondsInDay = 60 * 60 * 24;
minDensityTime = min(regriddedTime(:));
[minLat, maxLat] = findInterpolationLimits(magneticLatitude);
indicesToRemove = findMatrixIndicesInDatagap(regriddedTime, timestamps10s) | (regriddedTime < secondsInDay + minDensityTime);
regriddedTime(indicesToRemove) = nan(1);

regriddedTime = regriddedTime / secondsInDay + firstDatenum;
referenceDay = datestr(min(regriddedTime(:)), 'mmmm dd, yyyy');
regriddedTime = regriddedTime - datenum(referenceDay, 'mmmm dd, yyyy');

timestamps1min = timestamps1min / secondsInDay + firstDatenum - datenum(referenceDay, 'mmmm dd, yyyy');

minDensityTime = min(regriddedTime(:));
maxDensityTime = max(regriddedTime(:));
firstDayIndices = regriddedTime < minDensityTime + 1;
firstDaySat = regriddedSatDensity(firstDayIndices);
firstDayMsis = regriddedMsisDensity(firstDayIndices);
firstDayAe = regriddedAeProxy(firstDayIndices);
firstDayDtm = regriddedDtmDensity(firstDayIndices);
msisMultiplier = mean(firstDaySat(:) ./ firstDayMsis(:));
aeMultiplier = mean(firstDaySat(:) ./ firstDayAe(:));
dtmMultiplier = mean(firstDaySat(:) ./ firstDayDtm(:));

plotIndices = (regriddedTime >= minDensityTime + 0.0 & regriddedTime <= maxDensityTime);
[plotRows, ~] = find(plotIndices);
plotRows = unique(plotRows);

regriddedSatDensity = regriddedSatDensity(plotRows, :);
regriddedTime = regriddedTime(plotRows, :);
regriddedLatitude = regriddedLatitude(plotRows, :);
regriddedMsisDensity = msisMultiplier * regriddedMsisDensity(plotRows, :);
regriddedAeProxy = aeMultiplier * regriddedAeProxy(plotRows, :);
regriddedDtmDensity = dtmMultiplier * regriddedDtmDensity(plotRows, :);

plotHeight = maxDensity;

minDensityTime = minDensityTime + 0.0;
maxDensityTime = maxDensityTime - 0.0;
[~, minAeIndex] = min(abs(timestamps1min - minDensityTime));
[~, maxAeIndex] = min(abs(timestamps1min - maxDensityTime));
aeIndicesToPlot = minAeIndex:maxAeIndex;

figure(colormapFigHandle);

subplotAxesHandle = subplot(numPlotRows,2,satSubplot);
surf(subplotAxesHandle, regriddedTime, regriddedLatitude, regriddedSatDensity, 'EdgeColor', 'None')
xlim([minDensityTime maxDensityTime]);
ylim([minLat maxLat]);
caxis([minDensity maxDensity])
view(2);
colormap jet
colorbar('Location', 'EastOutside');
ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
title(['Measured ', timeOfDay,' density'], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
set(gca, 'fontsize', 12)

% hold all;
% aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
% aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
% set(aeLineHandle, 'LineWidth', 0.1)
% view(2);
% set(aeAxesHandle, 'yaxislocation', 'right');
% ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
% set(aeAxesHandle, 'Color', 'none', 'XTick', []);
% hold off;
% set(gca, 'fontsize', 12)


subplotAxesHandle = subplot(numPlotRows,2,msisDensitySubplot);
surf(subplotAxesHandle, regriddedTime, regriddedLatitude, regriddedMsisDensity, 'EdgeColor', 'None')
xlim([minDensityTime maxDensityTime]);
ylim([minLat maxLat]);
caxis([minDensity maxDensity])
colorbar('Location', 'EastOutside');
view(2);
colormap jet
ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
title(['NRLMSISE-00 ', timeOfDay,' density'], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
set(gca, 'fontsize', 12)

% hold all;
% aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
% aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
% set(aeLineHandle, 'LineWidth', 0.1)
% view(2);
% set(aeAxesHandle, 'yaxislocation', 'right');
% ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
% set(aeAxesHandle, 'Color', 'none', 'XTick', []);
% hold off;
% set(gca, 'fontsize', 12)


% subplotAxesHandle = subplot(numPlotRows,2,dtmDensitySubplot);
% surf(subplotAxesHandle, regriddedTime, regriddedLatitude, regriddedDtmDensity, 'EdgeColor', 'None')
% xlim([minDensityTime maxDensityTime]);
% ylim([minLat maxLat]);
% caxis([minDensity maxDensity])
% colorbar('Location', 'EastOutside');
% view(2);
% colormap jet
% ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
% title(['DTM-2013 ', timeOfDay,' density'], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
% set(gca, 'fontsize', 12)

% hold all;
% aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
% aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
% set(aeLineHandle, 'LineWidth', 0.1)
% view(2);
% set(aeAxesHandle, 'yaxislocation', 'right');
% ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
% set(aeAxesHandle, 'Color', 'none', 'XTick', []);
% hold off;
% set(gca, 'fontsize', 12)


subplotAxesHandle = subplot(numPlotRows,2,aeProxySubplot);
surf(subplotAxesHandle, regriddedTime, regriddedLatitude, regriddedAeProxy, 'EdgeColor', 'None')
xlim([minDensityTime maxDensityTime]);
ylim([minLat maxLat]);
caxis([minDensity maxDensity])
colorbar('Location', 'EastOutside');
view(2);
colormap jet
xlabel(['Days since the UTC beginning of ', referenceDay], 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
ylabel('Geomagnetic lat.', 'fontsize', 14, 'fontname', 'Courier', 'fontweight', 'bold')
title(['AE Predicted ', timeOfDay,' density'], 'fontsize', 13, 'fontname', 'courier', 'fontweight', 'bold')
set(gca, 'fontsize', 12)

% hold all;
% aeAxesHandle = axes('Position', get(subplotAxesHandle, 'Position'));
% aeLineHandle = plot3(aeAxesHandle, timestamps1min(aeIndicesToPlot), aeIntegral(aeIndicesToPlot), ones(size(aeIndicesToPlot)) * plotHeight, 'k');
% set(aeLineHandle, 'LineWidth', 0.1)
% view(2);
% set(aeAxesHandle, 'yaxislocation', 'right');
% ylabel(aeAxesHandle, 'AE 21-h Integral', 'fontsize', 12, 'fontname', 'courier', 'fontweight', 'bold')
% set(aeAxesHandle, 'Color', 'none', 'XTick', []);
% hold off;
% set(gca, 'fontsize', 12)

end

function [] = interpolateAndPlotByLatitude(firstDatenum, aeIntegral, timestamps1min, timestamps10s, magneticLatitude, ...
    correctedDensity, msisDensity, dtmDensity, aeProxyDensity, limitedLatitude, limitedTimestamps, minAllowedLatitude, maxAllowedLatitude, timeOfDay, fullscreenFigs)
%

oneQuarterDegreeStep = minAllowedLatitude+5:0.25:maxAllowedLatitude-5;
for i = 1:length(oneQuarterDegreeStep)
    regriddedTime(:,i) = latitudeCrossingTimes(limitedLatitude, limitedTimestamps, oneQuarterDegreeStep(i));    
end

regriddedSatDensity = interp1(timestamps10s, correctedDensity, regriddedTime, 'spline');
regriddedMsisDensity = interp1(timestamps10s, msisDensity, regriddedTime, 'spline');
regriddedDtmDensity = interp1(timestamps10s, dtmDensity, regriddedTime, 'spline');
regriddedAeProxy = interp1(timestamps10s, aeProxyDensity, regriddedTime, 'spline');

crossingTimes = regriddedTime(:,1:4:end);
satDensityByLatitude = regriddedSatDensity(:,1:4:end);
msisDensityByLatitude = regriddedMsisDensity(:,1:4:end);
dtmDensityByLatitude = regriddedDtmDensity(:,1:4:end);
aeProxyDensityByLatitude = regriddedAeProxy(:,1:4:end);


numOfOrbits = length(regriddedTime(:,1));
numOfValuesInOrbit = length(regriddedTime(1,:));
for i = 1:numOfValuesInOrbit
    timeThisLatitude = regriddedTime(:,i);
    goceDensityThisLatitude = regriddedSatDensity(:,i);
    msisDensityThisLatitude = regriddedMsisDensity(:,i);
    dtmDensityThisLatitude = regriddedDtmDensity(:,i);
    aeProxyThisLatitude = regriddedAeProxy(:,i);

    tInterp = interp1(1:numOfOrbits, timeThisLatitude, 1:1/20:numOfOrbits);
    interpolatedGoceDensity = interp1(timeThisLatitude, goceDensityThisLatitude, tInterp, 'spline');
    interpolatedMsisDensity = interp1(timeThisLatitude, msisDensityThisLatitude, tInterp, 'spline');
    interpolatedDtmDensity = interp1(timeThisLatitude, dtmDensityThisLatitude, tInterp, 'spline');
    interpolatedAeProxy = interp1(timeThisLatitude, aeProxyThisLatitude, tInterp, 'spline');

    latitudeMatrix(:,i) = ones(length(tInterp), 1) * oneQuarterDegreeStep(i); 
    goceDensityMatrix(:,i) = interpolatedGoceDensity;
    msisDensityMatrix(:,i) = interpolatedMsisDensity;
    dtmDensityMatrix(:,i) = interpolatedDtmDensity;
    aeProxyDensityMatrix(:,i) = interpolatedAeProxy;
    timeMatrix(:,i) = tInterp;
end

plotDensityLatitudeTimeSurf(firstDatenum, aeIntegral, timestamps1min, magneticLatitude, timestamps10s, latitudeMatrix, ...
    timeMatrix, goceDensityMatrix, msisDensityMatrix, dtmDensityMatrix, aeProxyDensityMatrix, timeOfDay, fullscreenFigs);

end