function [] = visualizeFit()
aeThreshold = 500;


% Check the existence of the data file.
if exist('ilData.mat', 'file')
    load('ilData.mat', 'rhoStruct')
else
    error('File ilData.mat not found!')
end

if exist('optCoeff.mat', 'file')
    load optCoeff.mat
else
    error('File optCoeff.mat not found!')
end

[rhoStruct] = removeAndFixData(rhoStruct, aeThreshold);

numBiasesStruct = struct('O', 4, 'N2', 5,...
    'He', 5, 'Ar', 2, 'O2', 0); % TODO: paivita arvot lisattyasi datasetteja

TexStruct.coeffInd = TexInd;
coeffStruct = struct('TexCoeff' , optCoeff(TexInd),... 
'OCoeff', optCoeff(OInd),...
'N2Coeff' , optCoeff(N2Ind),...
'HeCoeff' , optCoeff(HeInd),...
'O2Coeff' , optCoeff(O2Ind),...
'ArCoeff' , optCoeff(ArInd),...
'dTCoeff', dTCoeffs,...
'T0Coeff', T0Coeffs);

z = 800;
lat = -90:5:90;
lst = 0:0.5:24;
doy = 180;
F = 200;
FA = 200;
aeInt = 20*ones(1,9);
zonalMean = false;
latitudeMean = false;
devFromXmean = false;
sameColorBars = false;
%plotSurfs(z, lat, lst, doy, F, FA, aeInt, zonalMean, latitudeMean, devFromXmean, ...
%    sameColorBars, 'yx', 'rho', coeffStruct, numBiasesStruct);

if exist('comparisonRho.mat', 'file')
    load comparisonRho.mat
else
    [ilRho, msisRho, dtmRho] = computeComparisonData(rhoStruct, coeffStruct, numBiasesStruct);

    save('comparisonRho.mat', 'ilRho')
    save('comparisonRho.mat', 'msisRho', '-append')
    save('comparisonRho.mat', 'dtmRho', '-append')
end

modelStruct = struct('il', ilRho, 'msis', msisRho, 'dtm', dtmRho);

%plot3DOM(rhoStruct.aeInt(:,4), 50, rhoStruct.solarTime, 2, rhoStruct.data,...
%    modelStruct, 'O/M', 'AE 16h', 'lst');

%plot2DOM(rhoStruct.doy, 50, rhoStruct.data, modelStruct, 'O/M', 'DOY')

%computeStatistics(rhoStruct, ilRho, msisRho, dtmRho);

%plotStormFig(rhoStruct, ilRho, msisRho, dtmRho, '2003-10-27', '2003-11-02', 'CHAMP', TexCoeff, OCoeff, N2Coeff, HeCoeff);

end

function [Tex, dT0, T0] = findTempsForFit(varStruct, TexCoeffs, dTCoeffs, T0Coeffs)

Tex_est = evalTex(varStruct, TexCoeffs);

T0 = clamp(200, evalT0(varStruct, T0Coeffs), 1000);
dT0 = clamp(1, evalDT(varStruct, dTCoeffs), 30);
Tex = clamp(T0+1, Tex_est, 5000);

end

function [] = plotSurfs(altitude, lat, lst, doy, F, FA, aeInt, zonalMean, latitudeMean, ...
    devFromXmean, sameColorBars, paramOrder, paramName, coeffStruct, numBiasesStruct)

[X, Y, xname, yname] = findSurfXY(altitude, lat, lst, doy, paramOrder);

averLsts = 0:1:23;
averLats = -85:10:85;
if zonalMean
    N = length(X)*length(Y)*length(averLsts);
elseif latitudeMean
    N = length(X)*length(Y)*length(averLats);
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

if zonalMean
    [xmat, ymat, zmat] = meshgrid(X, Y, averLsts);
    S.solarTime = zmat(:);
elseif latitudeMean
    [xmat, ymat, zmat] = meshgrid(X, Y, averLats);
    S.latitude = zmat(:);
else
    [xmat, ymat] = meshgrid(X, Y);
    if length(lst) == 1; S.solarTime(:) = lst; end
end

if length(lat) == 1; S.latitude(:) = lat; end
S.longitude(:) = 0;
if length(doy) == 1; S.doy(:) = doy; end
if length(altitude) == 1; S.altitude(:) = altitude; end

[S, titleAddition] = assignPlotVars(S, xmat(:), ymat(:), xname, yname, zonalMean, latitudeMean);

%TODO: ap:t parametrina
S.Ap = 3*ones(N,1); S.apNow = 3*ones(N,1); S.ap3h = 3*ones(N,1); S.ap6h = 3*ones(N,1);
S.ap9h = 3*ones(N,1); S.ap12To33h = 3*ones(N,1); S.ap36To57h = 3*ones(N,1);

% [lstGrid, latGrid] = meshgrid(lst, lat);
% S = computeLatLstGrid(S, lat, lst);
% S.numBiases = 0;
S = computeVariablesForFit(S);
S = computeGeopotentialHeight(S);

TexCoeffs = coeffStruct.TexCoeff; dTCoeffs = coeffStruct.dTCoeff;
T0Coeffs = coeffStruct.T0Coeff;

if ~strcmpi(paramName, 'Tex') && ~strcmpi(paramName, 'T0') && ~strcmpi(paramName, 'dT')
    [Tex, dT0, T0] = findTempsForFit(S, TexCoeffs, dTCoeffs, T0Coeffs);
    OlbDens = evalMajorSpecies(S, coeffStruct.OCoeff, numBiasesStruct.O);
    N2lbDens = evalMajorSpecies(S, coeffStruct.N2Coeff, numBiasesStruct.N2);
    HelbDens = evalMajorSpecies(S, coeffStruct.HeCoeff, numBiasesStruct.He);
    ArlbDens = evalMinorSpecies(S, coeffStruct.ArCoeff, numBiasesStruct.Ar);
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

if zonalMean
    param = reshape(param, length(Y), length(X), length(averLsts));
    msisParam = reshape(msisParam, length(Y), length(X), length(averLsts));
    dtmParam = reshape(dtmParam, length(Y), length(X), length(averLsts));
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

figure;

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

function [S, titleAddition] = assignPlotVars(S, X, Y, xname, yname, zonalMean, latitudeMean)

if strcmpi(xname, 'altitude')
    S.altitude = X; xNum = 1;
elseif strcmpi(xname, 'lat')
    S.latitude = X; xNum = 2;
elseif strcmpi(xname, 'lst')
    S.solarTime = X; xNum = 3;
elseif strcmpi(xname, 'doy')
    S.doy = X; xNum = 4;
end

if strcmpi(yname, 'altitude')
    S.altitude = Y; yNum = 1;
elseif strcmpi(yname, 'lat')
    S.latitude = Y; yNum = 2;
elseif strcmpi(yname, 'lst')
    S.solarTime = Y; yNum = 3;
elseif strcmpi(yname, 'doy')
    S.doy = Y; yNum = 4;
end

names = {'z', 'lat', 'lst', 'doy'};
nonVecInd = setdiff(1:4, [xNum, yNum]);
nonVecNames = names(nonVecInd);
vals = [S.altitude(1), S.latitude(1), S.solarTime(1), S.doy(1)];
nonVecVals = vals(nonVecInd);
if zonalMean
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
    ldoy = length(varargin{4});
a = [lh, llat, llst, ldoy];
names = {'altitude', 'lat', 'lst', 'doy'};
if strcmpi(varargin{5}, 'xy')
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
ArlbDens = evalMinorSpecies(rhoStruct, coeffStruct.ArCoeff, numBiasesStruct.Ar);
O2lbDens = exp(coeffStruct.O2Coeff);

ilRho = computeRho(T0, dT0, Tex, rhoStruct.Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);

[~, msisRho] = computeMsis(rhoStruct);
[~, dtmRho] = computeDtm(rhoStruct);

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

function [] = computeStatistics(rhoStruct, ilRho, msisRho, dtmRho)

fprintf('IL O/C: %f \n', mean(rhoStruct.data./ilRho));
fprintf('MSIS O/C: %f \n', mean(rhoStruct.data./msisRho));
fprintf('DTM O/C: %f \n\n', mean(rhoStruct.data./dtmRho));

fprintf('IL RMS: %f \n', rms(rhoStruct.data./ilRho-1));
fprintf('MSIS RMS: %f \n', rms(rhoStruct.data./msisRho-1));
fprintf('DTM RMS: %f \n\n', rms(rhoStruct.data./dtmRho-1));

fprintf('IL STD: %f \n', std(rhoStruct.data./ilRho));
fprintf('MSIS STD: %f \n', std(rhoStruct.data./msisRho));
fprintf('DTM STD: %f \n\n', std(rhoStruct.data./dtmRho));

fprintf('IL CORR: %f \n', corr(rhoStruct.data,ilRho));
fprintf('MSIS CORR: %f \n', corr(rhoStruct.data,msisRho));
fprintf('DTM CORR: %f \n\n', corr(rhoStruct.data,dtmRho));

end

function [] = plot3DOM(x, dx, y, dy, obs, modelStruct, rmsOrOM, xName, yName)

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

figure;
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

end

function [] = plot2DOM(x, dx, obs, modelStruct, rmsOrOM, xName)

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

figure;
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

end

%function [] = plotLatCrossSection(paramName, AE, F, FA, lst, doy, )

function [] = plotStormFig(rhoStruct, ilRho, msisRho, dtmRho, date1, date2, satellite, TexCoeff, OCoeff, N2Coeff, HeCoeff)

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
end

ind = (tInd & logInd');
timestamps = t(ind);
lat = rhoStruct.latitude(ind);
lon = rhoStruct.longitude(ind);
alt = rhoStruct.altitude(ind);
measured = rhoStruct.data(ind);
ilRho = ilRho(ind);
msisRho = msisRho(ind);
dtmRho = dtmRho(ind);

os = 90*60 / (mode(round(diff(timestamps*86400))));
if mod(os, 2) == 0; os = os + 1; end

figure;
plot(timestamps, smooth(measured,os), timestamps, smooth(ilRho,os), timestamps, smooth(msisRho,os), timestamps, smooth(dtmRho,os));
legend(satellite, 'OUR', 'MSIS', 'DTM')
datetick('x')
ylabel('Rho [kg/m^3]')

fprintf('IL CORR: %f \n', corr(measured,ilRho));
fprintf('MSIS CORR: %f \n', corr(measured,msisRho));
fprintf('DTM CORR: %f \n\n', corr(measured,dtmRho));

aeIntegral = rhoStruct.aeInt(ind,5);
timestamps = (timestamps-t1)*86400;
removeInd = ~ind;
rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, true, false);

[ilRhoConstAlt, msisRhoConstAlt, dtmRhoConstAlt] = computeComparisonData(rhoStruct, TexCoeff, OCoeff, N2Coeff, HeCoeff);

correctedDensity = measured .* (msisRho./msisRhoConstAlt);

lat = convertToMagneticCoordinates(lat, lon, alt);

i = rhoStruct.solarTime <= 12;
[limitedTimestamps, limitedLatitude, minAllowedLatitude, maxAllowedLatitude] = giveExactOrbits(timestamps(i), lat(i), false);
interpolateAndPlotByLatitude(t1, aeIntegral(i), timestamps(i), timestamps(i), lat(i), ...
    correctedDensity(i), msisRhoConstAlt(i), dtmRhoConstAlt(i), ilRhoConstAlt(i), limitedLatitude, limitedTimestamps, minAllowedLatitude, maxAllowedLatitude, 'morning')

i = rhoStruct.solarTime > 12;
[limitedTimestamps, limitedLatitude, minAllowedLatitude, maxAllowedLatitude] = giveExactOrbits(timestamps(i), lat(i), false);
interpolateAndPlotByLatitude(t1, aeIntegral(i), timestamps(i), timestamps(i), lat(i), ...
    correctedDensity(i), msisRhoConstAlt(i), dtmRhoConstAlt(i), ilRhoConstAlt(i), limitedLatitude, limitedTimestamps, minAllowedLatitude, maxAllowedLatitude, 'evening')

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
    regriddedAeProxy, timeOfDay)
% plotDensityLatitudeTimeSurf(averagedDensity, averagedLatitude, timestamps

persistent colormapFigHandle
persistent minDensity
persistent maxDensity

numPlotRows = 3; % 2 or 3
if ~isempty(strfind(lower(timeOfDay), 'morning'))
    %colormapFigHandle = figure('Color', 'white', 'units','normalized','outerposition',[0 0 1 1]);
    colormapFigHandle = figure();
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
    correctedDensity, msisDensity, dtmDensity, aeProxyDensity, limitedLatitude, limitedTimestamps, minAllowedLatitude, maxAllowedLatitude, timeOfDay)
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
    timeMatrix, goceDensityMatrix, msisDensityMatrix, dtmDensityMatrix, aeProxyDensityMatrix, timeOfDay);

end