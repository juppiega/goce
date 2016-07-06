function [] = visualizeFit()
aeThreshold = 0;


if exist('optCoeff.mat', 'file')
    load optCoeff.mat
else
    error('File optCoeff.mat not found!')
end

% Check the existence of the data file.
if exist('ilData.mat', 'file')
    load ilData.mat
else
    error('File ilData.mat not found!')
end

[TempStruct, OStruct, N2Struct, HeStruct, rhoStruct, lbDTStruct, lbT0Struct] = ...
    removeAndFixData(TempStruct, OStruct, N2Struct, HeStruct, rhoStruct, lbDTStruct, lbT0Struct, aeThreshold);

% TODO: Ei huomioi numBiases:a!
TexCoeff = optCoeff(TexInd); TexStruct.coeffInd = TexInd;
OCoeff = optCoeff(OInd);
N2Coeff = optCoeff(N2Ind);
HeCoeff = optCoeff(HeInd);

plotSurfs(130, 140, 140, 1, 'dT', optCoeff, TexStruct, TexCoeff, dTCoeffs, T0Coeffs, OCoeff, N2Coeff, HeCoeff);

if exist('comparisonRho.mat', 'file')
    load comparisonRho.mat
else
    [ilRho, msisRho, dtmRho] = computeComparisonData(rhoStruct, optCoeff, TexStruct, dTCoeffs, T0Coeffs, OCoeff, N2Coeff, HeCoeff);

    save('comparisonRho.mat', 'ilRho')
    save('comparisonRho.mat', 'msisRho', '-append')
    save('comparisonRho.mat', 'dtmRho', '-append')
end

%computeStatistics(rhoStruct, ilRho, msisRho, dtmRho);

%plotStormFig(rhoStruct, ilRho, msisRho, dtmRho, '2003-10-27', '2003-11-02', 'CHAMP', TexCoeff, OCoeff, N2Coeff, HeCoeff);

end

function [] = plotSurfs(heights, F, FA, doy, paramName, allCoeff, TexStruct, TexCoeff, dTCoeffs, T0Coeffs, OCoeff, N2Coeff, HeCoeff)

lat = -90:1:90;
lst = 0:0.1:23.9;
N = length(lat)*length(lst);

S.aeInt = 20 * ones(N, 9);
S.latitude = zeros(N, 1);
S.longitude = zeros(N, 1);
S.solarTime = zeros(N,1);
S.F = F * ones(N,1);
S.FA = FA * ones(N,1);
S.doy = doy * ones(N,1);

S.Ap = 3*ones(N,1); S.apNow = 3*ones(N,1); S.ap3h = 3*ones(N,1); S.ap6h = 3*ones(N,1);
S.ap9h = 3*ones(N,1); S.ap12To33h = 3*ones(N,1); S.ap36To57h = 3*ones(N,1);

[lstGrid, latGrid] = meshgrid(lst, lat);
S = computeLatLstGrid(S, lat, lst);
S.numBiases = 0;
S = computeVariablesForFit(S);

if ~strcmpi(paramName, 'Tex') && ~strcmpi(paramName, 'T0') && ~strcmpi(paramName, 'dT')
    [Tex, dT0, T0] = findTempsForFit(S, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
    OlbDens = evalMajorSpecies(S, OCoeff);
    N2lbDens = evalMajorSpecies(S, N2Coeff);
    HelbDens = evalMajorSpecies(S, HeCoeff);
end

for a = 1:length(heights)
    S.altitude = ones(N,1) * heights(a);
    S = computeGeopotentialHeight(S);

    if strcmpi(paramName, 'Tex')
        param = evalTex(S, TexCoeff);
        msisParam = computeMsis(S);
        dtmParam = computeDtm(S);
    elseif strcmpi(paramName, 'T0')
        param = evalT0(S, T0Coeffs);
        param = (param - mean(param))*0.5 + 0.8*mean(param);     
        [msisParam, ~, dtmParam, ~] = computeMsisDtmLb(S);
        msisParam = (msisParam-mean(msisParam))*2 + mean(msisParam);
        dtmParam = (dtmParam-mean(dtmParam))*2 + mean(dtmParam);
    elseif strcmpi(paramName, 'dT')
        param = evalDT(S, dTCoeffs);
        [~, msisParam, ~, dtmParam] = computeMsisDtmLb(S);
    elseif strcmpi(paramName, 'rho')
        param = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens);
        [~,msisParam] = computeMsis(S);
        [~,dtmParam] = computeDtm(S);
    elseif strcmpi(paramName, 'O')
        [~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens);
        [~,~,msisParam] = computeMsis(S);
        [~,~,dtmParam] = computeDtm(S);
    elseif strcmpi(paramName, 'N2')
        [~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens);
        [~,~,~,msisParam] = computeMsis(S);
        [~,~,~,dtmParam] = computeDtm(S);
    elseif strcmpi(paramName, 'He')
        [~,~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens);   
        [~,~,~,~,msisParam] = computeMsis(S);
        [~,~,~,~,dtmParam] = computeDtm(S);
    else
        error(['Unknown variable: ', paramName])
    end
    
    param = reshape(param, length(lat), length(lst));
    msisParam = reshape(msisParam, length(lat), length(lst));
    dtmParam = reshape(dtmParam, length(lat), length(lst));
    
    figure;
    
    subplot(3,1,1);
    clims = plotSurfSubplot(lstGrid, latGrid, dtmParam, ['DTM ', paramName], FA, doy, heights(a), 16);
    %clims = plotSurfSubplot(lstGrid, latGrid, param, paramName, FA, doy, heights(a), 16);
    
    subplot(3,1,2);
    plotSurfSubplot(lstGrid, latGrid, msisParam, ['MSIS ', paramName], FA, doy, heights(a), 16, clims);
    
    subplot(3,1,3);
    %plotSurfSubplot(lstGrid, latGrid, dtmParam, ['DTM ', paramName], FA, doy, heights(a), 16, clims);
    plotSurfSubplot(lstGrid, latGrid, param, paramName, FA, doy, heights(a), 16, clims);
    
end

end

function [ilRho, msisRho, dtmRho] = computeComparisonData(rhoStruct, allCoeff, TexStruct, dTCoeffs, T0Coeffs, OCoeff, N2Coeff, HeCoeff)

rhoStruct = computeVariablesForFit(rhoStruct);

[Tex, dT0, T0] = findTempsForFit(rhoStruct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
OlbDens = evalMajorSpecies(rhoStruct, OCoeff);
N2lbDens = evalMajorSpecies(rhoStruct, N2Coeff);
HelbDens = evalMajorSpecies(rhoStruct, HeCoeff);

ilRho = computeRho(T0, dT0, Tex, rhoStruct.Z, OlbDens, N2lbDens, HelbDens);

[~, msisRho] = computeMsis(rhoStruct);
[~, dtmRho] = computeDtm(rhoStruct);

end

function clims = plotSurfSubplot(lstGrid, latGrid, param, paramName, FA, doy, height, fs, clims)

surf(lstGrid, latGrid, param, 'edgecolor', 'none');
view(2);
colorbar;
%colormap jet
axis tight;
xlabel('LST', 'fontsize', fs)
ylabel('lat', 'fontsize', fs)
if (nargin == 9)
    caxis(clims)
end
clims = caxis;
if strcmpi(paramName, 'Tex')
    title(['Tex: F10.7=', num2str(FA), ', doy=', num2str(doy)], 'fontsize', fs)
else
    title([paramName, ': F10.7=', num2str(FA), ', doy=', num2str(doy), ', Z=', num2str(height)], 'fontsize', fs)
end
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

function S = computeLatLstGrid(S, lat, lst)

for j = 1:length(lst)
    for i = 1:length(lat)
        ind = i + (j-1)*length(lat);
        S.latitude(ind) = lat(i);
        S.solarTime(ind) = lst(j);
    end
end

end

function [Tex, dT0, T0] = findTemps(varStruct, TexCoeff)

origNumBiases = varStruct.numBiases; varStruct.numBiases = 0;
Tex_est = evalTex(varStruct, TexCoeff);
varStruct.numBiases = origNumBiases;

T0 = 507;
dT0 = 12.6;
Tex = clamp(T0+1, Tex_est, 5000);

end

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

z0 = 130;
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
    
    [~,T1] = dtm2013_mex(S.doy(i), z0-dz, S.latitude(i), S.longitude(i), ...
    S.solarTime(i), S.F(i), S.FA(i), S.ap3h(i), S.Ap(i));
    [~,T2] = dtm2013_mex(S.doy(i), z0, S.latitude(i), S.longitude(i), ...
    S.solarTime(i), S.F(i), S.FA(i), S.ap3h(i), S.Ap(i));
    [~,T3] = dtm2013_mex(S.doy(i), z0+dz, S.latitude(i), S.longitude(i), ...
    S.solarTime(i), S.F(i), S.FA(i), S.ap3h(i), S.Ap(i));
    dT_dtm(i) = (T3-T1) / (2*dz);
    T0_dtm(i) = T2;
    
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