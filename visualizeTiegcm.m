function [] = visualizeTiegcm(tiegcmFiles)

if(matlabpool('size')==0)
    matlabpool;
end

[timeseriesFig, axes, datenums] = plotTimeseries(tiegcmFiles);
plotSurfs(timeseriesFig, axes, datenums, tiegcmFiles)

end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function [figHandle, hAx, datenums] = plotTimeseries(tiegcmFileNames)

tiegcmFiles = dir(tiegcmFileNames);
lat = ncread(tiegcmFiles(1).name, 'lat');
lon = ncread(tiegcmFiles(1).name, 'lon');
cpcp = [];
hemisphericPower = [];
jouleHeating = [];
datenums = [];
resolution = diff(lat); resolution = resolution(1);

for i = 1:length(tiegcmFiles)
    modelYear = ncread(tiegcmFiles(i).name, 'year')'; modelYear = modelYear(1);
    modelTimes = double(ncread(tiegcmFiles(i).name, 'mtime'));
    modelDatenums = repmat(datenum(num2str(modelYear),'yyyy'), 1, size(modelTimes,2));
    modelDatenums = modelDatenums + modelTimes(1,:) + modelTimes(2,:)/24 + modelTimes(3,:)/1440 - 1;
    datenums = [datenums; modelDatenums'];
    
    cpcp = [cpcp; ncread(tiegcmFiles(i).name, 'ctpoten')];
    %hemisphericPower = [hemisphericPower; ncread(tiegcmFiles(i).name, 'hpower')];
    
    jouleHeatingInteg = double(ncread(tiegcmFiles(i).name, 'QJOULE_INTEG'));
    hpInteg = double(ncread(tiegcmFiles(i).name, 'EFLUX'));
    for j = 1:size(jouleHeatingInteg,3)
        jhThisTime = 0;
        hpThisTime = 0;
        for k = 1:length(lat)
            ringArea = computeRingArea(lat(k), resolution);
            jhThisLat = 1E-12 * ringArea * sum(jouleHeatingInteg(:,k,j)) / length(lon);
            hpThisLat = 1E-12 * ringArea * sum(hpInteg(:,k,j)) / length(lon);
            jhThisTime = jhThisTime + jhThisLat;
            hpThisTime = hpThisTime + hpThisLat;
        end
        jouleHeating = [jouleHeating; jhThisTime];
        hemisphericPower = [hemisphericPower; hpThisTime];
    end
end

amieFile = 'timegcm_amie_dres_v2_apr5_qamie.nc';
amieJoule = [];
if exist(amieFile,'file')
    modelYear = 2010;
    modelTimes = double(ncread(amieFile, 'mtime'));
    modelDatenums = repmat(datenum(num2str(modelYear),'yyyy'), 1, size(modelTimes,2));
    amieDatenums = (modelDatenums + modelTimes(1,:) + modelTimes(2,:)/24 + modelTimes(3,:)/1440 - 1)';
    
    jouleHeatingInteg = double(ncread(amieFile, 'QAMIE'));
    jouleHeatingInteg = reshape(jouleHeatingInteg, size(jouleHeatingInteg,1), ...
                                size(jouleHeatingInteg,2), size(jouleHeatingInteg,4));
    for j = 1:size(jouleHeatingInteg,3)
        jhThisTime = 0;
        for k = 1:length(lat)
            ringArea = computeRingArea(lat(k), resolution);
            jhThisLat = 1E-12 * ringArea * sum(jouleHeatingInteg(:,k,j)) / length(lon);
            jhThisTime = jhThisTime + jhThisLat;
        end
        amieJoule = [amieJoule; jhThisTime];
    end
end

jhInfs = jouleHeating >= 1E36;
jhNoInfs = jouleHeating(~jhInfs);
dateNumsNoInfs = datenums(~jhInfs);
jouleHeating = interp1(dateNumsNoInfs, jhNoInfs, datenums, 'linear', 'extrap');
t1 = datenum('2010-04-05');
t2 = datenum('2010-04-06');

load('goceVariables.mat', 'densityNoBg', 'latitude', 'timestampsDensityDatenum');
goceTimestamps = timestampsDensityDatenum;
goceDens = densityNoBg;

ind = (goceTimestamps >= datenums(1) & goceTimestamps <= datenums(end));
goceDens = goceDens(ind);
goceTimestamps = goceTimestamps(ind);
goceDens = smooth(goceDens, 540);
firstInd = find(goceTimestamps <= t1, 1, 'last');
calmDens = goceDens(firstInd);
gocePlot = goceDens;%((goceDens / calmDens) - 1) * 100;
latitude = latitude(ind);

timestampsEquator = goceTimestamps(latitude >= 10);
times = timestampsEquator(find(diff(timestampsEquator) > 45 * 60/86400) + 1);
orbAverInd = ismember(goceTimestamps, times);
gocePlot = gocePlot(orbAverInd) * 1.23;

figHandle = figure('Color', 'w');
subplot(2,1,1);
hAx(1) = plot(goceTimestamps(orbAverInd), gocePlot, 'linewidth', 2.0);
%ylabel('\Delta \rho [%]');
ylabel('Density')
datetick('x', 'dd/mm');
xlim([t1 t2])
%ylim([0 100])

hold all;
subplot(2,1,2);
hAx(2) = plot(datenums, jouleHeating, 'k');
hold all;
plot(datenums, hemisphericPower, 'r')
datetick('x', 'dd/mm');
xlim([t1 t2])
if ~isempty(amieJoule)
    hAx(2) = plot(amieDatenums, amieJoule, 'k--');
    legend('Weimer JH', 'R&R HP', 'AMIE JH');
else
    legend('Weimer JH', 'R&R HP');
end
ylabel('Power [GW]');
hold off;

% subplot(2,1,2);
% plot(datenums, cpcp)
% ylabel('kV');
% datetick('x', 'dd/mm');
% xlim([t1 t2])
% title('CPCP')

end

function [ringArea] = computeRingArea(lat, resolution)

northLat = lat + resolution / 2;
southLat = lat - resolution / 2;
earthRad = 6371e3;

ringArea = 2 * pi * earthRad.^2 * (sind(northLat) - sind(southLat));

end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function [] = plotSurfs(timeseriesFig, hAx, datenums, tiegcmFileNames)

load('goceVariables.mat', 'latitude', 'magneticLatitude', 'longitude', 'altitude', 'solarTime', ...
    'densityNoBg', 'timestampsDensityDatenum', 'crwindEast', 'crwindNorth');
ind = (timestampsDensityDatenum >= datenums(1) & timestampsDensityDatenum <= datenums(end));
goceDatenums = timestampsDensityDatenum(ind);

goceTimestamps = (goceDatenums - datenums(1)) * 86400;
goceLon = longitude(ind);
goceLat = latitude(ind);
goceAlt = altitude(ind);
goceSolarTime = solarTime(ind);
goceU = crwindEast(ind);
goceV = crwindNorth(ind);

tiegcmTime = (datenums - datenums(1)) * 86400;
tiegcmFiles = dir(tiegcmFileNames);
lat = ncread(tiegcmFiles(1).name, 'lat');
lon = ncread(tiegcmFiles(1).name, 'lon');
lev = ncread(tiegcmFiles(1).name, 'lev');
tiegcmAlt = zeros(length(lon), length(lat), length(lev), length(datenums));
k = 1;
for i = 1:length(tiegcmFiles)
    mtimes = double(ncread(tiegcmFiles(i).name, 'mtime')); numTimes = size(mtimes,2);
    tiegcmAlt(:,:,:,k:k+numTimes-1) = double(ncread(tiegcmFiles(i).name, 'ZG'));
    k = k + numTimes;
end
tiegcmAlt = tiegcmAlt / 100;

%tiegcmAlt(:,:,:,1) = [];
%tiegcmTime = tiegcmTime(2:end);

if ~exist('tiegcmDens.mat', 'file')
    tiegcmDens = zeros(length(lon), length(lat), length(lev), length(datenums));
    k = 1;
    for i = 1:length(tiegcmFiles)
        mtimes = double(ncread(tiegcmFiles(i).name, 'mtime')); numTimes = size(mtimes,2);
        tiegcmDens(:,:,:,k:k+numTimes-1) = double(ncread(tiegcmFiles(i).name, 'DEN'));
        k = k + numTimes;
    end
    %tiegcmDens(:,:,:,1) = [];
    
    tiegcmGoceInterp = nan(size(goceTimestamps));
    tiegcmGoce270km = nan(size(goceTimestamps));
    targetCount = round(length(goceTimestamps) / 500);
    barWidth = 50;
    warning('off', 'MATLAB:qhullmx:InternalWarning');
    p = TimedProgressBar( targetCount, barWidth, ...
                        'Interpolating GOCE, ETA ', ...
                        '. Now at ', ...
                        'Completed in ' );

    for i = 1:length(goceTimestamps)
        tiegcmGoceInterp(i) = interpSatellite(lon, lat, tiegcmAlt, tiegcmTime, tiegcmDens,...
                                                goceLon(i), goceLat(i), goceAlt(i), goceTimestamps(i), 1);
        tiegcmGoce270km(i) = interpSatellite(lon, lat, tiegcmAlt, tiegcmTime, tiegcmDens,...
                                                goceLon(i), goceLat(i), 270E3, goceTimestamps(i), 1);
        if mod(i, 500) == 0
            p.progress;
        end
    end
    p.stop;

    tgNoNans = tiegcmGoceInterp(~isnan(tiegcmGoceInterp));
    tNoNans = goceTimestamps(~isnan(tiegcmGoceInterp));
    tiegcmGoceInterp = interp1(tNoNans, tgNoNans, goceTimestamps, 'nearest', 'extrap');

    tgNoNans = tiegcmGoce270km(~isnan(tiegcmGoce270km));
    tNoNans = goceTimestamps(~isnan(tiegcmGoce270km));
    tiegcmGoce270km = interp1(tNoNans, tgNoNans, goceTimestamps, 'nearest', 'extrap');
    
    if exist('timegcm_dres.s_apr2010_amie_rho.nc', 'file')
        amieDens = ncread('timegcm_dres.s_apr2010_amie_rho.nc','RHO');
        timegcmAlt = ncread('timegcm_dres.s_apr2010_amie_rho.nc','ZG') * 1000;
        modelYear = 2010;
        modelTimes = double(ncread('timegcm_dres.s_apr2010_amie_rho.nc', 'mtime'));
        modelDatenums = repmat(datenum(num2str(modelYear),'yyyy'), 1, size(modelTimes,2));
        modelDatenums = modelDatenums + modelTimes(1,:) + modelTimes(2,:)/24 + modelTimes(3,:)/1440 - 1;
        timegcmTime = (modelDatenums - datenums(1)) * 86400;
        
%         amieDens(:,:,:,end) = [];
%         timegcmAlt(:,:,:,end) = [];
        amieDens(end,:,:,:) = [];
        timegcmAlt(end,:,:,:) = [];
%         zeroCube = zeros(size(timegcmAlt,1), size(timegcmAlt,2), size(timegcmAlt,3));
%         amieDens = cat(4, zeroCube, amieDens);
%         timegcmAlt = cat(4, zeroCube, timegcmAlt);
        
        amieGoce270km = nan(size(goceTimestamps));
        targetCount = round(length(goceTimestamps) / 500);
        barWidth = 50;
        warning('off', 'MATLAB:qhullmx:InternalWarning');
        p = TimedProgressBar( targetCount, barWidth, ...
                            'Interpolating GOCE/AMIE, ETA ', ...
                            '. Now at ', ...
                            'Completed in ' );

        for i = 1:length(goceTimestamps)
            amieGoce270km(i) = interpSatellite(lon, lat, timegcmAlt, timegcmTime, amieDens,...
                                                    goceLon(i), goceLat(i), 270E3, goceTimestamps(i), 1);
            if mod(i, 500) == 0
                p.progress;
            end
        end
        p.stop;
        
        tgNoNans = amieGoce270km(~isnan(amieGoce270km));
        tNoNans = goceTimestamps(~isnan(amieGoce270km));
        amieGoce270km = interp1(tNoNans, tgNoNans, goceTimestamps, 'nearest', 'extrap');
        
        save('tiegcmDens.mat', 'amieGoce270km', '-v7.3')       
    end
    
    tiegcmDatenums = goceDatenums;
    if exist('amieDens', 'var')
        save('tiegcmDens.mat', 'tiegcmGoceInterp', '-append')
    else
        save('tiegcmDens.mat', 'tiegcmGoceInterp', '-v7.3')
    end
    save('tiegcmDens.mat', 'tiegcmGoce270km', '-append')
    save('tiegcmDens.mat', 'tiegcmDatenums', '-append')
end

load tiegcmDens.mat
tiegcmPlot = smooth(tiegcmGoce270km * 1e14, 540);
tiegcmPlot = tiegcmPlot(ismember(tiegcmDatenums, goceDatenums));
t1 = datenum('2010-04-05');
t2 = datenum('2010-04-06');
firstInd = find(goceDatenums <= t1, 1, 'last');
calmDens = tiegcmPlot(firstInd);
%tiegcmPlot = ((tiegcmPlot / calmDens) - 1) * 100;

latitude = latitude(ind);
timestampsEquator = goceDatenums(latitude >= 80);
times = timestampsEquator(find(diff(timestampsEquator) > 45 * 60/86400) + 1);
orbAverInd = ismember(goceDatenums, times);
tiegcmPlot = tiegcmPlot(orbAverInd);

figure(timeseriesFig);
hold all;
subplot(2,1,1)
plot(goceDatenums(orbAverInd), tiegcmPlot, 'linewidth', 2.0)
if exist('amieGoce270km','var')
    timegcmPlot = smooth(amieGoce270km * 1e14, 540);
    timegcmPlot = timegcmPlot(orbAverInd);
    plot(goceDatenums(orbAverInd), timegcmPlot, 'linewidth', 2.0)
    legend('1.23 x GOCE', 'TIEGCM', 'TIMEGCM/AMIE')
else
    legend('1.23 x GOCE', 'TIEGCM')
end

hold off;

tiegcm10minDatenums = datenums;
stormBegin = datenum('2010-04-05 07:00:00');
stormEnd = datenum('2010-04-05 12:00:00');

% plotLonCrossSections(tiegcmFileNames, tiegcm10minDatenums, lon, lat, tiegcmAlt, {'QJOULE','NO_COOL','WN','VN'}, ...
%     goceDatenums, goceSolarTime, goceLat, goceAlt, stormBegin, stormEnd)

plotAltCrossSections(tiegcmFileNames, tiegcm10minDatenums, lon, lat, lev, tiegcmAlt, {'QJOULE_INTEG','AMIE_QJOULE','DEN','AMIE_DEN'}, ...
    goceDatenums, goceLon, goceLat, goceAlt, goceU, goceV, stormBegin, stormEnd, 270)

end

function plotLonCrossSections(tiegcmFileNames, tiegcm10minDatenums, lon, lat, tiegcmAlt, fieldNames, ...
    goceDatenums, goceSolarTime, goceLat, goceAlt, stormBegin, stormEnd)

fields = cell(length(fieldNames),1);
tiegcmFiles = dir(tiegcmFileNames);
for i = 1:length(fieldNames)
    if ~isempty(strfind(fieldNames{i},'AMIE'))
        if exist('timegcm_amie_dres_v2_apr5_qamie.nc','file') && ...
            exist('timegcm_dres.s_apr2010_amie_rho.nc','file')
            if ~isempty(strfind(fieldNames{i},'QJOULE'))
                temporaryField = double(ncread('timegcm_amie_dres_v2_apr5_qamie.nc', 'QAMIE'));
                amieAlt = double(ncread('timegcm_amie_dres_v2_apr5_qamie.nc', 'Z'))/100;
            elseif ~isempty(strfind(fieldNames{i},'DEN'))
                timegcmAlt = double(ncread('timegcm_dres.s_apr2010_amie_rho.nc', 'Z'))*1000;
                temporaryField = double(ncread('timegcm_dres.s_apr2010_amie_rho.nc', 'RHO'));
                temporaryField(:,:,:,end) = [];
            end
            if size(temporaryField,4) ~= size(tiegcmAlt,4)
                error('TIEGCM and AMIE results must contain the same number of time steps!')
            end
        else
            error('AMIE Files do not exist!')
        end
    else
        temporaryField = zeros(size(tiegcmAlt));
        k = 1;
        for j = 1:length(tiegcmFiles)
            mtimes = double(ncread(tiegcmFiles(j).name, 'mtime')); numTimes = size(mtimes,2);
            temporaryField(:,:,:,k:k+numTimes-1) = double(ncread(tiegcmFiles(j).name, fieldNames{i}));
            k = k + numTimes;
        end
    end

    if any(temporaryField(:,:,end,:) >= 1e35)
        temporaryField(:,:,end,:) = temporaryField(:,:,end-1,:);
    end
    if ~isempty(strfind(fieldNames{i}, 'DEN')) || ~isempty(strfind(fieldNames{i}, 'O_N2'))
        temporaryField = log(temporaryField);
    end
    if ~isempty(strfind(fieldNames{i}, 'VN')) || ~isempty(strfind(fieldNames{i}, 'WN'))
        temporaryField = temporaryField / 100;
    end
    
    fields{i} = temporaryField;
end

numPlotRows = ceil(length(fieldNames) / 2);
[latGrid,~] = meshgrid(lat, 1:size(tiegcmAlt,3));
plotLocalTime = 7; % hours
figHandle = figure('units','normalized','outerposition',[0 0 1 1]);

[~,beginInd] = min(abs(tiegcm10minDatenums - stormBegin));
[~,endInd] = min(abs(tiegcm10minDatenums - stormEnd));

axisLim = repmat([realmax, -realmax], length(fieldNames), 1);

for t = beginInd:endInd
    [~,~,~,hours,minutes,seconds] = datevec(tiegcm10minDatenums(t));
    hourNow = hours + minutes/60 + seconds/3600;
    plotLon = 15 * (plotLocalTime - hourNow);
    if plotLon < -180; plotLon = plotLon + 360; end
    if plotLon >= 180; plotLon = plotLon - 360; end
    [~,plotLonInd(t)] = min(abs(lon - plotLon));
    
    for i = 1:length(fieldNames)
        fieldSlice = fields{i}(plotLonInd(t),:,:,t);
        if min(fieldSlice(:)) < axisLim(i,1)
            axisLim(i,1) = min(fieldSlice(:));
        end
        if max(fieldSlice(:)) > axisLim(i,2)
            axisLim(i,2) = max(fieldSlice(:));
        end
    end
end

writerObj = VideoWriter('LatCrossSect.avi');
writerObj.FrameRate = 1;
open(writerObj);

for t = beginInd:endInd
    for i = 1:length(fieldNames)
        
        if ~isempty(strfind(fieldNames{i},'AMIE_QJOULE'))
            alt = amieAlt;
        elseif ~isempty(strfind(fieldNames{i},'AMIE_DEN'))
            alt = timegcmAlt;
        else
            alt = tiegcmAlt;
        end
        altSlice = alt(plotLonInd(t),:,:,t);
        altSlice = reshape(altSlice,size(alt,2),size(alt,3));
        altSlice = flipud(altSlice') / 1E3;
        latForGeopHeight = repmat(lat, 1, length(altSlice(1:4:end,1)));
        geopHeight = altSlice(1:4:end,:)';
        
        subplot(numPlotRows,2,i)
        fieldSlice = fields{i}(plotLonInd(t),:,:,t);
        fieldSlice = reshape(fieldSlice,size(alt,2),size(alt,3));
        fieldSlice = flipud(fieldSlice');
        surf(latGrid, altSlice, fieldSlice, 'linestyle', 'none', 'edgecolor', 'none')
        view(2);
        grid off;
        shading interp
        xlim([min(lat) max(lat)])
        ylim([100 600])
        colorbar;
        caxis(axisLim(i,:))
        xlabel('Geographic Latitude', 'fontsize', 14)
        ylabel('Geometric Height [km]', 'fontsize', 14)
        fieldLongName = ncreadatt(tiegcmFiles(1).name, fieldNames{i}, 'long_name');
        fieldUnit = ncreadatt(tiegcmFiles(1).name, fieldNames{i}, 'units');
        if strcmpi(fieldNames{i}, 'VN') || strcmpi(fieldNames{i}, 'WN')
            fieldUnit = 'm/s';
        end
        title([fieldLongName, ' [', fieldUnit, ']'], 'fontsize', 14)
        set(gca,'FontSize',14)
        
        hold all;
        [~,nearestGoce] = min(abs(goceDatenums - tiegcm10minDatenums(t)));
        if goceSolarTime(nearestGoce) > plotLocalTime-1 && goceSolarTime(nearestGoce) < plotLocalTime+1
            plot(goceLat(nearestGoce), goceAlt(nearestGoce)/1E3, 'kx');
        end
        
        geopPlotHeight = repmat(max(fieldSlice(:)),size(geopHeight,1),size(geopHeight,2));
        plot3(latForGeopHeight, geopHeight, geopPlotHeight, 'color', 'w', 'linewidth', 1.5)
        hold off;
    end
    
    supertitle = ['UT ', datestr(tiegcm10minDatenums(t),'HH:MM')];
    [~,h] = suplabel(supertitle  ,'t');
    set(h,'FontSize',16)
    
    frame = getframe(figHandle);
    writeVideo(writerObj,frame);
end
close(writerObj);

end

function plotAltCrossSections(tiegcmFileNames, tiegcm10minDatenums, lon, lat, lev, tiegcmAlt, fieldNames, ...
    goceDatenums, goceLon, goceLat, goceAlt, goceU, goceV, stormBegin, stormEnd, plotAlt)

fields = cell(length(fieldNames),1);
tiegcmFiles = dir(tiegcmFileNames);

windPlot = 'DEN';
wnFound = 0;
for i = 1:length(fieldNames)
    if strcmpi(fieldNames{i},'WN')
        wnFound = 1;
        break;
    end
end
if wnFound == 1
    windPlot = 'WN';
end

newFieldNames = fieldNames;
plotFieldNames = fieldNames;
for i = 1:length(fieldNames)
    if strcmpi(fieldNames{i},'QJOULE_INTEG')
        newFieldNames = [newFieldNames, {'UI_ExB', 'VI_ExB'}];
        ionDriftInd = length(newFieldNames) - 1;
    elseif strcmpi(fieldNames{i},windPlot)
        newFieldNames = [newFieldNames, {'UN', 'VN'}];
        windInd = length(newFieldNames) - 1;
    end
end
fieldNames = newFieldNames;

for i = 1:length(fieldNames)
    
    if ~isempty(strfind(fieldNames{i},'AMIE'))
        if exist('timegcm_amie_dres_v2_apr5_qamie.nc','file') && ...
            exist('timegcm_dres.s_apr2010_amie_rho.nc','file')
            if ~isempty(strfind(fieldNames{i},'QJOULE'))
                temporaryField = double(ncread('timegcm_amie_dres_v2_apr5_qamie.nc', 'QAMIE'));
                temporaryField = reshape(temporaryField, size(temporaryField,1), size(temporaryField,2),...
                                            size(temporaryField,4));
            elseif ~isempty(strfind(fieldNames{i},'DEN'))
                timegcmAlt = double(ncread('timegcm_dres.s_apr2010_amie_rho.nc', 'ZG'))*1000;
                temporaryField = double(ncread('timegcm_dres.s_apr2010_amie_rho.nc', 'RHO'));
                temporaryField(:,:,:,end) = [];
                timegcmAlt(:,:,:,end) = [];
                temporaryField(end,:,:,:) = [];
                timegcmAlt(end,:,:,:) = [];
                zeroCube = zeros(size(timegcmAlt,1), size(timegcmAlt,2), size(timegcmAlt,3));
                temporaryField = cat(4, zeroCube, temporaryField);
                timegcmAlt = cat(4, zeroCube, timegcmAlt);
                if size(temporaryField,4) ~= size(tiegcmAlt,4)
                    error('TIEGCM and AMIE results must contain the same number of time steps!')
                end
            end
        else
            error('AMIE Files do not exist!')
        end
    else
        temporaryField = double(ncread(tiegcmFiles(1).name, fieldNames{i}));
        if length(size(temporaryField)) == 3
            temporaryField = zeros(size(tiegcmAlt,1),size(tiegcmAlt,2),size(tiegcmAlt,4));
        else
            temporaryField = zeros(size(tiegcmAlt));
        end
        k = 1;
        for j = 1:length(tiegcmFiles)
            mtimes = double(ncread(tiegcmFiles(j).name, 'mtime')); numTimes = size(mtimes,2);
            if length(size(temporaryField)) == 3
                temporaryField(:,:,k:k+numTimes-1) = double(ncread(tiegcmFiles(j).name, fieldNames{i}));
            else
                temporaryField(:,:,:,k:k+numTimes-1) = double(ncread(tiegcmFiles(j).name, fieldNames{i}));
            end
            k = k + numTimes;
        end
    end
%     if length(size(temporaryField)) == 4 && any(temporaryField(:,:,end,:) >= 1e35)
%         temporaryField(:,:,end,:) = temporaryField(:,:,end-1,:);
%     end
    
    fields{i} = temporaryField;
end

endLat = 40;
if endLat > 0
    latInd = lat > endLat;
else
    latInd = lat < endLat;
end
latToPlot = lat(latInd);

[lonGrid,latGrid] = meshgrid(lon, latToPlot);
[x, y, ~] = geod2ecef(latGrid(:), lonGrid(:), plotAlt*1E3*ones(numel(latGrid),1));
x = reshape(x, size(latGrid,1), size(latGrid,2));
y = reshape(y, size(latGrid,1), size(latGrid,2));
radius = x.^2 + y.^2;
plotRadius = mean(radius, 2);

axisLim = zeros(length(fieldNames),2);
for i = 1:length(fieldNames)
    if length(size(fields{i})) == 4
        interpolatedField = zeros(length(lon), length(latToPlot), length(tiegcm10minDatenums));
        fourDimField = fields{i};
        
        
        targetCount = length(tiegcm10minDatenums);
        barWidth = 50;
        progressString = ['Interpolating ', fieldNames{i}, ', ETA '];
        p = TimedProgressBar( targetCount, barWidth, ...
                        progressString, ...
                        '. Now at ', ...
                        'Completed in ' );
        if ~isempty(strfind(fieldNames{i},'AMIE_DEN'))
            alt = timegcmAlt;
        else
            alt = tiegcmAlt;
        end
        nlev = size(alt,3);
        
        parfor t = 1:length(tiegcm10minDatenums)
            tempField = zeros(length(lon), length(latToPlot));
            for latIter = 1:size(latGrid,1)
                for lonIter = 1:size(lonGrid,2)
                    altCol = alt(lonIter,latIter,:,t);
                    levelBelow = find(altCol < plotAlt*1E3 & altCol > 0, 1, 'last');
                    if isempty(levelBelow) || levelBelow == nlev
                        tempField(lonIter, latIter) = nan(1);
                    else
                        zInterp = [altCol(levelBelow), altCol(levelBelow+1)];
                        valInterp = [fourDimField(lonIter, latIter, levelBelow, t),...
                                     fourDimField(lonIter, latIter, levelBelow+1, t)];
                        tempField(lonIter, latIter) = interp1(zInterp, valInterp, plotAlt*1E3);
                    end
                end             
            end
            interpolatedField(:,:,t) = tempField;
            p.progress;
        end
        p.stop;
        interpolatedField(isnan(interpolatedField)) = 0;
        fields{i} = interpolatedField;
    else
        fields{i} = fields{i}(:,latInd,:);
    end
    
    if i > length(plotFieldNames)
        scaleSpeed = 500; % m/s
        scaleFac = 10;
        fields{i} = scaleFac * (fields{i} / (scaleSpeed*100));
    end
    
    if strcmpi(fieldNames{i}, 'WN')
        fields{i} = fields{i} / 100;
    end
%     if strcmpi(fieldNames{i}, 'DEN')
%         fields{i} = (fields{i} ./ repmat(fields{i}(:,:,1),[1,1,length(tiegcm10minDatenums)]) - 1) * 100;
%     end
    
    axisLim(i,:) = [min(fields{i}(:)), max(fields{i}(:))];
    if strcmpi(fieldNames{i}, 'QJOULE_INTEG')
        weimerJouleInd = i;
    elseif strcmpi(fieldNames{i}, 'DEN')
        denInd = i;
    end
    %fields{i}(1,1:10,1) = ones(1,10)*axisLim(i,2);
    
end

[~,beginInd] = min(abs(tiegcm10minDatenums - stormBegin));
[~,endInd] = min(abs(tiegcm10minDatenums - stormEnd));
numPlotRows = ceil(length(plotFieldNames)/2);

if length(fieldNames) == 1
    numPlotCols = 1;
else
    numPlotCols = 2;
end

figHandle = figure('units','normalized','outerposition',[0 0 1 1]);

writerObj = VideoWriter('AltCrossSect.avi');
writerObj.FrameRate = 1;
open(writerObj);

for t = beginInd:endInd
    [~,~,~,hours,minutes,seconds] = datevec(tiegcm10minDatenums(t));
    hourNow = hours + minutes/60 + seconds/3600;
    lst = hourNow + lon/15;
    lst(lst < 0) = lst(lst < 0) + 24;
    lst(lst >= 24) = lst(lst >= 24) - 24;
    plotAngle = (lst*2*pi/24)';
    midnightInd = find(abs(diff(lst)) > 12);
    if ~isempty(midnightInd)
        plotAngle = [plotAngle(midnightInd+1:end), plotAngle(1:midnightInd)];
    end
    
    
    for i = 1:length(plotFieldNames)
        hAx = subplot(numPlotRows,numPlotCols,i);
        plotField = rot90(fields{i}(:,:,t));
        
        opengl neverselect
        if ~isempty(midnightInd)
            plotField = cat(2, plotField(:,midnightInd+1:end), plotField(:,1:midnightInd));           
        end
        plotDist = 90-latToPlot;
        polarplot3d(plotField, 'RadialRange', plotDist, 'TickSpacing', 45, 'InterpMethod', 'spline', 'MeshScale', [0.5 1.0],...
            'PolarGrid', {1 1});
        view(2);
        shading interp;
        %grid off;
        axis off;
        if numPlotRows == 1; cbarLoc = 'northoutside'; else cbarLoc = 'eastoutside'; end
        colorbar(cbarLoc)
        if strcmpi(fieldNames{i}, 'AMIE_QJOULE')
            caxis(axisLim(weimerJouleInd,:))
        elseif strcmpi(fieldNames{i}, 'AMIE_DEN')
            caxis(axisLim(denInd,:))
        else
            caxis(axisLim(i,:))
        end
        
        if isempty(strfind(fieldNames{i},'AMIE'))
            fieldLongName = ncreadatt(tiegcmFiles(1).name, fieldNames{i}, 'long_name');
            fieldUnit = ncreadatt(tiegcmFiles(1).name, fieldNames{i}, 'units');
        else
            if strcmpi(fieldNames{i}, 'AMIE_QJOULE')
                fieldLongName = 'AMIE Joule Heating';
                fieldUnit = 'mW/m^2';
            elseif strcmpi(fieldNames{i}, 'AMIE_DEN')
                fieldLongName = 'TIMEGCM/AMIE density';
                fieldUnit = 'g/cm^3';
            end
        end
        
        if strcmpi(fieldNames{i}, 'WN')
            fieldUnit = 'm/s';
            velocityFieldName = ' and Neutral Wind';
        elseif strcmpi(fieldNames{i}, windPlot)
            velocityFieldName = ' and Neutral Wind';
        elseif strcmpi(fieldNames{i}, 'QJOULE_INTEG')
            velocityFieldName = ' and E x B';
        else
            velocityFieldName = [];
        end
        title([fieldLongName, ' [', fieldUnit, ']', velocityFieldName], 'fontsize', 14)
        set(gca,'FontSize',14)
        
        if strcmpi(fieldNames{i}, 'QJOULE_INTEG')
            uMat = rot90(fields{ionDriftInd}(:,:,t));
            vMat = rot90(fields{ionDriftInd+1}(:,:,t));
%             uMat = zeros(size(uMat));
%             uMat(9,1:100) = ones(1, length(uMat(9,1:100))) * 5;
%             vMat = uMat;
        elseif strcmpi(fieldNames{i}, windPlot)
            uMat = rot90(fields{windInd}(:,:,t));
            vMat = rot90(fields{windInd+1}(:,:,t));
%             vMat = zeros(size(uMat));
%             uMat = ones(size(vMat)) * 6;
        end
        if strcmpi(fieldNames{i}, 'QJOULE_INTEG') || strcmpi(fieldNames{i}, windPlot)
            if ~isempty(midnightInd)
                uMat = cat(2, uMat(:,midnightInd+1:end), uMat(:,1:midnightInd));  
                vMat = cat(2, vMat(:,midnightInd+1:end), vMat(:,1:midnightInd)); 
            end
            lonWind = 1:4:size(uMat,2);
            latWind = 1:2:size(uMat,1);
            uMat = uMat(latWind,lonWind); vMat = vMat(latWind,lonWind);
            vertVel = zeros(size(uMat));
            
            plotHeight = ones(size(uMat)) * axisLim(i,2);
            [theta, rho] = meshgrid(plotAngle(lonWind), plotDist(latWind)');
            rho = flipud(rho);
            [plotX, plotY] = pol2cart(theta, rho);
            uRotated = uMat.*cos(theta) - vMat.*sin(theta);
            vRotated = uMat.*sin(theta) + vMat.*cos(theta);
            rotatedWind = [uRotated(:), vRotated(:), zeros(size(uRotated(:)))];
            zDir = repmat([0,0,1], size(rotatedWind,1),1);
            finalWind = cross(zDir, rotatedWind);
            uFinal = reshape(finalWind(:,1), size(uMat)); vFinal = reshape(finalWind(:,2), size(vMat));
            
            hold all;
            
            quiver3(hAx, plotX, plotY, plotHeight, uFinal, vFinal, vertVel, 0, 'Color', 'r', 'lineWidth', 1.0);
            if strcmpi(fieldNames{i}, windPlot)
                fiveMin = 5/1440;
                td = tiegcm10minDatenums;
                goceInd = (goceDatenums >= td(t)-fiveMin & goceDatenums < td(t)+fiveMin) & goceLat > min(latToPlot);
                
                if any(goceInd)
                    nearestLats = goceLat(goceInd);
                    nearestLons = goceLon(goceInd);
                    nearestU = scaleFac * goceU(goceInd) / scaleSpeed;
                    nearestV = scaleFac * goceV(goceInd) / scaleSpeed;
                    consInd = [1:3:length(nearestLats)];
                    nearestLats = nearestLats(consInd);
                    nearestLons = nearestLons(consInd);
                    nearestU = nearestU(consInd);
                    nearestV = nearestV(consInd);
                    
                    goceLst = hourNow + nearestLons/15;
                    goceLst(goceLst < 0) = goceLst(goceLst < 0) + 24;
                    goceLst(goceLst >= 24) = goceLst(goceLst >= 24) - 24;
                    theta = goceLst*2*pi/24;
                    goceDist = 90 - nearestLats;
                    
                    [goceX, goceY] = pol2cart(theta, goceDist);
                    uRotated = nearestU.*cos(theta) - nearestV.*sin(theta);
                    vRotated = nearestU.*sin(theta) + nearestV.*cos(theta);
                    
                    rotatedWind = [uRotated(:), vRotated(:), zeros(size(uRotated(:)))];
                    zDir = repmat([0,0,1], size(rotatedWind,1),1);
                    finalWind = cross(zDir, rotatedWind);
                    uFinal = reshape(finalWind(:,1), size(nearestU)); vFinal = reshape(finalWind(:,2), size(nearestV));
                    vertVel = zeros(length(uFinal),1);
                    plotHeight = ones(size(goceX)) * axisLim(i,2);
                    
                    quiver3(hAx, goceX, goceY, plotHeight, uFinal, vFinal, vertVel, 0, 'Color', 'k', 'lineWidth', 1.5);                  
                end
                quiverscale(10, 0, hAx, scaleFac, scaleSpeed);
            end
            hold off;
        end
    end
    
    supertitle = ['UT ', datestr(tiegcm10minDatenums(t),'HH:MM')];
    [~,h] = suplabel(supertitle  ,'t');
    set(h,'FontSize',16)
    
    frame = getframe(figHandle);
    writeVideo(writerObj,frame);
end
close(writerObj);

end

function quiverscale(px,py,H,scaleFac,scaleSpeed)
%QUIVERSCALE creates a scale at the bottom of the quiver plot, 
%which plots the longest vector of the data set(px,py) and displays 
%the vector's corresponding magnitude. The H input argument is the 
%handle to the quiver plot axes.
%
%Please note that the QUIVERSCALE function returns 
%accurate results only when the vectors are plotted 
%to actual scale (e.g.  QUIVER(U,V,S), where S=0).
%
%Vincent Hodges 02/03/98

h=get(H,'position');
h2=h;
h1=h(2)/4;
h(2)=h(2)+h1;
h(4)=h(4)-h1;
set(gca,'position',h);
hlim=get(gca,'xlim');
h2(end)=h2(end)/20;
h2(2)=h2(2)/3;
ax=axes('position',h2);
hlim(1)=0;
set(ax,'xlim',hlim);

xsq=px.^2;
ysq=py.^2;
arrowmag=max(unique(sqrt([xsq+ysq])));
size1=prod(size(arrowmag));
arrowmag2=reshape(arrowmag,size1,1);
X=zeros(size1,1);
Y=1:size1;
null=zeros(size1,1);
quiver(X,Y,arrowmag2,null,0);
set(ax,'xlim',hlim);
xTicks = [0 scaleFac];
xTickLabels = {'',num2str(scaleSpeed)};
set(ax,'ytick',[]);
set(ax,'xtick',xTicks);
set(ax,'xtickLabel',xTickLabels);
figColor = get(gcf,'color');
set(ax,'XColor',figColor,'YColor',figColor)
set(ax,'color', 'none')


end

function [valAtSatLoc] = interpSatellite(lon, lat, altGrid, modelTime, modelField, satLon, satLat, satAlt, satTime, interpOption)

if satTime < modelTime(1) || satTime > modelTime(end)
    valAtSatLoc = nan(1);
    return
end

i = find(lon <= satLon, 1, 'last');
j = find(lat <= satLat, 1, 'last');
if isempty(i)
    i = 1;
end
if isempty(j)
    j = 1;
end

if interpOption == 1
    if satLat < min(lat)
        jNext = 0;
    else jNext = j + 1; 
    end
    if i == length(lon)
        iNext = 1;
    elseif satLon < min(lon)
        iNext = length(lon);
    else iNext = i + 1;
    end

    if jNext == 0
        interpLon = [lon; lon; lon; lon];
        interpLat = ones(size(interpLon)) * min(lat);
        lonInd = (1:length(lon))'; latInd = ones(size(lonInd))*1;
    elseif jNext == length(lat) + 1
        interpLon = [lon; lon; lon; lon];
        interpLat = ones(size(interpLon)) * max(lat);
        lonInd = (1:length(lon))'; latInd = ones(size(lonInd))*length(lat);
    else
        interpLon = repmat([lon(i); lon(i); lon(iNext); lon(iNext)], 4, 1);
        interpLat = repmat([lat(j); lat(jNext); lat(jNext); lat(j)], 4, 1);
        lonInd = [i;i;iNext;iNext];
        latInd = [j;jNext;jNext;j];
    end
end

nearestTime = find(modelTime <= satTime, 1, 'last');
if nearestTime == length(modelTime)
    nextTime = nearestTime;
else nextTime = nearestTime + 1;
end

vals = zeros(2,1);
for n = 1:2
    if n == 1
        t = nearestTime;
    else t = nextTime;
    end
    
    if satAlt < altGrid(i,j,1,t) || satAlt > altGrid(i,j,end-1,t)
        valAtSatLoc = nan(1);
        return
    end
    
    altCol = altGrid(i,j,:,t);
    k = find(altCol(1:end-1) <= satAlt, 1, 'last');
    
    if interpOption == 1
        tInd = ones(size(lonInd)) * t;
        nlev = size(altGrid, 3);
        interpAlt = [];
        interpField = [];
        for m = -1:2
            if k + m < 1 || k + m > nlev-1
                k_Ind = ones(size(lonInd)) * k;
                km_Ind = ones(size(lonInd)) * (k-m);
                linInd1 = sub2ind(size(altGrid), lonInd, latInd, k_Ind, tInd);
                linInd2 = sub2ind(size(altGrid), lonInd, latInd, km_Ind, tInd);
                extrapAlt = 2 * altGrid(linInd1) - altGrid(linInd2);
                extrapField = modelField(linInd1).^2 ./ modelField(linInd2);
                interpAlt = [interpAlt; extrapAlt];
                interpField = [interpField; extrapField];
            else
                km_Ind = ones(size(lonInd)) * (k+m);
                linInd = sub2ind(size(altGrid), lonInd, latInd, km_Ind, tInd);
                interpAlt = [interpAlt; altGrid(linInd)];
                interpField = [interpField; modelField(linInd)];
            end
        end

        tiegcmCoord = geod2ecef(interpLat, interpLon, interpAlt) / 1E6;
        satCoord = geod2ecef(satLat, satLon, satAlt) / 1E6;
        vals(n) = griddatan(double(tiegcmCoord), double(interpField), satCoord, 'linear');
    else
        vals(n) = modelField(i,j,k,t);
    end
end

tInterp = modelTime([nearestTime; nextTime]);
if nearestTime ~= nextTime
    valAtSatLoc = interp1(tInterp, vals, satTime);
else
    valAtSatLoc = vals(1);
end

end

function [x, y, z] = geod2ecef(latitude, longitude, altitude)

% GEOD2ECEF Convert geodetic coordinates to ECEF coordinates.
% 
% Usage: [X, Y, Z] = GEOD2ECEF(LATITUDE, LONGITUDE, ALTITUDE)
%     or [X, Y, Z] = GEOD2ECEF(LLA)
%     or XYZ = GEOD2ECEF(LATITUDE, LONGITUDE, ALTITUDE)
%     or XYZ = GEOD2ECEF(LLA)
% 
% Converts geodetic coordinates LATITUDE, LONGITUDE, and ALTITUDE to
% Earth-centered, Earth fixed (ECEF) coordinates X, Y, and Z. The inputs
% can either be three separate arguments or 1 matrix. For a matrix input,
% the first dimension with length 3 is assumed to have the three separate
% LATITUDE, LONGITUDE, and ALTITUDE inputs across it. The World Geodetic
% System 1984 (WGS84) ellipsoid model of the Earth is assumed.
% 
% Inputs:
%   -LATITUDE: Geodetic latitude in degrees.
%   -LONGITUDE: Geodetic longitude in degrees.
%   -ALTITUDE: Height above the Earth in meters.
%   -LLA: Matrix with at least one dimension with length 3, the first of
%   which corresponding to the dimension across which the three inputs
%   above go.
% 
% Ouputs:
%   -X: x coordinates of the point in meters.
%   -Y: y coordinates of the point in meters.
%   -Z: z coordinates of the point in meters.
%   -XYZ: When just one output is requested, the three outputs above are
%   returned as a row vector for scalar inputs, an M-by-3 matrix for column
%   vector inputs, a 3-by-M matrix for row vector inputs, or the three
%   outputs concatenated either along the next largest dimension when the
%   inputs are separate arguments or the same dimension that the inputs
%   went across when a single matrix is input.
% 
% Copyright (c) 2011, Drew Compston
% All rights reserved.
 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


% Input checking/conversion.
narginchk(1, 3);
if nargin == 1
    sizelatitude = size(latitude);
    first3 = find(sizelatitude == 3, 1, 'first');
    latitude = reshape(permute(latitude, [first3, 1:(first3 - 1), ...
        (first3 + 1):ndims(latitude)]), 3, []);
    sizelatitude(first3) = 1;
    longitude = reshape(latitude(2, :), sizelatitude);
    altitude = reshape(latitude(3, :), sizelatitude);
    latitude = reshape(latitude(1, :), sizelatitude);
end
latitude = latitude*pi/180; longitude = longitude*pi/180;

% WGS84 parameters.
a = 6378137; f = 1/298.257223563; b = a*(1 - f); e2 = 1 - (b/a)^2;

% Conversion from:
% en.wikipedia.org/wiki/Geodetic_system#Conversion_calculations
Nphi = a ./ sqrt(1 - e2*sin(latitude).^2);
x = (Nphi + altitude).*cos(latitude).*cos(longitude);
y = (Nphi + altitude).*cos(latitude).*sin(longitude);
z = (Nphi.*(1 - e2) + altitude).*sin(latitude);

% Shape output according to number of arguments.
if nargout <= 1
    if nargin == 1
        x = cat(first3, x, y, z);
    else
        dims = ndims(x);
        if dims == 2
            if size(x, 2) == 1
                x = [x, y, z];
            else
                x = [x; y; x];
            end
        else
            x = cat(dims + 1, x, y, z);
        end
    end
end

end

function [ax,h]=suplabel(text,whichLabel,supAxes)
% PLaces text as a title, xlabel, or ylabel on a group of subplots.
% Returns a handle to the label and a handle to the axis.
%  [ax,h]=suplabel(text,whichLabel,supAxes)
% returns handles to both the axis and the label.
%  ax=suplabel(text,whichLabel,supAxes)
% returns a handle to the axis only.
%  suplabel(text) with one input argument assumes whichLabel='x'
%
% whichLabel is any of 'x', 'y', 'yy', or 't', specifying whether the 
% text is to be the xlable, ylabel, right side y-label, 
% or title respectively.
%
% supAxes is an optional argument specifying the Position of the 
%  "super" axes surrounding the subplots. 
%  supAxes defaults to [.08 .08 .84 .84]
%  specify supAxes if labels get chopped or overlay subplots
%
% EXAMPLE:
%  subplot(2,2,1);ylabel('ylabel1');title('title1')
%  subplot(2,2,2);ylabel('ylabel2');title('title2')
%  subplot(2,2,3);ylabel('ylabel3');xlabel('xlabel3')
%  subplot(2,2,4);ylabel('ylabel4');xlabel('xlabel4')
%  [ax1,h1]=suplabel('super X label');
%  [ax2,h2]=suplabel('super Y label','y');
%  [ax3,h2]=suplabel('super Y label (right)','yy');
%  [ax4,h3]=suplabel('super Title'  ,'t');
%  set(h3,'FontSize',30)
%
% SEE ALSO: text, title, xlabel, ylabel, zlabel, subplot,
%           suptitle (Matlab Central)

% Author: Ben Barrowes <barrowes@alum.mit.edu>

%modified 3/16/2010 by IJW to make axis behavior re "zoom" on exit same as
%at beginning. Requires adding tag to the invisible axes


currax=findobj(gcf,'type','axes','-not','tag','suplabel');

if nargin < 3
 supAxes=[.08 .08 .84 .84];
 ah=findall(gcf,'type','axes');
 if ~isempty(ah)
  supAxes=[inf,inf,0,0];
  leftMin=inf;  bottomMin=inf;  leftMax=0;  bottomMax=0;
  axBuf=.04;
  set(ah,'units','normalized')
  ah=findall(gcf,'type','axes');
  for ii=1:length(ah)
   if strcmp(get(ah(ii),'Visible'),'on')
    thisPos=get(ah(ii),'Position');
    leftMin=min(leftMin,thisPos(1));
    bottomMin=min(bottomMin,thisPos(2));
    leftMax=max(leftMax,thisPos(1)+thisPos(3));
    bottomMax=max(bottomMax,thisPos(2)+thisPos(4));
   end
  end
  supAxes=[leftMin-axBuf,bottomMin-axBuf,leftMax-leftMin+axBuf*2,bottomMax-bottomMin+axBuf*2];
 end
end
if nargin < 2, whichLabel = 'x';  end
if nargin < 1, help(mfilename); return; end

if ~isstr(text) | ~isstr(whichLabel)
  error('text and whichLabel must be strings')
end
whichLabel=lower(whichLabel);

ax=axes('Units','Normal','Position',supAxes,'Visible','off','tag','suplabel');
if strcmp('t',whichLabel)
  set(get(ax,'Title'),'Visible','on')
  title(text);
elseif strcmp('x',whichLabel)
  set(get(ax,'XLabel'),'Visible','on')
  xlabel(text);
elseif strcmp('y',whichLabel)
  set(get(ax,'YLabel'),'Visible','on')
  ylabel(text);
elseif strcmp('yy',whichLabel)
  set(get(ax,'YLabel'),'Visible','on')
  ylabel(text);
  set(ax,'YAxisLocation','right')
end

for k=1:length(currax), axes(currax(k));end % restore all other axes

if (nargout < 2)
  return
end
if strcmp('t',whichLabel)
  h=get(ax,'Title');
  set(h,'VerticalAlignment','middle')
elseif strcmp('x',whichLabel)
  h=get(ax,'XLabel');
elseif strcmp('y',whichLabel) | strcmp('yy',whichLabel)
  h=get(ax,'YLabel');
end

end