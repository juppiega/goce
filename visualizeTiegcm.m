function [] = visualizeTiegcm(tiegcmFiles)

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

jhInfs = jouleHeating >= 1E36;
jhNoInfs = jouleHeating(~jhInfs);
dateNumsNoInfs = datenums(~jhInfs);
jouleHeating = interp1(dateNumsNoInfs, jhNoInfs, datenums, 'linear', 'extrap');
t1 = datenum('2010-04-04');
t2 = datenum('2010-04-09');

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
legend('Joule Heating', 'Energy Flux');
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

load('goceVariables.mat', 'latitude', 'magneticLatitude', 'longitude', 'altitude', 'solarTime', 'densityNoBg', 'timestampsDensityDatenum');
ind = (timestampsDensityDatenum >= datenums(1) & timestampsDensityDatenum <= datenums(end));
goceDatenums = timestampsDensityDatenum(ind);

if ~exist('tiegcmDens.mat', 'file')
    goceTimestamps = (goceDatenums - goceDatenums(1)) * 86400;
    goceLon = longitude(ind);
    goceLat = latitude(ind);
    goceAlt = altitude(ind);

    tiegcmTime = (datenums - datenums(1)) * 86400;
    tiegcmFiles = dir(tiegcmFileNames);
    lat = ncread(tiegcmFiles(1).name, 'lat');
    lon = ncread(tiegcmFiles(1).name, 'lon');
    lev = ncread(tiegcmFiles(1).name, 'lev');
    tiegcmAlt = zeros(length(lon), length(lat), length(lev), length(datenums));
    tiegcmDens = zeros(length(lon), length(lat), length(lev), length(datenums));
    k = 1;
    for i = 1:length(tiegcmFiles)
        mtimes = double(ncread(tiegcmFiles(i).name, 'mtime')); numTimes = size(mtimes,2);
        tiegcmAlt(:,:,:,k:k+numTimes-1) = double(ncread(tiegcmFiles(i).name, 'ZG'));
        tiegcmDens(:,:,:,k:k+numTimes-1) = double(ncread(tiegcmFiles(i).name, 'DEN'));
        k = k + numTimes;
    end
    tiegcmAlt = tiegcmAlt / 100;

    tiegcmAlt(:,:,:,1) = [];
    tiegcmDens(:,:,:,1) = [];
    tiegcmTime = tiegcmTime(2:end);

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
                                                goceLon(i), goceLat(i), goceAlt(i), goceTimestamps(i));
        tiegcmGoce270km(i) = interpSatellite(lon, lat, tiegcmAlt, tiegcmTime, tiegcmDens,...
                                                goceLon(i), goceLat(i), 270E3, goceTimestamps(i));
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
    
    tiegcmDatenums = goceDatenums;

    save('tiegcmDens.mat', 'tiegcmGoceInterp', '-v7.3')
    save('tiegcmDens.mat', 'tiegcmGoce270km', '-append')
    save('tiegcmDens.mat', 'tiegcmDatenums', '-append')
end

load tiegcmDens.mat
tiegcmPlot = smooth(tiegcmGoce270km * 1e14, 540);
t1 = datenum('2010-04-04');
t2 = datenum('2010-04-09');
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
%datetick('x', 'dd/mm');
%xlim([t1 t2])
legend('1.23 x GOCE', 'TIEGCM')

hold off;

goceDens270km = densityNoBg(ind);
magLat = magneticLatitude(ind);
solarTime = solarTime(ind);
morningIndices = solarTime <= 12;
eveningIndices = solarTime > 12;
tiegcmMorning = tiegcmGoce270km(morningIndices);
tiegcmEvening = tiegcmGoce270km(eveningIndices);
goceMorning = goceDens270km(morningIndices);
goceEvening = goceDens270km(eveningIndices);
timeMorning = (goceDatenums(morningIndices) - t1) * 1440;
timeEvening = (goceDatenums(eveningIndices) - t1) * 1440;
latMorning = magLat(morningIndices);
latEvening = magLat(eveningIndices);

[timeGrid, latGrid] = meshgrid(1440:5:max(timeMorning)-3*1440, -70:1:70);

goceMorningSurf = griddata(timeMorning, latMorning, goceMorning, timeGrid, latGrid);
goceEveningSurf = griddata(timeEvening, latEvening, goceEvening, timeGrid, latGrid);
tiegcmMorningSurf = griddata(timeMorning, latMorning, tiegcmMorning, timeGrid, latGrid);
tiegcmEveningSurf = griddata(timeEvening, latEvening, tiegcmEvening, timeGrid, latGrid);

% timeGrid = timeGrid / 1440;
% figure('Color', 'w');
% subplot(2,2,1)
% surf(timeGrid, latGrid, tiegcmMorningSurf, 'linestyle', 'none')
% view(2)
% axis('tight')
% colorbar
% xlabel(['Days from the beginning of ', datestr(t1, 'mmmm dd yyyy')]);
% ylabel('Geomagnetic latitude')
% title('Goce Morning Density')

end

function [valAtSatLoc] = interpSatellite(lon, lat, altGrid, modelTime, modelField, satLon, satLat, satAlt, satTime)

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
    tInd = ones(size(lonInd)) * t;
    
    if satAlt < altGrid(i,j,1,t) || satAlt > altGrid(i,j,end,t)
        valAtSatLoc = nan(1);
        return
    end
    
    k = find(altGrid(i,j,:,t) <= satAlt, 1, 'last'); % Toimiiko
    
    nlev = size(altGrid, 3);
    interpAlt = [];
    interpField = [];
    for m = -1:2
        if k + m < 1 || k + m > nlev
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
    
    vals(n) = griddatan(tiegcmCoord, interpField, satCoord, 'linear');
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