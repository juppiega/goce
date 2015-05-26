function readFiles()
%

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool();
end

tic;
[ae, timestampsAeDatenum] = readAeFiles();

[F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, datenumToJulian] = readSolarIndexFiles();

[absB, timestampsAbsBDatenum, akasofuEpsilon, epsilonQualityFlag, timestampsEpsilonDatenum, vBz] ...
    = readSolarWindFiles(timestampsAeDatenum);

[density, longitude, latitude, altitude, solarTime, magneticLatitude, magneticLongitude, magneticLocalTime, crwindEast, crwindNorth, crwindUp, ...
 timestamps10sFixed, timestamps1min, timestamps1minFixed, timestampsDensityDatenum, doy, timestampsAbsB, timestampsEpsilon, firstDatenum, absDensityError, ...
 absWindError, noiseAffected, eclipseAffected, isMorningPass, ionThrusterActive] = readDensityFile(timestampsAeDatenum, timestampsAbsBDatenum, timestampsEpsilonDatenum);

[aeIntegrals] = computeAeIntegrals(ae, timestamps1min, timestamps10sFixed, absB, timestampsAbsB);

[apAll, amAll, dtc] = readApAndDtcFiles(timestampsDensityDatenum);

[F10, F81A, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A] = giveSolarInputForModels(F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, timestampsDensityDatenum);

[apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, ApDaily, am3h, amAver24h, ap, timestamps3h, timestamps3hFixed] =...
    giveApValsForMSIS(apAll, amAll, timestamps10sFixed, timestamps1minFixed, timestamps1min, timestampsDensityDatenum);

[densityNoBg, msisDensityVariableAlt, msisDensity270km, msisDensity270kmNoAp, jb2008DensityVariableAlt, jb2008Density270km, jb2008Density270kmNoDtc, dtm2013Density270km, dtm2013DensityVariableAlt, dtm2013Density270kmNoAm, densityIndex, densityIndex1min, densityIndexNoBg, averagedDensity, averagedDensityNoBg, density3h]  = relateMsisToDensity(density, altitude, datenumToJulian, timestampsDensityDatenum, doy,...
    timestampsAeDatenum, timestamps1minFixed, timestamps3hFixed, solarTime, latitude, longitude, F10, F81A, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, am3h, amAver24h, dtc, ApDaily, apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h);
measuredDensity = density;

[morningTimestamps10s, eveningTimestamps10s, morningLatitude, eveningLatitude, morningDoy, eveningDoy] = ...
    splitBySolarTime(timestamps10sFixed, latitude, doy, solarTime);

[morningGrid, morningBins] = computeTimeCells(morningTimestamps10s, firstDatenum, morningLatitude);
[eveningGrid, eveningBins] = computeTimeCells(eveningTimestamps10s, firstDatenum, eveningLatitude);

%msisSimulatedResidue = msisDensity270km - msisDensity270kmNoAp;

[morningFourierGrid, eveningFourierGrid, morningLatDoy, eveningLatDoy] = computeFourierGrids(morningGrid, morningBins, eveningGrid, ...
    eveningBins, doy, morningLatitude, eveningLatitude, morningDoy, eveningDoy, timestamps10sFixed, densityIndex, aeIntegrals(:,9), median(F10));

if fclose('all') ~= 0
    display('File close unsuccesful - check if some of the files are reserved by another editor.')
end

fprintf('%s\n', 'Saving results to "goceVariables.mat" file')

save('goceVariables.mat', 'ae', '-v7.3')
save('goceVariables.mat', 'aeIntegrals', '-append')
save('goceVariables.mat', 'ap', '-append')
save('goceVariables.mat', 'absB', '-append')
save('goceVariables.mat', 'vBz', '-append')
save('goceVariables.mat', 'densityIndex', '-append')
save('goceVariables.mat', 'densityIndex1min', '-append')
save('goceVariables.mat', 'densityIndexNoBg', '-append')
save('goceVariables.mat', 'averagedDensity', '-append')
save('goceVariables.mat', 'averagedDensityNoBg', '-append')
save('goceVariables.mat', 'density3h', '-append')
save('goceVariables.mat', 'akasofuEpsilon', '-append')
save('goceVariables.mat', 'epsilonQualityFlag', '-append')
save('goceVariables.mat', 'longitude', '-append')
save('goceVariables.mat', 'latitude', '-append')
save('goceVariables.mat', 'altitude', '-append')
save('goceVariables.mat', 'solarTime', '-append')
save('goceVariables.mat', 'crwindEast', '-append')
save('goceVariables.mat', 'crwindNorth', '-append')
save('goceVariables.mat', 'crwindUp', '-append')
save('goceVariables.mat', 'magneticLatitude', '-append')
save('goceVariables.mat', 'magneticLongitude', '-append')
save('goceVariables.mat', 'magneticLocalTime', '-append')
save('goceVariables.mat', 'densityNoBg', '-append')
save('goceVariables.mat', 'msisDensityVariableAlt', '-append')
save('goceVariables.mat', 'msisDensity270km', '-append')
save('goceVariables.mat', 'msisDensity270kmNoAp', '-append')
save('goceVariables.mat', 'dtm2013Density270km', '-append')
save('goceVariables.mat', 'dtm2013Density270kmNoAm', '-append')
save('goceVariables.mat', 'dtm2013DensityVariableAlt', '-append')
save('goceVariables.mat', 'jb2008Density270km', '-append')
save('goceVariables.mat', 'jb2008Density270kmNoDtc', '-append')
save('goceVariables.mat', 'jb2008DensityVariableAlt', '-append')
save('goceVariables.mat', 'morningFourierGrid', '-append')
save('goceVariables.mat', 'eveningFourierGrid', '-append')
save('goceVariables.mat', 'morningLatDoy', '-append')
save('goceVariables.mat', 'eveningLatDoy', '-append')
save('goceVariables.mat', 'morningBins', '-append')
save('goceVariables.mat', 'eveningBins', '-append')
save('goceVariables.mat', 'measuredDensity', '-append')
save('goceVariables.mat', 'absDensityError', '-append')
save('goceVariables.mat', 'absWindError', '-append')
save('goceVariables.mat', 'noiseAffected', '-append')
save('goceVariables.mat', 'eclipseAffected', '-append')
save('goceVariables.mat', 'isMorningPass', '-append')
save('goceVariables.mat', 'ionThrusterActive', '-append')
save('goceVariables.mat', 'timestamps10sFixed', '-append')
save('goceVariables.mat', 'timestamps1min', '-append')
save('goceVariables.mat', 'timestamps1minFixed', '-append')
save('goceVariables.mat', 'timestamps3h', '-append')
save('goceVariables.mat', 'timestamps3hFixed', '-append')
save('goceVariables.mat', 'timestampsAbsB', '-append')
save('goceVariables.mat', 'timestampsEpsilon', '-append')
save('goceVariables.mat', 'firstDatenum', '-append')
save('goceVariables.mat', 'timestampsDensityDatenum', '-append')
save('goceVariables.mat', 'doy', '-append')
save('goceVariables.mat', 'timestampsAeDatenum', '-append')
save('goceVariables.mat', 'timestampsEpsilonDatenum', '-append')

toc;
end

function [ae, timestampsAeDatenum] = readAeFiles()
%
fprintf('%s\n', 'Began reading AE files')

ae = [];
timestampsAeDatenum = [];

aeFiles = dir('ae*');
parfor i = 1:length(aeFiles)
    aeFile = fopen(aeFiles(i).name);
    if aeFile == -1
        error('aeFile open unsuccesful')
    end

    aeData = textscan(aeFile, '%s %s %f %f %f %f %f', 'MultipleDelimsAsOne',1, 'HeaderLines',15);

    ae = [ae; aeData{4}];
    timestampsAeDatenum = [timestampsAeDatenum; datenum(strcat(aeData{1}, aeData{2}), 'yyyy-mm-ddHH:MM:SS.FFF')];
end
[timestampsAeDatenum, indicesToConserve, ~] = unique(timestampsAeDatenum);
ae = ae(indicesToConserve);

interpIndices = ae > 50000;
aeInterp = ae(~interpIndices);
tInterp = timestampsAeDatenum(~interpIndices);
ae = interp1(tInterp, aeInterp, timestampsAeDatenum, 'linear', 'extrap');

end

function [ap, am, dtc] = readApAndDtcFiles(timestampsDensityDatenum)
%
fprintf('%s\n', 'Began reading ap file')

apFile = fopen('apData');
if apFile == -1
    error('ap File open unsuccesful! Check that you have a file "apData" in your WORKING DIRECTORY.')
end

apData = textscan(apFile, '%s %s %s %s %s %s %s %s %s %s %f %f %f %f %f %f %f %f %f', ...
    'delimiter', {' ', '+', '-'}, 'MultipleDelimsAsOne',1, 'HeaderLines',1);
apValues = horzcat(apData{11:18});

[row, ~] = find(isnan(apValues) | apValues >= 1000);
row = unique(row);
for i = 1:length(row)
    frewind(apFile);
    nanLine = textscan(apFile, '%s', 1, 'delimiter', '\n', 'HeaderLines', row(i));
    nanLine = nanLine{1}(1);
    correctValue = nan(52, 1);
    for k = 31:3:52
     correctValue(k) = str2double(horzcat(nanLine{1}(k-2), nanLine{1}(k-1), nanLine{1}(k)));
    end
    correctValue = correctValue(~isnan(correctValue));
    apValues(row(i), :) = correctValue;
end

timestampsAp = datenum(apData{1}, 'yyyymmdd');
apRows = find(timestampsAp > timestampsDensityDatenum(1) - 4 & timestampsAp <= timestampsDensityDatenum(end));
ap = reshape(apValues(apRows, :)', numel(apValues(apRows, :)), 1);

firstHour = round(str2double(datestr(timestampsDensityDatenum(1), 'HH')));
ap(1:firstHour/3) = [];


amFile = fopen('am_file_spider.dat');
if amFile == -1
    error('am File open unsuccesful! Check that you have a file "m_file_spider.data" in your WORKING DIRECTORY.')
end

amData = textscan(amFile, '%f %f %f %f %f %f %f %f %f %f %f %f %f', 'MultipleDelimsAsOne',1);
amValues = horzcat(amData{end-8:end});
[row, ~] = find(isnan(amValues));
amValues(row,2:9) = amValues(row,1:8);
amValues = amValues(:,2:9);

year = mod(amData{2}, 10000);
month = amData{1};
day = floor(amData{2}/10000);
timestampsAm = datenum(year, month, day);
amValues(amValues < 0) = 0.0;
amRows = find(timestampsAm > timestampsDensityDatenum(1) - 4 & timestampsAm <= timestampsDensityDatenum(end));
am = reshape(amValues(amRows, :)', numel(amValues(amRows, :)), 1);
am(1:firstHour/3) = [];


dtcFile = fopen('DTCFILE.TXT');
if dtcFile == -1
    error('Dtc file open unsuccesful! Check that you have a file "DTCFILE.TXT" in your WORKING DIRECTORY.')
end

dtcData = textscan(dtcFile, '%s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d', 'MultipleDelimsAsOne',1, 'EmptyValue', 9999);
dtcData = horzcat(dtcData{2:end});
row = find(dtcData(:,end) > 999);
for i = 1:length(row)
    frewind(dtcFile);
    nanLine = textscan(dtcFile, '%s', 1, 'delimiter', '\n', 'HeaderLines', row(i) - 1);
    nanLine = nanLine{1}(1);
    yy = str2double(nanLine{1}(5:6));
    ddd = str2double(nanLine{1}(7:9));
    dtcData(row(i), :) = [yy ddd dtcData(row(i), 2:end-1)];
end

year = dtcData(:,1);
lastCentury = year >= 97;
year(lastCentury) = year(lastCentury) + 1900;
year(~lastCentury) = year(~lastCentury) + 2000;
dayOfYear = dtcData(:,2);
dtcDatenum = datenum(double(year), 1, 1) + double(dayOfYear) - 1;
dtcRows = find(dtcDatenum > timestampsDensityDatenum(1) - 1 & dtcDatenum <= timestampsDensityDatenum(end));
dtcValues = dtcData(:,3:end);
dtc = reshape(dtcValues(dtcRows, :)', numel(dtcValues(dtcRows, :)), 1);

numOfDays = length(dtcRows);
firstDtcDay = dtcDatenum(dtcRows(1));
timestampsDtc = (0 : 3600 : numOfDays * 86400 - 1)';
timestampsFixed = round((timestampsDensityDatenum - firstDtcDay) * 86400);
dtc = interp1(timestampsDtc, double(dtc), timestampsFixed, 'linear', 'extrap');

end

function [F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, datenumToJulian] = readSolarIndexFiles
%
fprintf('%s\n', 'Began reading solar index files')

solarFile = fopen('SOLFSMY.TXT');
if solarFile == -1
    error('Solar index file open unsuccesful! Check that you have a file "SOLFSMY.TXT" in your WORKING DIRECTORY.')
end

solarData = textscan(solarFile, '%d %d %f %f %f %f %f %f %f %f %f %s', 'MultipleDelimsAsOne',1, 'CommentStyle','#');
F10 = solarData{4};
F81A = solarData{5};
S10 = solarData{6};
S81A = solarData{7};
M10 = solarData{8};
M81A = solarData{9};
Y10 = solarData{10};
Y81A = solarData{11};

SnoNans = find(S10 > 1);
S81noNans = find(S81A > 1);
S10 = interp1(SnoNans, S10(SnoNans), 1:length(S10), 'linear', 'extrap')';
S81A = interp1(S81noNans, S81A(S81noNans), 1:length(S81A), 'linear', 'extrap')';

year = solarData{1};
dayOfYear = solarData{2};
indexDatenums = datenum(double(year), 1, 1) + double(dayOfYear) - 1;
julianDate = solarData{3};
datenumToJulian = julianDate(1) - indexDatenums(1) - 0.5;


solarFile = fopen('proxies_unadjusted.dat');
if solarFile == -1
    error('Solar index file open unsuccesful! Check that you have a file "proxies_unadjusted.dat" in your WORKING DIRECTORY.')
end

solarData = textscan(solarFile, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'MultipleDelimsAsOne',1, 'CommentStyle','//');
F30 = solarData{14};
F30A = solarData{15};

year = solarData{1};
month = solarData{2};
day = solarData{3};
F30datenum = datenum(year, month, day);
F30 = interp1(F30datenum, F30, indexDatenums, 'nearest', 'extrap');
F30A = interp1(F30datenum, F30A, indexDatenums, 'nearest', 'extrap');

end

function [absB, timestampsAbsBDatenum, akasofuEpsilon, epsilonQualityFlag, timestampsEpsilonDatenum, vBz] = ...
    readSolarWindFiles(timestampsAeDatenum)
%
fprintf('%s\n', 'Began reading IMF and SW files')

imfFiles = dir('ac_h0s_mfi*');

if isempty(imfFiles)
    error('No IMF B files found! Make sure there are files beginning "ac_h0s_mfi_*" in the WORKING DIRECTORY')
end

timestampsBDatenum = [];
Bx = [];
By = [];
Bz = [];
for i = 1:length(imfFiles)
     data = cdfread(imfFiles(i).name, 'CombineRecords',true, 'ConvertEpochToDatenum',true);
     timestampsBDatenum = [timestampsBDatenum; data{1}];
     Bx = [Bx; data{2}(1:3:end)'];
     By = [By; data{2}(2:3:end)'];
     Bz = [Bz; data{2}(3:3:end)'];
end
[timestampsBDatenum, indicesToConserve, ~] = unique(timestampsBDatenum);
Bx = Bx(indicesToConserve);
By = By(indicesToConserve);
Bz = Bz(indicesToConserve);
B2 = Bx.^2 + By.^2 + Bz.^2;
theta = atan2(By, Bz);

indicesToConserve = find(Bx > -1e30 & By > -1e30 & Bz > -1e30);

tBNoNans = timestampsBDatenum(indicesToConserve);
thetaNoNans = theta(indicesToConserve);
B2NoNans = B2(indicesToConserve);
BzNoNans = Bz(indicesToConserve);

tBToInterpolate = timestampsAeDatenum(timestampsAeDatenum>=min(timestampsBDatenum) & timestampsAeDatenum<=max(timestampsBDatenum));
tThetaToInterpolate = tBToInterpolate;
B2Interpolated = interp1(tBNoNans, B2NoNans, tBToInterpolate, 'linear', 'extrap');
BzInterpolated = interp1(tBNoNans, BzNoNans, tBToInterpolate, 'linear', 'extrap');
thetaInterpolated = interp1(tBNoNans, thetaNoNans, tThetaToInterpolate, 'linear', 'extrap');

secondsInDay = 24 * 60 * 60;
absB = sqrt(B2Interpolated);
timestampsAbsBDatenum = tBToInterpolate;

solarWindFiles = dir('ac_h0s_swe*');

if isempty(solarWindFiles)
    error('No Solar Wind files found! Make sure there are files beginning "ac_h0s_swe_*" in the WORKING DIRECTORY')
end

timestampsVDatenum = [];
Vx = [];
Vy = [];
Vz = [];
for i = 1:length(solarWindFiles)
     data = cdfread(solarWindFiles(i).name, 'CombineRecords',true, 'ConvertEpochToDatenum',true);
     timestampsVDatenum = [timestampsVDatenum; data{1}];
     Vx = [Vx; data{2}(1:3:end)'];
     Vy = [Vy; data{2}(2:3:end)'];
     Vz = [Vz; data{2}(3:3:end)'];
end
[timestampsVDatenum, indicesToConserve, ~] = unique(timestampsVDatenum);
Vx = Vx(indicesToConserve);
Vy = Vy(indicesToConserve);
Vz = Vz(indicesToConserve);
V = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
indicesToConserve = find(Vx > -1e30 & Vy > -1e30 & Vz > -1e30);
tVNoNans = timestampsVDatenum(indicesToConserve);
VNoNans = V(indicesToConserve);

tVToInterpolate = timestampsAeDatenum(timestampsAeDatenum>=min(timestampsVDatenum) & timestampsAeDatenum<=max(timestampsVDatenum));
VInterpolated = interp1(tVNoNans, VNoNans, tVToInterpolate, 'linear', 'extrap');
% timesToNeareastNonNan = arrayfun(@(x) min(abs(tVNoNans - x)), tVToInterpolate);
% epsilonQualityFlag = ones(length(tVToInterpolate), 1);
% epsilonQualityFlag(timesToNeareastNonNan > oneMinuteInDays) = 0;  
commonIndices = timestampsVDatenum >= min(tVToInterpolate) & timestampsVDatenum <= max(tVToInterpolate);
badDataIndices = ~ismember(timestampsVDatenum(commonIndices), tVNoNans);
epsilonQualityFlag = ones(length(tVToInterpolate), 1);
epsilonQualityFlag(badDataIndices) = 0; 

VInterpolated = VInterpolated(ismember(tVToInterpolate, tBToInterpolate)) * 1000;
epsilonQualityFlag = epsilonQualityFlag(ismember(tVToInterpolate, tBToInterpolate));
B2Interpolated = B2Interpolated(ismember(tBToInterpolate, tVToInterpolate)) * 1e-9;
BzInterpolated = BzInterpolated(ismember(tBToInterpolate, tVToInterpolate)) * 1e-9;
thetaInterpolated = thetaInterpolated(ismember(tThetaToInterpolate, tVToInterpolate));

akasofuEpsilon = VInterpolated .* B2Interpolated .* sin(thetaInterpolated * 0.5).^4;
vBz = VInterpolated .* BzInterpolated;

tVToInterpolate = tVToInterpolate(ismember(tVToInterpolate, tBToInterpolate));
timestampsEpsilonDatenum = tVToInterpolate;

end

function [density, longitude, latitude, altitude, solarTime, magneticLatitude, magneticLongitude, magneticLocalTime, crwindEast, crwindNorth, crwindUp,...
    timestamps10sFixed, timestamps1min, timestamps1minFixed, timestampsDensityDatenum, doy, timestampsAbsB, timestampsEpsilon, firstDatenum,...
    absDensityError, absWindError, noiseAffected, eclipseAffected, isMorningPass, ionThrusterActive] = ...
    readDensityFile(timestampsAeDatenum, timestampsAbsBDatenum, timestampsEpsilonDatenum)
%
fprintf('%s\n', 'Began reading density files')

densityFiles = dir('goce_dens*');

density = [];
longitude = [];
latitude = [];
altitude = [];
solarTime = [];
crwindEast = [];
crwindNorth = [];
crwindUp = [];
absDensityError = [];
absWindError = [];
noiseAffected = [];
eclipseAffected = [];
isMorningPass = [];
ionThrusterActive = [];

timestampsDensityDatenum = [];
doy = [];

secondsInDay = 60 * 60 * 24;
parfor i = 1:length(densityFiles)
    densityFile = fopen(densityFiles(i).name);

    if densityFile == -1
        error('densityFile open unsuccesful')
    end

    densityData = textscan(densityFile, '%s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'commentStyle','#');

    density = [density; densityData{9} * power(10, 11)];
    longitude = [longitude; densityData{5}];
    latitude = [latitude; densityData{6}];
    altitude = [altitude; densityData{4}];
    solarTime = [solarTime; densityData{7}];
    crwindEast = [crwindEast; densityData{10}];
    crwindNorth = [crwindNorth; densityData{11}];
    crwindUp = [crwindUp; densityData{12}];
    absDensityError = [absDensityError; densityData{13} * power(10, 11)];
    absWindError = [absWindError; densityData{14}];
    noiseAffected = [noiseAffected; densityData{15}];
    eclipseAffected = [eclipseAffected; densityData{16}];
    isMorningPass = [isMorningPass; densityData{17}];
    ionThrusterActive = [ionThrusterActive; densityData{18}];

    timestampsThisFile = datenum(strcat(densityData{1}, densityData{2}), 'yyyy-mm-ddHH:MM:SS.FFF');
    timestampsDensityDatenum = [timestampsDensityDatenum; timestampsThisFile];
    
    doyThisFile = timestampsThisFile - datenum(datestr(timestampsThisFile, 'yyyy'), 'yyyy');
    doy = [doy; doyThisFile];
end

indicesToRemove = roundTimestampsToBeginFromNearestThreeHourMark(timestampsDensityDatenum);

[~, uniqueIndices, ~] = unique(timestampsDensityDatenum);
uniqueIndices = ismember(1:length(timestampsDensityDatenum), uniqueIndices)';
indicesToConserve = ~indicesToRemove & uniqueIndices;

timestampsDensityDatenum = timestampsDensityDatenum(indicesToConserve);
doy = doy(indicesToConserve);
density = density(indicesToConserve);
longitude = longitude(indicesToConserve);
latitude = latitude(indicesToConserve);
altitude = altitude(indicesToConserve);
solarTime = solarTime(indicesToConserve);
crwindEast = crwindEast(indicesToConserve);
crwindNorth = crwindNorth(indicesToConserve);
crwindUp = crwindUp(indicesToConserve);
absDensityError = absDensityError(indicesToConserve);
absWindError = absWindError(indicesToConserve);
noiseAffected = logical(round(noiseAffected(indicesToConserve)));
eclipseAffected = logical(round(eclipseAffected(indicesToConserve)));
isMorningPass = logical(round(isMorningPass(indicesToConserve)));
ionThrusterActive = ~logical(round(ionThrusterActive(indicesToConserve)));

[magneticLatitude, magneticLongitude] = convertToMagneticCoordinates(latitude, longitude, altitude);

timestamps10sFixed = timestampsDensityDatenum * secondsInDay;
timestamps1minFixed = timestamps10sFixed(ismember(timestampsDensityDatenum, timestampsAeDatenum));
firstDatenum = datenum(datestr(timestampsDensityDatenum(1), 'yyyy-mm-dd'), 'yyyy-mm-dd');
timestamps10sFixed = timestamps10sFixed - firstDatenum * secondsInDay;
timestamps1minFixed = timestamps1minFixed - firstDatenum * secondsInDay;
timestamps10sFixed = round(timestamps10sFixed);
timestamps1minFixed = round(timestamps1minFixed);
timestamps1min = round((timestampsAeDatenum - firstDatenum) * secondsInDay);

timestampsAbsB = round((timestampsAbsBDatenum - firstDatenum) * secondsInDay);
timestampsEpsilon = round((timestampsEpsilonDatenum - firstDatenum) * secondsInDay);

[magneticLocalTime] = computeMagneticTime(magneticLongitude, doy, timestampsDensityDatenum);

end

function indicesToRemove = roundTimestampsToBeginFromNearestThreeHourMark(timestampsDensityDatenum)
%

[y, m, d, hours, minutes, seconds] = datevec(timestampsDensityDatenum(1));
hours = hours + minutes / 60 + seconds / 3600;
hours = 3 * ceil(hours / 3);
nearest3hTime = datenum(y, m, d, hours, 0, 0);

if isempty(find(timestampsDensityDatenum == nearest3hTime, 1))
    while 1
        
        hours = hours + 3;
        nearest3hTime = datenum(y, m, d, hours, 0, 0);
        
        if ~isempty(find(timestampsDensityDatenum == nearest3hTime, 1))
            break;
        end
    end
end

indicesToRemove = timestampsDensityDatenum < nearest3hTime;

end

function [aeIntegrals] = computeAeIntegrals(ae, timestampsAe, timestamps10s, absB, timestampsAbsB)
%
absB = interp1(timestampsAbsB, absB, timestampsAe, 'linear', 0);
aeIntegrals = zeros(length(timestamps10s), 122);
cumulativeAe = cumsum(ae);
cumulativeAbsB = cumsum(absB);

% thirtyMinutes = 30;
% firstAeInt = cumulativeAe(thirtyMinutes + 1 : end) - cumulativeAe(1 : end - thirtyMinutes);
% firstAeTime = timestampsAe(thirtyMinutes + 1 : end);
% aeIntegrals(:,1) = interp1(firstAeTime, firstAeInt, timestampsAe, 'linear', 0);
% 
% firstAbsBint = cumulativeAbsB(thirtyMinutes + 1 : end) - cumulativeAbsB(1 : end - thirtyMinutes);
% aeIntegrals(:,62) = interp1(firstAeTime, firstAbsBint, timestampsAe, 'linear', 0);

oneHour = 60;
for i = [2 4 8 16 21 30 40 50 60]
    lag = i * oneHour;
    aeInt = cumulativeAe(lag + 1 : end) - cumulativeAe(1 : end - lag);
    aeTime = timestampsAe(lag + 1 : end);
    aeIntegrals(:,i) = interp1(aeTime, aeInt, timestamps10s, 'linear', 0);
end

% for i = [2 4 8 16 30 37 50 60]
%     lag = i * oneHour;
%     absBInt = cumulativeAbsB(lag + 1 : end) - cumulativeAbsB(1 : end - lag);
%     absBTime = timestampsAe(lag + 1 : end);
%     aeIntegrals(:,i + 61) = interp1(absBTime, absBInt, timestampsAe, 'linear', 0);
% end

columnsToConserve = sum(aeIntegrals) > 1;
aeIntegrals = aeIntegrals(:,columnsToConserve);

end

function [F10, F81A, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A] = giveSolarInputForModels(F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, timestampsDensityDatenum)
%
fprintf('%s\n', 'Computing solar index values for models')

msisF81A = interp1(indexDatenums, F81A, timestampsDensityDatenum, 'linear', 100);
F30A = interp1(indexDatenums, F30A, timestampsDensityDatenum, 'linear', 100);

F10 = interp1(indexDatenums + 1, F10, timestampsDensityDatenum, 'linear', 100);
F81A = interp1(indexDatenums + 1, F81A, timestampsDensityDatenum, 'linear', 100);
S10 = interp1(indexDatenums + 1, S10, timestampsDensityDatenum, 'linear', 100);
S81A = interp1(indexDatenums + 1, S81A, timestampsDensityDatenum, 'linear', 100);
F30 = interp1(indexDatenums + 1, F30, timestampsDensityDatenum, 'linear', 100);

M10 = interp1(indexDatenums + 2, M10, timestampsDensityDatenum, 'linear', 100);
M81A = interp1(indexDatenums + 2, M81A, timestampsDensityDatenum, 'linear', 100);

Y10 = interp1(indexDatenums + 5, Y10, timestampsDensityDatenum, 'linear', 100);
Y81A = interp1(indexDatenums + 5, Y81A, timestampsDensityDatenum, 'linear', 100);

% indexBeginDay = floor(timestampsDensityDatenum(1)) - 6;
% indexEndDay = floor(timestampsDensityDatenum(end));
% indicesToConserve = indexDatenums >= indexBeginDay & indexDatenums <= indexEndDay;
% F10 = F10(indicesToConserve);
% F81A = F81A(indicesToConserve);
% S10 = S10(indicesToConserve);
% S81A = S81A(indicesToConserve);
% M10 = M10(indicesToConserve);
% M81A = M81A(indicesToConserve);
% Y10 = Y10(indicesToConserve);
% Y81A = Y81A(indicesToConserve);
% F30 = F30(indicesToConserve);
% F30A = F30A(indicesToConserve);

% secondsInDay = 86400;
% timestampsFixed = round((timestampsDensityDatenum - indexBeginDay(1)) * secondsInDay);
% numOfDays = indexEndDay - indexBeginDay + 1;
% timestampsMsisF81 = (0 : 10 : numOfDays * secondsInDay - 1)';
% timestampsFS = (0 : 10 : numOfDays * secondsInDay - 1)' + secondsInDay;
% timestampsM = (0 : 10 : numOfDays * secondsInDay - 1)' + 2 * secondsInDay;
% timestampsY = (0 : 10 : numOfDays * secondsInDay - 1)' + 5 * secondsInDay;

% pointsInDay = 8640;
% F10 = reshape(repmat((F10)', pointsInDay, 1), [], 1);
% msisF81A = reshape(repmat((F81A)', pointsInDay, 1), [], 1);
% F81A = reshape(repmat((F81A)', pointsInDay, 1), [], 1);
% S10 = reshape(repmat((S10)', pointsInDay, 1), [], 1);
% S81A = reshape(repmat((S81A)', pointsInDay, 1), [], 1);
% M10 = reshape(repmat((M10)', pointsInDay, 1), [], 1);
% M81A = reshape(repmat((M81A)', pointsInDay, 1), [], 1);
% Y10 = reshape(repmat((Y10)', pointsInDay, 1), [], 1);
% Y81A = reshape(repmat((Y81A)', pointsInDay, 1), [], 1);
% F30 = reshape(repmat((F30)', pointsInDay, 1), [], 1);
% F30A = reshape(repmat((F30A)', pointsInDay, 1), [], 1);
% 
% F81Indices = ismember(timestampsMsisF81, timestampsFixed);
% msisF81A = msisF81A(F81Indices);
% F30A = F30A(F81Indices);
% 
% FSindices = ismember(timestampsFS, timestampsFixed);
% F10 = F10(FSindices);
% F81A = F81A(FSindices);
% S10 = S10(FSindices);
% S81A = S81A(FSindices);
% F30 = F30(FSindices);
% 
% Mindices = ismember(timestampsM, timestampsFixed);
% M10 = M10(Mindices);
% M81A = M81A(Mindices);
% 
% Yindices = ismember(timestampsY, timestampsFixed);
% Y10 = Y10(Yindices);
% Y81A = Y81A(Yindices);

end

function [apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, ApDaily, am3h, amAver24h, apAllFixed, timestamps3h, timestamps3hFixed] = ...
    giveApValsForMSIS(apAll, amAll, timestamps10sFixed, timestamps1minFixed, timestamps1min, timestampsDensityDatenum)
% 
fprintf('%s\n', 'Computing ap values for msis')

threeHinSec = 3 * 60 * 60;
secondsInDay = 24 * 60 * 60;
firstDatenum = datenum(datestr(timestampsDensityDatenum(1), 'yyyy-mm-dd'));
lastDatenum = datenum(datestr(timestampsDensityDatenum(end), 'yyyy-mm-dd'));
timestamps3h = -3 * 8 * threeHinSec + threeHinSec : threeHinSec : ceil(lastDatenum - firstDatenum + 1) * secondsInDay;
firstHour = round(str2double(datestr(timestampsDensityDatenum(1), 'HH')));
timestamps3h(1: firstHour / 3) = [];
timestamps10s = min(timestamps10sFixed) : 10 : max(timestamps10sFixed);
nextThreeHStamp = threeHinSec * ceil(timestamps10s / threeHinSec);

ap10s = reshape(repmat(apAll', 180 * 6, 1),[],1);
apTime10s = reshape(repmat(timestamps3h, 180 * 6, 1),[],1);
ap10s(1) = [];
apTime10s(1) = [];
apNow = ap10s(ismember(apTime10s, nextThreeHStamp));
ap3h = ap10s(ismember(apTime10s, nextThreeHStamp - threeHinSec));
ap6h = ap10s(ismember(apTime10s, nextThreeHStamp - 2 * threeHinSec));
ap9h = ap10s(ismember(apTime10s, nextThreeHStamp - 3 * threeHinSec));
ap10sSmoothed = reshape(repmat(smooth(apAll, 9)', 180 * 6, 1),[],1);
ap10sSmoothed(1) = [];
apAver12To33h = ap10sSmoothed(ismember(apTime10s, nextThreeHStamp - 8 * threeHinSec));
apAver36To57h = ap10sSmoothed(ismember(apTime10s, nextThreeHStamp - 16 * threeHinSec));
ApDaily = ap10sSmoothed(ismember(apTime10s, nextThreeHStamp - 4 * threeHinSec));

finalIndices = ismember(timestamps10s, timestamps10sFixed);
apNow = apNow(finalIndices);
ap3h = ap3h(finalIndices);
ap6h = ap6h(finalIndices);
ap9h = ap9h(finalIndices);
apAver12To33h = apAver12To33h(finalIndices);
apAver36To57h = apAver36To57h(finalIndices);
ApDaily = ApDaily(finalIndices);

am3h = interp1(timestamps3h + threeHinSec, amAll, timestamps10sFixed, 'linear', 'extrap');

averagingMatrix = nan(8, length(amAll) - 7);

for i = 1:8
    lastIndex = length(amAll) - 8 + i;
    averagingMatrix(i,:) = amAll(i:lastIndex)';
end
averagingMatrix = [repmat(averagingMatrix(:,1),1,7) averagingMatrix];
amAver24h = trapz(averagingMatrix)' / 8;

amAver24h = interp1(timestamps3h, amAver24h, timestamps10sFixed, 'linear', 'extrap');

% apNow = interp1(timestamps3h, apAll, timestamps10sFixed, 'linear', 'extrap');
% ap3h = interp1(timestamps3h + threeHinSec, apAll, timestamps10sFixed, 'linear', 'extrap');
% ap6h = interp1(timestamps3h + 2 * threeHinSec, apAll, timestamps10sFixed, 'linear', 'extrap');
% ap9h = interp1(timestamps3h + 3 * threeHinSec, apAll, timestamps10sFixed, 'linear', 'extrap');
% 
% ap24hSmooth = smooth(apAll,8);
% ApDaily = interp1(timestamps3h, ap24hSmooth, timestamps10sFixed, 'linear', 'extrap');
% apAver12To33h = interp1(timestamps3h + 7.5 * threeHinSec, ap24hSmooth, timestamps10sFixed, 'linear', 'extrap');
% apAver36To57h = interp1(timestamps3h + 15.5 * threeHinSec, ap24hSmooth, timestamps10sFixed, 'linear', 'extrap');

timestamps3h(timestamps3h <= firstHour * threeHinSec) = nan(1);
apAllFixed = apAll(ismember(timestamps3h, timestamps1min));
timestamps3hFixed = timestamps3h(ismember(timestamps3h, timestamps1minFixed))';
timestamps3h = timestamps3h(ismember(timestamps3h, timestamps1min))';

end

function [correctedDensity, msisDensityVariableAlt, msisDensity270km, msisDensity270kmNoAp, jb2008DensityVariableAlt, jb2008Density270km, jb2008Density270kmNoDtc, dtm2013Density270km, dtm2013DensityVariableAlt, dtm2013Density270kmNoAm, densityIndex, densityIndex1min, densityIndexNoBg, averagedDensity, averagedDensityNoBg, density3h] = relateMsisToDensity(density, altitude, datenumToJulian, timestampsDatenum, doy,...
          timestampsAeDatenum, timestamps1minFixed, timestamps3hFixed, solarTime, latitude, longtitude, F10, F81A, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, am3h, amAver24h, dtc, ApDaily, apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h)
      
fprintf('%s\n', 'Computing normalized densities with models. This may take even half an hour.')

msisDensity270km = nan(size(density));
msisDensityVariableAlt = nan(size(density));
msisDensity270kmNoAp = nan(size(density));

modelingIndices = 1:length(density);

altitudeInKm = altitude / 1000;

doyDecimal = doy;
doy = ceil(doy);
doy(doy == 0) = 1;

% second of the day
secondsInDay = 24 * 60 * 60;
seconds = (timestampsDatenum - floor(timestampsDatenum)) * secondsInDay;

julianDay = timestampsDatenum + datenumToJulian;

targetCount = round(length(modelingIndices) / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running MSIS, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

parfor i = modelingIndices
    [~,~,~,~,~,msisDensityVariableAlt(i),~,~,~,~,~]...
      =nrlmsise_mex(doy(i),seconds(i),altitudeInKm(i),latitude(i),longtitude(i),solarTime(i),msisF81A(i),F10(i),...
      ApDaily(i),apNow(i),ap3h(i),ap6h(i),ap9h(i),apAver12To33h(i),apAver36To57h(i));

    [~,~,~,~,~,msisDensity270km(i),~,~,~,~,~]...
      =nrlmsise_mex(doy(i),seconds(i),270,latitude(i),longtitude(i),solarTime(i),msisF81A(i),F10(i),...
      ApDaily(i),apNow(i),ap3h(i),ap6h(i),ap9h(i),apAver12To33h(i),apAver36To57h(i));

    [~,~,~,~,~,msisDensity270kmNoAp(i),~,~,~,~,~]...
      =nrlmsise_mex(doy(i),seconds(i),270,latitude(i),longtitude(i),solarTime(i),msisF81A(i),F10(i), 3);

    if mod(i, 10000) == 0
     p.progress;
    end
end
p.stop;

jb2008Density270km = nan(size(density));
jb2008DensityVariableAlt = nan(size(density));
jb2008Density270kmNoDtc = nan(size(density));

p = TimedProgressBar( targetCount, barWidth, ...
                    'Running JB2008, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
                
parfor i = modelingIndices
    [~,~,jb2008DensityVariableAlt(i)] = jb2008_mex(julianDay(i), altitudeInKm(i), latitude(i), longtitude(i), F10(i), F81A(i), S10(i),...
        S81A(i), M10(i), M81A(i), Y10(i), Y81A(i), dtc(i));
    
    [~,~,jb2008Density270km(i)] = jb2008_mex(julianDay(i), 270, latitude(i), longtitude(i), F10(i), F81A(i), S10(i),...
        S81A(i), M10(i), M81A(i), Y10(i), Y81A(i), dtc(i));
    
    [~,~,jb2008Density270kmNoDtc(i)] = jb2008_mex(julianDay(i), 270, latitude(i), longtitude(i), F10(i), F81A(i), S10(i),...
        S81A(i), M10(i), M81A(i), Y10(i), Y81A(i), 0);
    
    if mod(i, 10000) == 0
     p.progress;
    end
end
p.stop;

dtm2013Density270km = nan(size(density));
dtm2013DensityVariableAlt = nan(size(density));
dtm2013Density270kmNoAm = nan(size(density));

p = TimedProgressBar( targetCount, barWidth, ...
                    'Running DTM-2013, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
                
parfor i = modelingIndices
    [~, ~, dtm2013DensityVariableAlt(i),~,~,~] = dtm2013_mex(doyDecimal(i), altitudeInKm(i), latitude(i), longtitude(i), ...
        solarTime(i), F30(i), F30A(i), am3h(i), amAver24h(i));
    
    [~, ~, dtm2013Density270km(i),~,~,~] = dtm2013_mex(doyDecimal(i), 270, latitude(i), longtitude(i), ...
        solarTime(i), F30(i), F30A(i), am3h(i), amAver24h(i));
    
    [~, ~, dtm2013Density270kmNoAm(i),~,~,~] = dtm2013_mex(doyDecimal(i), 270, latitude(i), longtitude(i), ...
        solarTime(i), F30(i), F30A(i), 0.0, 0.0);
    
    if mod(i, 10000) == 0
     p.progress;
    end
end
p.stop;

dtm2013DensityVariableAlt = dtm2013DensityVariableAlt * power(10, 14);
dtm2013Density270km = dtm2013Density270km * power(10, 14);
dtm2013Density270kmNoAm = dtm2013Density270kmNoAm * power(10, 14);

jb2008DensityVariableAlt = jb2008DensityVariableAlt * power(10, 11);
jb2008Density270km = jb2008Density270km * power(10, 11);
jb2008Density270kmNoDtc = jb2008Density270kmNoDtc * power(10, 11);

correctedDensity = density .* msisDensity270km ./ msisDensityVariableAlt;

msisDensityVariableAlt = msisDensityVariableAlt * power(10, 14);
msisDensity270km = msisDensity270km * power(10, 14);
msisDensity270kmNoAp = msisDensity270kmNoAp * power(10, 14);

correctionFactor = mean(msisDensityVariableAlt ./ density);
densityIndex = (correctionFactor .* correctedDensity) - msisDensity270kmNoAp;
densityIndex1min = smooth(densityIndex, 7);
densityIndex1min = densityIndex1min(ismember(timestampsDatenum, timestampsAeDatenum));
densityIndexNoBg = removePeriodicBackground(densityIndex1min, 30, 1, 0);
densityIndexNoBg = normalize(densityIndexNoBg, densityIndex1min);

averagedDensity = smooth(correctedDensity, 7);
averagedDensity = averagedDensity(ismember(timestampsDatenum, timestampsAeDatenum));
averagedDensityNoBg = removePeriodicBackground(averagedDensity, 125, 1, 0);
averagedDensityNoBg = normalize(averagedDensityNoBg, averagedDensity);

density3h = smooth(averagedDensityNoBg, 179);
density3h = density3h(find(ismember(timestamps1minFixed, timestamps3hFixed)) - 90);

end

function [morningTimestamps10s, eveningTimestamps10s, morningLatitude, eveningLatitude, morningDoy, eveningDoy] = ...
    splitBySolarTime(timestamps10s, latitude, doy, solarTime)
%

morningIndices = find(solarTime <= 12);
eveningIndices = find(solarTime > 12);

morningTimestamps10s = timestamps10s(morningIndices);
morningLatitude = latitude(morningIndices);
morningDoy = doy(morningIndices);

eveningTimestamps10s = timestamps10s(eveningIndices);
eveningLatitude = latitude(eveningIndices);
eveningDoy = doy(eveningIndices);

end

function [timeByLatitude, latitudeBins] = computeTimeCells(timestamps10s, firstDatenum, latitude)
%

analysisLimit = datenum('2013-06-16', 'yyyy-mm-dd');

timeInDays = timestamps10s / 86400 + firstDatenum;
indicesToConserve = timeInDays < analysisLimit;
allTimestamps10s = timestamps10s;
timestamps10s = timestamps10s(indicesToConserve);
latitude = latitude(ismember(allTimestamps10s, timestamps10s));

[minLatitude, maxLatitude] = findInterpolationLimits(latitude);


latitudeBins = minLatitude + 1 : 3 : maxLatitude - 1;
if maxLatitude - 1 > max(latitudeBins)
    latitudeBins(end + 1) = max(latitudeBins) + 3;
end

timeByLatitude = cell(length(latitudeBins), 1);
indicesToRemove = false(length(latitudeBins),1);
parfor i = 1:length(latitudeBins)
    indices = (latitudeBins(i) - 1.5 < latitude & latitude <= latitudeBins(i) + 1.5);
    if isempty(find(indices, 1))
        indicesToRemove(i) = 1;
        continue;
    end
    timeByLatitude{i} = timestamps10s(indices);
end

latitudeBins = latitudeBins(~indicesToRemove);
timeByLatitude = timeByLatitude(~indicesToRemove);

end

function [morningFourierGrid, eveningFourierGrid, morningLatDoy, eveningLatDoy] = computeFourierGrids(morningGrid, morningBins, eveningGrid, eveningBins, doy, morningLatitude, eveningLatitude, morningDoy, eveningDoy, timestamps10sFixed, densityIndex, aeIntegral, F107median)
%

fprintf('%s\n', 'Computing fourier fits.')

targetCount = 2*length(morningBins);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running Parametrizations, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
%aeLimit = aeIntegral < quantile(aeIntegral, 0.99);

binWidth = 10;
doy1day = (1:365)';

parfor i = 1:length(morningBins)
    timeThisLat = morningGrid{i};
    thisLatIndices = ismember(timestamps10sFixed, timeThisLat);
    densityThisLat = densityIndex(thisLatIndices);
    
    doyThisLat = ceil(doy(thisLatIndices));
    
    thisLatDens = zeros(length(doy1day),1);
    for k = 1:length(doy1day)
        indices = doyThisLat == doy1day(k);
        thisLatDens(k,:) = mean(densityThisLat(indices));
    end
   
    morningDensities(:,i) = thisLatDens; 
    
    p.progress;
end

parfor i = 1:length(eveningBins)
    timeThisLat = eveningGrid{i};
    thisLatIndices = ismember(timestamps10sFixed, timeThisLat);
    densityThisLat = densityIndex(thisLatIndices);

    doyThisLat = ceil(doy(thisLatIndices));
    
    thisLatDens = zeros(length(doy1day),1);
    for k = 1:length(doy1day)
        indices = doyThisLat == doy1day(k);
        thisLatDens(k,:) = mean(densityThisLat(indices));
    end

    eveningDensities(:,i) = thisLatDens;
    
    p.progress;
end

% latitude = morningBins;
% morningDensities = zeros(length(doy1day), length(latitude));
% eveningDensities = zeros(length(doy1day), length(latitude));
% morningSolarTime = 7.0;
% eveningSolarTime = 19.0;
% seconds = 0 : 3 * 60 * 60 : 24 * 60 * 60 - 1;
% morningMsis270km = zeros(length(seconds), 1);
% morningMsisNoAp = zeros(length(seconds), 1);
% eveningMsis270km = zeros(length(seconds), 1);
% eveningMsisNoAp = zeros(length(seconds), 1);
% 
% F107A = F107median; % = 107
% F107 = F107median;
% ApDaily = 27;
% apNow = 80;
% ap3h = 207;
% ap6h = 48;
% ap9h = 15;
% apAver12To33h = 7;
% apAver36To57h = 7;
% for i = 1:length(latitude)
%     for j = 1:length(doy1day)
%         for k = 1:length(seconds)
%             morningLongitude = 180 * (morningSolarTime - (seconds(k)/3600)) / 12;
%             eveningLongitude = 180 * (eveningSolarTime - (seconds(k)/3600)) / 12;
%             morningLongitude(morningLongitude < -180) = morningLongitude + 360;
%             eveningLongitude(eveningLongitude > 180) = eveningLongitude - 360;
%             
%             [~,~,~,~,~,morningMsis270km(k),~,~,~,~,~]...
%             =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),morningLongitude,morningSolarTime,F107A,F107,...
%             ApDaily,apNow,ap3h,ap6h,ap9h,apAver12To33h,apAver36To57h);
% 
%             [~,~,~,~,~,morningMsisNoAp(k),~,~,~,~,~]...
%             =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),morningLongitude,morningSolarTime,F107A,F107,3);
%         
%             [~,~,~,~,~,eveningMsis270km(k),~,~,~,~,~]...
%             =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),eveningLongitude,eveningSolarTime,F107A,F107,...
%             ApDaily,apNow,ap3h,ap6h,ap9h,apAver12To33h,apAver36To57h);
% 
%             [~,~,~,~,~,eveningMsisNoAp(k),~,~,~,~,~]...
%             =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),eveningLongitude,eveningSolarTime,F107A,F107,3);            
%         end
%         
%         morningResidue = mean(morningMsis270km - morningMsisNoAp);
%         eveningResidue = mean(eveningMsis270km - eveningMsisNoAp);
%         morningDensities(j,i) = morningResidue * 1e14;        
%         eveningDensities(j,i) = eveningResidue * 1e14;
%     end
%     p.progress;    
% end
% 
% p.stop;

morningDensities = morningDensities';
eveningDensities = eveningDensities';
morningDensities = bsxfun(@minus, morningDensities, mean(morningDensities));
eveningDensities = bsxfun(@minus, eveningDensities, mean(eveningDensities));
morningDensities = bsxfun(@rdivide, morningDensities, std(morningDensities));
eveningDensities = bsxfun(@rdivide, eveningDensities, std(eveningDensities));
morningDensities = bsxfun(@minus, morningDensities, mean(morningDensities,2));
eveningDensities = bsxfun(@minus, eveningDensities, mean(eveningDensities,2));
% h = 1/9*ones(3);
% morningDensities = filter2(h,morningDensities);
% eveningDensities = filter2(h,eveningDensities);
% morningDensities = morningDensities - min(morningDensities(:));
% morningDensities = morningDensities / max(morningDensities(:));
% eveningDensities = eveningDensities - min(eveningDensities(:));
% eveningDensities = eveningDensities / max(eveningDensities(:));

x = [doy1day - 365; doy1day; doy1day + 365];
fitCoeffs(10) = 2 * pi / 365;
parfor i = 1:length(morningBins)
    y = repmat(morningDensities(i,:)', 3, 1);
    fourierFit = fit(x, y, 'fourier4', 'StartPoint', fitCoeffs);
    morningFourierGrid{i} = fourierFit;
    morningPlot(i,:) = fourierFit(doy1day');
end

parfor i = 1:length(eveningBins)
    y = repmat(eveningDensities(i,:)', 3, 1);
    fourierFit = fit(x, y, 'fourier4', 'StartPoint', fitCoeffs);
    eveningFourierGrid{i} = fourierFit;
    eveningPlot(i,:) = fourierFit(doy1day');
end

[morningDoyGrid, morningLatGrid] = meshgrid(1:365, morningBins);
[eveningDoyGrid, eveningLatGrid] = meshgrid(1:365, eveningBins);
%figure('units','normalized','outerposition',[0 0 1 1]);
figure;
subplot(2,1,1)
surf(morningDoyGrid, morningLatGrid, morningPlot)
title('Morning Fourier Grid')
view(2);
shading flat
colorbar;
subplot(2,1,2)
surf(eveningDoyGrid, eveningLatGrid, eveningPlot)
title('Evening Fourier Grid')
view(2);
shading flat
colorbar;

figure;
subplot(2,1,1)
surf(morningDoyGrid, morningLatGrid, morningDensities)
title('Morning Residue Grid')
view(2);
shading flat
colorbar;
subplot(2,1,2)
surf(eveningDoyGrid, eveningLatGrid, eveningDensities)
title('Evening Residue Grid')
view(2);
shading flat
colorbar;

morningLatDoy = interp2(morningDoyGrid, morningLatGrid, morningDensities, morningDoy, morningLatitude, 'spline');
eveningLatDoy = interp2(eveningDoyGrid, eveningLatGrid, eveningDensities, eveningDoy, eveningLatitude, 'spline');

end

function [magneticLocalTime] = computeMagneticTime(magneticLongitude, doy, timestampsDatenum)
%

hours = (timestampsDatenum - floor(timestampsDatenum)) * 24;

subSolarLat = 23.5 * sin(0.0172 * doy - 1.405);
subSolarLon = 15 * (12 - hours);

altitude = ones(size(doy)) * 270e3;
[~, magneticSubSolarLon] = convertToMagneticCoordinates(subSolarLat, subSolarLon, altitude);

magneticLocalTime = 12 + (magneticLongitude - magneticSubSolarLon) / 15;
tooBigTimes = magneticLocalTime >= 24;
tooSmallTimes = magneticLocalTime < 0;
magneticLocalTime(tooBigTimes) = magneticLocalTime(tooBigTimes) - 24;
magneticLocalTime(tooSmallTimes) = magneticLocalTime(tooSmallTimes) + 24;

end

function [magneticLatitude, magneticLongitude] = convertToMagneticCoordinates(latitude, longitude, altitude)
% [magneticLatitude] = convertToMagneticCoordinates(latitude, longitude, altitude)

ecefXYZ = geod2ecef(latitude, longitude, altitude)';
ecefToMagTransform = [0.33907 -0.91964 -0.19826; ...
                      0.93826  0.34594  0      ; ...
                      0.06859  0.18602  0.98015];
magXYZ = ecefToMagTransform * ecefXYZ;
r = sqrt(magXYZ(1,:).^2 + magXYZ(2,:).^2 + magXYZ(3,:).^2);
magneticLatitude = pi/2 - acos(magXYZ(3,:) ./ r);
magneticLatitude = magneticLatitude' * 180 / pi;

magneticLongitude = 180 * atan2(magXYZ(2,:), magXYZ(1,:))' / pi;

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


