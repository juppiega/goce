function readFiles()
%

% poolobj = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(poolobj)
%     parpool();
% end

if(matlabpool('size')==0)
    matlabpool;
end


tic;
[ae, timestampsAeDatenum] = readAeFiles();

[F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, F10datenum, datenumToJulian] = readSolarIndexFiles();

[absB, timestampsAbsBDatenum, akasofuEpsilon, epsilonQualityFlag, timestampsEpsilonDatenum, vBz] ...
    = readSolarWindFiles(timestampsAeDatenum);

[density, longitude, latitude, altitude, solarTime, magneticLatitude, magneticLongitude, magneticLocalTime, crwindEast, crwindNorth, crwindUp, ...
 timestamps10sFixed, timestamps1min, timestamps1minFixed, timestampsDensityDatenum, doy, timestampsAbsB, timestampsEpsilon, firstDatenum, absDensityError, ...
 absWindError, noiseAffected, eclipseAffected, isMorningPass, ionThrusterActive] = readDensityFile(timestampsAeDatenum, timestampsAbsBDatenum, timestampsEpsilonDatenum);

goceData = struct('density', density, 'timestamps', timestampsDensityDatenum, 'latitude', latitude, 'longitude', longitude,...
                  'altitude', altitude, 'solarTime', solarTime);

[champData, graceData, de2Data, aeData] = computeOtherDensities();

[tiegcmDensityVariableAlt, tiegcmDensity270km] = readTiegcmFile(timestampsDensityDatenum);

[aeIntegrals, goceData, champData, graceData, de2Data, aeData] = computeAeIntegrals(ae, timestampsAeDatenum, timestampsDensityDatenum, goceData, champData, graceData, de2Data, aeData);

[apGoce, ap, timestampsAp, amGoce, dtc] = readApAndDtcFiles(timestampsDensityDatenum);

[F10, F81A, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, goceData, champData, graceData, de2Data, aeData] = giveSolarInputForModels(F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, F10datenum, timestampsDensityDatenum, goceData, champData, graceData, de2Data, aeData);

[apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, ApDaily, am3h, amAver24h, ap, timestamps3h, timestamps3hFixed, goceData, champData, graceData, de2Data, aeData] =...
    giveApValsForMSIS(apGoce, amGoce, ap, timestampsAp, timestamps10sFixed, timestamps1minFixed, timestamps1min, timestampsDensityDatenum, goceData, champData, graceData, de2Data, aeData);

[densityNoBg, msisDensityVariableAlt, msisDensity270km, msisDensity270kmNoAp, jb2008DensityVariableAlt, jb2008Density270km, jb2008Density270kmNoDtc, ...
    dtm2013Density270km, dtm2013DensityVariableAlt, dtm2013Density270kmNoAm, hwmU, hwmV, densityIndex, densityIndex1min, densityIndexNoBg, averagedDensity, averagedDensityNoBg, density3h, goceData, champData, graceData, de2Data, aeData]  = ...
    relateMsisToDensity(density, altitude, datenumToJulian, timestampsDensityDatenum, doy,...
    timestampsAeDatenum, timestamps1minFixed, timestamps3hFixed, solarTime, latitude, longitude, F10, F81A, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, am3h, amAver24h, dtc, ...
    ApDaily, apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, goceData, champData, graceData, de2Data, aeData);
measuredDensity = density;

[morningTimestamps10s, eveningTimestamps10s, morningLatitude, eveningLatitude, morningDoy, eveningDoy] = ...
    splitBySolarTime(timestamps10sFixed, latitude, doy, solarTime);

[morningGrid, morningBins] = computeTimeCells(morningTimestamps10s, firstDatenum, morningLatitude);
[eveningGrid, eveningBins] = computeTimeCells(eveningTimestamps10s, firstDatenum, eveningLatitude);

%msisSimulatedResidue = msisDensity270km - msisDensity270kmNoAp;

[morningFourierGrid, eveningFourierGrid, morningLatDoy, eveningLatDoy] = computeFourierGrids(morningGrid, morningBins, eveningGrid, ...
    eveningBins, doy, morningLatitude, eveningLatitude, morningDoy, eveningDoy, timestamps10sFixed, densityIndex, aeIntegrals(:,9), median(F10));

[hwmU, hwmV] = computeCrossTrackWind(hwmU, hwmV, crwindEast, crwindNorth);

if fclose('all') ~= 0
    display('File close unsuccesful - check if some of the files are reserved by another editor.')
end

fprintf('%s\n', 'Saving results to "goceVariables.mat" file')

save('goceVariables.mat', 'ae', '-v7')
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
save('goceVariables.mat', 'tiegcmDensity270km', '-append')
save('goceVariables.mat', 'tiegcmDensityVariableAlt', '-append')
save('goceVariables.mat', 'hwmU', '-append')
save('goceVariables.mat', 'hwmV', '-append')
save('goceVariables.mat', 'goceData', '-append')
save('goceVariables.mat', 'champData', '-append')
save('goceVariables.mat', 'graceData', '-append')
save('goceVariables.mat', 'de2Data', '-append')
save('goceVariables.mat', 'aeData', '-append')
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

aeFiles = dir('ae2*');
parfor i = 1:length(aeFiles)
    aeFile = fopen(aeFiles(i).name);
    if aeFile == -1
        error('aeFile open unsuccesful')
    end

    aeData = textscan(aeFile, '%s %s %f %f %f %f %f', 'MultipleDelimsAsOne',1, 'HeaderLines',15);

    ae = [ae; aeData{4}];
    strVec = strcat(aeData{1}, aeData{2});
    timestampsAeDatenum = [timestampsAeDatenum; datenum(strcat(aeData{1}, aeData{2}), 'yyyy-mm-ddHH:MM:SS.FFF')];
end
[timestampsAeDatenum, indicesToConserve, ~] = unique(timestampsAeDatenum);
ae = ae(indicesToConserve);

interpIndices = ae > 50000;
aeInterp = ae(~interpIndices);
tInterp = timestampsAeDatenum(~interpIndices);
ae = interp1(tInterp, aeInterp, timestampsAeDatenum, 'linear', 'extrap');

end

function [apGoce, ap, timestampsAp, am, dtc] = readApAndDtcFiles(timestampsDensityDatenum)
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
apGoce = reshape(apValues(apRows, :)', numel(apValues(apRows, :)), 1);

timestampsAp = repmat(timestampsAp', 8, 1); 
timestampsAp = timestampsAp + repmat((3:3:24)'/24, 1, size(timestampsAp, 2));
timestampsAp = timestampsAp(:);
apValues = apValues'; ap = apValues(:);

firstHour = round(str2double(datestr(timestampsDensityDatenum(1), 'HH')));
apGoce(1:firstHour/3) = [];


amFile = fopen('am_file_spider.dat');
if amFile == -1
    error('am File open unsuccesful! Check that you have a file "am_file_spider.data" in your WORKING DIRECTORY.')
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

function [F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, F10datenum, datenumToJulian] = readSolarIndexFiles
%
fprintf('%s\n', 'Began reading solar index files')

solarFile = fopen('SOLFSMY.TXT');
if solarFile == -1
    error('Solar index file open unsuccesful! Check that you have a file "SOLFSMY.TXT" in your WORKING DIRECTORY.')
end

solarData = textscan(solarFile, '%d %d %f %f %f %f %f %f %f %f %f %s', 'MultipleDelimsAsOne',1, 'CommentStyle','#');
%F10 = solarData{4};
%F81A = solarData{5};
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
F10 = solarData{5};
F81A = solarData{6};
F30toF10 = mean(F10 ./ F30);

year = solarData{1};
month = solarData{2};
day = solarData{3};
F10datenum = datenum(year, month, day);
F30 = F30toF10 * F30; %interp1(F10datenum, F30, indexDatenums, 'nearest', 'extrap');
F30A = F30toF10 * F30A;% interp1(F10datenum, F30A, indexDatenums, 'nearest', 'extrap');

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
fprintf('%s\n', 'Began reading GOCE files')

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

function [champData, graceData, de2Data, aeData] = computeOtherDensities()

fprintf('%s\n', 'Began reading CHAMP files')

champFiles1 = dir('champ/Density_3deg_0*');
champFiles2 = dir('champ/Density_10_*');

champDensity = [];
champTimestamps = [];
champLatitude = [];
champLongitude = [];
champLst = [];
champAltitude = [];
for i = 1:length(champFiles1)
    load(['champ/', champFiles1(i).name])
    yearVec = repmat([Year.data+2000,1,1,0,0,0], length(Sec.data), 1);
    thisFileTimestamps = datenum(yearVec) + Doy.data + Sec.data/86400 - 1;
    champTimestamps = [champTimestamps; thisFileTimestamps];
    champDensity = [champDensity; Density.data];
    champLatitude = [champLatitude; Lat.data];
    champLongitude = [champLongitude; Lon.data];
    champAltitude = [champAltitude; Height.data];
    champLst = [champLst; LocTim.data];
end
for i = 1:length(champFiles2)
    load(['champ/',champFiles2(i).name])
    yearVec = [loop.GPSyy'+2000, repmat([1,1,0,0,0], length(loop.GPSsec), 1)];
    thisFileTimestamps = datenum(yearVec) + loop.GPSdoy' + loop.GPSsec'/86400 - 1;
    champTimestamps = [champTimestamps; thisFileTimestamps];
    champDensity = [champDensity; loop.Dstar'];
    champLatitude = [champLatitude; loop.lat'];
    champLongitude = [champLongitude; loop.lon'];
    champAltitude = [champAltitude; loop.height'];
    champLst = [champLst; loop.slt'];
end
[champTimestamps, order] = unique(champTimestamps);
champDensity = champDensity(order);
champLatitude = champLatitude(order);
champLongitude = champLongitude(order);
champAltitude = champAltitude(order);
champLst = champLst(order);
champData = struct('density', champDensity, 'timestamps', champTimestamps, 'latitude', champLatitude, 'longitude', champLongitude,...
                  'altitude', champAltitude, 'solarTime', champLst);

fprintf('%s\n', 'Began reading GRACE files')
graceFiles1 = dir('grace/Density_graceA_3deg_0*');
graceFiles2 = dir('grace/Density_graceA_10_*');

graceDensity = [];
graceTimestamps = [];
graceLatitude = [];
graceLongitude = [];
graceLst = [];
graceAltitude = [];
for i = 1:length(graceFiles1)
    load(['grace/',graceFiles1(i).name])
    yearVec = repmat([Year.data+2000,1,1,0,0,0], length(Sec.data), 1);
    thisFileTimestamps = datenum(yearVec) + Doy.data + Sec.data/86400 - 1;
    graceTimestamps = [graceTimestamps; thisFileTimestamps];
    graceDensity = [graceDensity; Density.data];
    graceLatitude = [graceLatitude; Lat.data];
    graceLongitude = [graceLongitude; Lon.data];
    graceAltitude = [graceAltitude; Height.data];
    graceLst = [graceLst; LocTim.data];
end
for i = 1:length(graceFiles2)
    try
        load(['grace/',graceFiles2(i).name])
    catch
        continue
    end
    yearVec = [loop.GPSyy'+2000, repmat([1,1,0,0,0], length(loop.GPSsec'), 1)];
    thisFileTimestamps = datenum(yearVec) + loop.GPSdoy' + loop.GPSsec'/86400 - 1;
    graceTimestamps = [graceTimestamps; thisFileTimestamps];
    graceDensity = [graceDensity; loop.Dstar'];
    graceLatitude = [graceLatitude; loop.lat'];
    graceLongitude = [graceLongitude; loop.lon'];
    graceAltitude = [graceAltitude; loop.height'];
    graceLst = [graceLst; loop.slt'];
end
[graceTimestamps, order] = unique(graceTimestamps);
graceDensity = graceDensity(order);
graceLatitude = graceLatitude(order);
graceLongitude = graceLongitude(order);
graceAltitude = graceAltitude(order);
graceLst = graceLst(order);
graceData = struct('density', graceDensity, 'timestamps', graceTimestamps, 'latitude', graceLatitude, 'longitude', graceLongitude,...
                  'altitude', graceAltitude, 'solarTime', graceLst);

fprintf('%s\n', 'Began reading DE-2 files')       
de2Files = dir('ae_de2_3*');
de2Density = [];
de2Timestamps = [];
de2Latitude = [];
de2Longitude = [];
de2Lst = [];
de2Altitude = [];
for i = 1:length(de2Files)
    de2File = fopen(de2Files(i).name);
    data = textscan(de2File, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f','MultipleDelimsAsOne',1);
    thisFileTimestamps = datenum([data{1}+1900, data{2}, data{3}, data{4}, data{5}, data{6}]);
    de2Timestamps = [de2Timestamps; thisFileTimestamps];
    de2Altitude = [de2Altitude; data{7}];
    de2Latitude = [de2Latitude; data{8}];
    de2Longitude = [de2Longitude; data{9}];
    utHour = data{4} + data{5}/60 + data{6}/3600;
    lst = utHour + data{9} / 15;
    lst(lst >= 24) = lst(lst >= 24) - 24;
    lst(lst < 0) = lst(lst < 0) + 24;
    de2Lst = [de2Lst; lst];
    
    N2 = data{10}; N2(N2 > 1E30) = 0;
    O = data{11}; O(O > 1E30) = 0;
    He = data{12}; He(He > 1E30) = 0;
    Ar = data{13}; Ar(Ar > 1E30) = 0;
    N = data{14}; N(N > 1E30) = 0;
    thisFileDens = (28*N2 + 16*O + 4*He + 40*Ar + 14*N) * 1.661E-21; % kg/m3
    lowAlt = data{7} < 300; highAlt = ~lowAlt;
    acceptedVals = false(length(O),1);
    acceptedVals(lowAlt) = O(lowAlt) > 0 & N2(lowAlt) > 0;
    acceptedVals(highAlt) = O(highAlt) > 0 & He(highAlt) > 0;
    thisFileDens(~acceptedVals) = 0;
    de2Density = [de2Density; thisFileDens];
end
[de2Timestamps, order] = unique(de2Timestamps);
de2Density = de2Density(order);
de2Latitude = de2Latitude(order);
de2Longitude = de2Longitude(order);
de2Altitude = de2Altitude(order);
de2Lst = de2Lst(order);

nonZeroTimes = de2Density > 0;
de2Timestamps = de2Timestamps(nonZeroTimes);
de2Density = de2Density(nonZeroTimes);
de2Latitude = de2Latitude(nonZeroTimes);
de2Longitude = de2Longitude(nonZeroTimes);
de2Altitude = de2Altitude(nonZeroTimes);
de2Lst = de2Lst(nonZeroTimes);
de2Data = struct('density', de2Density, 'timestamps', de2Timestamps, 'latitude', de2Latitude, 'longitude', de2Longitude,...
                  'altitude', de2Altitude, 'solarTime', de2Lst);

fprintf('%s\n', 'Began reading AE-C/E files')
aeFiles = [dir('ae_c_*'); dir('ae_e_2*')];
aeDensity = [];
fromNace = [];
fromOss = [];
aeTimestamps = [];
aeLatitude = [];
aeLongitude = [];
aeLst = [];
aeAltitude = [];
for i = 1:length(aeFiles)
    aeFile = fopen(aeFiles(i).name);
    data = textscan(aeFile, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','MultipleDelimsAsOne',1);
    thisFileTimestamps = datenum([data{1}+1900, data{2}, data{3}, data{4}, data{5}, data{6}]);
    utHour = data{4} + data{5}/60 + data{6}/3600;
    lst = utHour + data{9} / 15;
    lst(lst >= 24) = lst(lst >= 24) - 24;
    lst(lst < 0) = lst(lst < 0) + 24;
    
    aeLst = [aeLst; lst];
    aeTimestamps = [aeTimestamps; thisFileTimestamps];
    aeAltitude = [aeAltitude; data{7}];
    aeLatitude = [aeLatitude; data{8}];
    aeLongitude = [aeLongitude; data{9}];
    
    N2 = data{11}; N2(N2 > 1E30) = 0;
    O = data{12}; O(O > 1E30) = 0;
    He = data{13}; He(He > 1E30) = 0;
    Ar = data{14}; Ar(Ar > 1E30) = 0;
    NO = data{15}; NO(NO > 1E30) = 0;
    naceDens = (28*N2 + 16*O + 4*He + 40*Ar + 30*NO) * 1.661E-21; % kg/m3  
    lowAlt = data{7} < 300; highAlt = ~lowAlt;
    acceptedValsNace = false(length(O),1);
    acceptedValsNace(lowAlt) = O(lowAlt) > 0 & N2(lowAlt) > 0;
    acceptedValsNace(highAlt) = O(highAlt) > 0 & He(highAlt) > 0;
    
    naceInd = naceDens > 0;
    
    N2 = data{16}; N2(N2 > 1E30) = 0;
    O2 = data{17}; O2(O2 > 1E30) = 0;
    O = data{19}; O(O > 1E30) = 0;
    He = data{18}; He(He > 1E30) = 0;
    Ar = data{20}; Ar(Ar > 1E30) = 0;
    N = data{21}; N(N > 1E30) = 0;
    ossDens = (28*N2 + 32*O2 + 16*O + 4*He + 40*Ar + 14*N) * 1.661E-21; % kg/m3 
    acceptedValsOss = false(length(O),1);
    acceptedValsOss(lowAlt) = O(lowAlt) > 0 & N2(lowAlt) > 0;
    acceptedValsOss(highAlt) = O(highAlt) > 0 & He(highAlt) > 0;
    ossInd = ossDens > 0;
    
    thisFileDens = zeros(length(naceDens), 1);
    thisFileDens(naceInd) = naceDens(naceInd);
    thisFileDens(ossInd) = ossDens(ossInd); % OSS priority if both
    acceptedVals = acceptedValsNace | acceptedValsOss;
    thisFileDens(~acceptedVals) = 0;
    
    aeDensity = [aeDensity; thisFileDens];   
    fromNace = [fromNace; naceInd];
    fromOss = [fromOss; ossInd];
end
[aeTimestamps, order] = unique(aeTimestamps);
aeDensity = aeDensity(order);
aeLatitude = aeLatitude(order);
aeLongitude = aeLongitude(order);
aeAltitude = aeAltitude(order);
aeLst = aeLst(order);
fromNace = fromNace(order);
fromOss = fromOss(order);

nonZeroTimes = aeDensity > 0;
aeTimestamps = aeTimestamps(nonZeroTimes);
aeDensity = aeDensity(nonZeroTimes);
aeLatitude = aeLatitude(nonZeroTimes);
aeLongitude = aeLongitude(nonZeroTimes);
aeAltitude = aeAltitude(nonZeroTimes);
aeLst = aeLst(nonZeroTimes);
fromNace = fromNace(nonZeroTimes);
fromOss = fromOss(nonZeroTimes);
aeData = struct('density', aeDensity, 'timestamps', aeTimestamps, 'latitude', aeLatitude, 'longitude', aeLongitude,...
                  'altitude', aeAltitude, 'solarTime', aeLst, 'fromNace', fromNace, 'fromOss', fromOss);

end

function [tiegcmDensityVariableAlt, tiegcmDensity270km] = readTiegcmFile(timestampsDensityDatenum)

tiegcmDensityVariableAlt = ones(size(timestampsDensityDatenum));
tiegcmDensity270km = ones(size(timestampsDensityDatenum));

if exist('tiegcmDens.mat', 'file')
    load tiegcmDens.mat
    ind1 = ismember(timestampsDensityDatenum, tiegcmGoceDatenums);
    ind2 = ismember(tiegcmGoceDatenums, timestampsDensityDatenum);
    tiegcmDensityVariableAlt(ind1) = tiegcmGoceInterp(ind2) * 1e14;
    tiegcmDensity270km(ind1) = tiegcmGoce270km(ind2) * 1e14;
end

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

function [aeIntGoce, goceData, champData, graceData, de2Data, aeData] = computeAeIntegrals(ae, timestampsAe, timestampsGoce, goceData, champData, graceData, de2Data, aeData)
%
lags = [2 4 8 16 21 30 40 50 60];
aeIntGoce = zeros(length(timestampsGoce), length(lags));
aeIntChamp = zeros(length(champData.timestamps), length(lags));
aeIntGrace = zeros(length(graceData.timestamps), length(lags));
aeIntDe2 = zeros(length(de2Data.timestamps), length(lags));
aeIntAE = zeros(length(aeData.timestamps), length(lags));
cumulativeAe = cumsum(ae);

oneHour = 60;
for i = 1:length(lags)
    lag = lags(i) * oneHour;
    aeInt = cumulativeAe(lag + 1 : end) - cumulativeAe(1 : end - lag);
    aeTime = timestampsAe(lag + 1 : end);
    aeIntGoce(:,i) = interp1(aeTime, aeInt, timestampsGoce, 'linear', 0);
    aeIntChamp(:,i) = interp1(aeTime, aeInt, champData.timestamps, 'linear', 0);
    aeIntGrace(:,i) = interp1(aeTime, aeInt, graceData.timestamps, 'linear', 0);
    aeIntDe2(:,i) = interp1(aeTime, aeInt, de2Data.timestamps, 'linear', 0);
    aeIntAE(:,i) = interp1(aeTime, aeInt, aeData.timestamps, 'linear', 0);
end

goceData.aeInt = aeIntGoce;
champData.aeInt = aeIntChamp;
graceData.eInt = aeIntGrace;
de2Data.aeInt = aeIntDe2;
aeData.aeInt = aeIntAE;

end

function [F10out, F81Aout, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, goceData, champData, graceData, de2Data, aeData] = giveSolarInputForModels(F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, F10datenum, timestampsDensityDatenum, goceData, champData, graceData, de2Data, aeData)
%
fprintf('%s\n', 'Computing solar index values for models')

msisF81A = interp1(F10datenum, F81A, timestampsDensityDatenum, 'linear', 100);
F30A = interp1(F10datenum, F30A, timestampsDensityDatenum, 'linear', 100);

F10out = interp1(F10datenum + 1, F10, timestampsDensityDatenum, 'linear', 100);
F81Aout = interp1(F10datenum + 1, F81A, timestampsDensityDatenum, 'linear', 100);
S10 = interp1(indexDatenums + 1, S10, timestampsDensityDatenum, 'linear', 100);
S81A = interp1(indexDatenums + 1, S81A, timestampsDensityDatenum, 'linear', 100);
F30 = interp1(F10datenum + 1, F30, timestampsDensityDatenum, 'linear', 100);

M10 = interp1(indexDatenums + 2, M10, timestampsDensityDatenum, 'linear', 100);
M81A = interp1(indexDatenums + 2, M81A, timestampsDensityDatenum, 'linear', 100);

Y10 = interp1(indexDatenums + 5, Y10, timestampsDensityDatenum, 'linear', 100);
Y81A = interp1(indexDatenums + 5, Y81A, timestampsDensityDatenum, 'linear', 100);

goceData.F10 = interp1(F10datenum + 1, F10, timestampsDensityDatenum, 'previous', 100);
goceData.F81A = interp1(F10datenum, F81A, timestampsDensityDatenum, 'previous', 100);

champData.F10 = interp1(F10datenum + 1, F10, champData.timestamps, 'previous', 100);
champData.F81A = interp1(F10datenum, F81A, champData.timestamps, 'previous', 100);

graceData.F10 = interp1(F10datenum + 1, F10, graceData.timestamps, 'previous', 100);
graceData.F81A = interp1(F10datenum, F81A, graceData.timestamps, 'previous', 100);

de2Data.F10 = interp1(F10datenum + 1, F10, de2Data.timestamps, 'previous', 100);
de2Data.F81A = interp1(F10datenum, F81A, de2Data.timestamps, 'previous', 100);

aeData.F10 = interp1(F10datenum + 1, F10, aeData.timestamps, 'previous', 100);
aeData.F81A = interp1(F10datenum, F81A, aeData.timestamps, 'previous', 100);

end

function [apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, ApDaily, am3h, amAver24h, apAllFixed, timestamps3h, timestamps3hFixed, goceData, champData, graceData, de2Data, aeData] = ...
    giveApValsForMSIS(apGoce, amGoce, ap, timestampsAp, timestamps10sFixed, timestamps1minFixed, timestamps1min, timestampsDensityDatenum, goceData, champData, graceData, de2Data, aeData)
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

ap10s = reshape(repmat(apGoce', 180 * 6, 1),[],1);
apTime10s = reshape(repmat(timestamps3h, 180 * 6, 1),[],1);
ap10s(1) = [];
apTime10s(1) = [];
apNow = ap10s(ismember(apTime10s, nextThreeHStamp));
ap3h = ap10s(ismember(apTime10s, nextThreeHStamp - threeHinSec));
ap6h = ap10s(ismember(apTime10s, nextThreeHStamp - 2 * threeHinSec));
ap9h = ap10s(ismember(apTime10s, nextThreeHStamp - 3 * threeHinSec));
ap10sSmoothed = reshape(repmat(smooth(apGoce, 9)', 180 * 6, 1),[],1);
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

am3h = interp1(timestamps3h + threeHinSec, amGoce, timestamps10sFixed, 'linear', 'extrap');

averagingMatrix = nan(8, length(amGoce) - 7);

for i = 1:8
    lastIndex = length(amGoce) - 8 + i;
    averagingMatrix(i,:) = amGoce(i:lastIndex)';
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
apAllFixed = apGoce(ismember(timestamps3h, timestamps1min));
timestamps3hFixed = timestamps3h(ismember(timestamps3h, timestamps1minFixed))';
timestamps3h = timestamps3h(ismember(timestamps3h, timestamps1min))';

dt = 3 / 24;
ap24hSmooth = smooth(ap, 8);
champData.apNow = interp1(timestampsAp, ap, champData.timestamps, 'previous', 'extrap');
champData.ap3h = interp1(timestampsAp + dt, ap, champData.timestamps, 'previous', 'extrap');
champData.ap6h = interp1(timestampsAp + 2*dt, ap, champData.timestamps, 'previous', 'extrap');
champData.ap9h = interp1(timestampsAp + 3*dt, ap, champData.timestamps, 'previous', 'extrap');
champData.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, champData.timestamps, 'previous', 'extrap');
champData.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, champData.timestamps, 'previous', 'extrap');
champData.ApDaily = interp1(timestampsAp, ap24hSmooth, champData.timestamps, 'previous', 'extrap');

graceData.apNow = interp1(timestampsAp, ap, graceData.timestamps, 'previous', 'extrap');
graceData.ap3h = interp1(timestampsAp + dt, ap, graceData.timestamps, 'previous', 'extrap');
graceData.ap6h = interp1(timestampsAp + 2*dt, ap, graceData.timestamps, 'previous', 'extrap');
graceData.ap9h = interp1(timestampsAp + 3*dt, ap, graceData.timestamps, 'previous', 'extrap');
graceData.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, graceData.timestamps, 'previous', 'extrap');
graceData.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, graceData.timestamps, 'previous', 'extrap');
graceData.ApDaily = interp1(timestampsAp, ap24hSmooth, graceData.timestamps, 'previous', 'extrap');

goceData.apNow = interp1(timestampsAp, ap, goceData.timestamps, 'previous', 'extrap');
goceData.ap3h = interp1(timestampsAp + dt, ap, goceData.timestamps, 'previous', 'extrap');
goceData.ap6h = interp1(timestampsAp + 2*dt, ap, goceData.timestamps, 'previous', 'extrap');
goceData.ap9h = interp1(timestampsAp + 3*dt, ap, goceData.timestamps, 'previous', 'extrap');
goceData.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, goceData.timestamps, 'previous', 'extrap');
goceData.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, goceData.timestamps, 'previous', 'extrap');
goceData.ApDaily = interp1(timestampsAp, ap24hSmooth, goceData.timestamps, 'previous', 'extrap');

de2Data.apNow = interp1(timestampsAp, ap, de2Data.timestamps, 'previous', 'extrap');
de2Data.ap3h = interp1(timestampsAp + dt, ap, de2Data.timestamps, 'previous', 'extrap');
de2Data.ap6h = interp1(timestampsAp + 2*dt, ap, de2Data.timestamps, 'previous', 'extrap');
de2Data.ap9h = interp1(timestampsAp + 3*dt, ap, de2Data.timestamps, 'previous', 'extrap');
de2Data.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, de2Data.timestamps, 'previous', 'extrap');
de2Data.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, de2Data.timestamps, 'previous', 'extrap');
de2Data.ApDaily = interp1(timestampsAp, ap24hSmooth, de2Data.timestamps, 'previous', 'extrap');

aeData.apNow = interp1(timestampsAp, ap, aeData.timestamps, 'previous', 'extrap');
aeData.ap3h = interp1(timestampsAp + dt, ap, aeData.timestamps, 'previous', 'extrap');
aeData.ap6h = interp1(timestampsAp + 2*dt, ap, aeData.timestamps, 'previous', 'extrap');
aeData.ap9h = interp1(timestampsAp + 3*dt, ap, aeData.timestamps, 'previous', 'extrap');
aeData.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, aeData.timestamps, 'previous', 'extrap');
aeData.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, aeData.timestamps, 'previous', 'extrap');
aeData.ApDaily = interp1(timestampsAp, ap24hSmooth, aeData.timestamps, 'previous', 'extrap');

end

function [correctedDensity, msisDensityVariableAlt, msisDensity270km, msisDensity270kmNoAp, jb2008DensityVariableAlt, jb2008Density270km, jb2008Density270kmNoDtc, dtm2013Density270km, dtm2013DensityVariableAlt, dtm2013Density270kmNoAm, hwmU, hwmV, densityIndex, densityIndex1min, densityIndexNoBg, averagedDensity, averagedDensityNoBg, density3h, goceData, champData, graceData, de2Data, aeData] = relateMsisToDensity(density, altitude, datenumToJulian, timestampsDatenum, doy,...
          timestampsAeDatenum, timestamps1minFixed, timestamps3hFixed, solarTime, latitude, longitude, F10, F81A, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, am3h, amAver24h, dtc, ApDaily, apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, goceData, champData, graceData, de2Data, aeData)
      
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

obsYear = datevec(timestampsDatenum);
obsYear = obsYear(:,1);

goceData = computeRelativeChanges(goceData);
champData = computeRelativeChanges(champData);
graceData = computeRelativeChanges(graceData);
de2Data = computeRelativeChanges(de2Data);
aeData = computeRelativeChanges(aeData);

targetCount = round(length(modelingIndices) / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running MSIS, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

parfor i = modelingIndices
    [~,~,~,~,~,msisDensityVariableAlt(i),~,~,~,~,~]...
      =nrlmsise_mex(doy(i),seconds(i),altitudeInKm(i),latitude(i),longitude(i),solarTime(i),msisF81A(i),F10(i),...
      ApDaily(i),apNow(i),ap3h(i),ap6h(i),ap9h(i),apAver12To33h(i),apAver36To57h(i));

    [~,~,~,~,~,msisDensity270km(i),~,~,~,~,~]...
      =nrlmsise_mex(doy(i),seconds(i),270,latitude(i),longitude(i),solarTime(i),msisF81A(i),F10(i),...
      ApDaily(i),apNow(i),ap3h(i),ap6h(i),ap9h(i),apAver12To33h(i),apAver36To57h(i));

    [~,~,~,~,~,msisDensity270kmNoAp(i),~,~,~,~,~]...
      =nrlmsise_mex(doy(i),seconds(i),270,latitude(i),longitude(i),solarTime(i),msisF81A(i),F10(i), 3);

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
    [~,~,jb2008DensityVariableAlt(i)] = jb2008_mex(julianDay(i), altitudeInKm(i), latitude(i), longitude(i), F10(i), F81A(i), S10(i),...
        S81A(i), M10(i), M81A(i), Y10(i), Y81A(i), dtc(i));
    
    [~,~,jb2008Density270km(i)] = jb2008_mex(julianDay(i), 270, latitude(i), longitude(i), F10(i), F81A(i), S10(i),...
        S81A(i), M10(i), M81A(i), Y10(i), Y81A(i), dtc(i));
    
    [~,~,jb2008Density270kmNoDtc(i)] = jb2008_mex(julianDay(i), 270, latitude(i), longitude(i), F10(i), F81A(i), S10(i),...
        S81A(i), M10(i), M81A(i), Y10(i), Y81A(i), 0);
    
    if mod(i, 10000) == 0
     p.progress;
    end
end
p.stop;

dtm2013Density270km = nan(size(density));
dtm2013DensityVariableAlt = nan(size(density));
dtm2013Density270kmNoAm = nan(size(density));
dtm2013_mex();

p = TimedProgressBar( targetCount, barWidth, ...
                    'Running DTM-2013, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
                
for i = modelingIndices
    [~, ~, dtm2013DensityVariableAlt(i),~,~,~] = dtm2013_mex(doyDecimal(i), altitudeInKm(i), latitude(i), longitude(i), ...
        solarTime(i), F30(i), F30A(i), am3h(i), amAver24h(i));
    
    [~, ~, dtm2013Density270km(i),~,~,~] = dtm2013_mex(doyDecimal(i), 270, latitude(i), longitude(i), ...
        solarTime(i), F30(i), F30A(i), am3h(i), amAver24h(i));
    
    [~, ~, dtm2013Density270kmNoAm(i),~,~,~] = dtm2013_mex(doyDecimal(i), 270, latitude(i), longitude(i), ...
        solarTime(i), F30(i), F30A(i), 0.0, 0.0);
    
    if mod(i, 10000) == 0
     p.progress;
    end
end
p.stop;

hwmU = ones(size(density));
hwmV = ones(size(density));

% p = TimedProgressBar( targetCount, barWidth, ...
%                     'Running HWM07, ETA ', ...
%                     '. Now at ', ...
%                     'Completed in ' );
%                 
% parfor i = modelingIndices
%     [hwmU(i), hwmV(i)] = hwm07_mex(obsYear(i), doyDecimal(i), altitudeInKm(i), latitude(i), longitude(i), apNow(i));
%     
%     if mod(i, 10000) == 0
%      p.progress;
%     end
% end
% p.stop;

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

function satelliteData = computeRelativeChanges(satelliteData)

t = satelliteData.timestamps;
seconds = (t - floor(t)) * 86400;
[tYear,~,~,~,~,~] = datevec(t);
yearVec = [tYear, ones(length(t),2), zeros(length(t), 3)];
doy = t - datenum(yearVec) + 1;
altitudeInKm = satelliteData.altitude;
latitude = satelliteData.latitude;
longitude = satelliteData.longitude;
solarTime = satelliteData.solarTime;
F81A = satelliteData.F81A;
F10 = satelliteData.F10;
ApDaily = satelliteData.ApDaily;
apNow = satelliteData.apNow;
ap3h = satelliteData.ap3h;
ap6h = satelliteData.ap3h;
ap9h = satelliteData.ap3h;
apAver12To33h = satelliteData.apAver12To33h;
apAver36To57h = satelliteData.apAver36To57h;
density = satelliteData.density;

modelingIndices = 1:length(t);

targetCount = round(length(modelingIndices) / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Computing relative change, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

parfor i = modelingIndices
    [~,~,~,~,~,msisDensityVariableAlt(i),~,~,~,~,~]...
      =nrlmsise_mex(doy(i),seconds(i),altitudeInKm(i),latitude(i),longitude(i),solarTime(i),F81A(i),F10(i),...
      ApDaily(i),apNow(i),ap3h(i),ap6h(i),ap9h(i),apAver12To33h(i),apAver36To57h(i));

    [~,~,~,~,~,msisDensityNoAp(i),~,~,~,~,~]...
      =nrlmsise_mex(doy(i),seconds(i),altitudeInKm(i),latitude(i),longitude(i),solarTime(i),F81A(i),F10(i), 3);

    if mod(i, 10000) == 0
     p.progress;
    end
end
p.stop;

msisDensityVariableAlt = msisDensityVariableAlt';
msisDensityNoAp = msisDensityNoAp';

density = mean(msisDensityVariableAlt ./ density) * density;
relativeChange = density ./ msisDensityNoAp;
satelliteData.relativeChange = relativeChange;

end

function [crwindEast, crwindNorth] = computeCrossTrackWind(fullU, fullV, goceCrwindEast, goceCrwindNorth)

goceMag = sqrt(goceCrwindEast.^2 + goceCrwindNorth.^2);
crossTrackX = goceCrwindEast ./ goceMag;
crossTrackY = goceCrwindNorth ./ goceMag;

projectionMag = fullU .* crossTrackX + fullV .* crossTrackY;

crwindEast = projectionMag .* crossTrackX;
crwindNorth = projectionMag .* crossTrackY;

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


