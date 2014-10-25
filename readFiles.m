function readFiles()
%

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(16);
end

tic;
[ae, timestampsAeDatenum] = readAeFiles();

[F107values, F107datenum] = readF107File();

[absB, timestampsAbsBDatenum, akasofuEpsilon, epsilonQualityFlag, timestampsEpsilonDatenum, vBz] ...
    = readSolarParameterFiles(timestampsAeDatenum);

[density, longitude, latitude, altitude, solarTime, magneticLatitude, crwindEast, crwindNorth, crwindUp, ...
 timestamps10sFixed, timestamps1min, timestamps1minFixed, timestampsDensityDatenum, doy, timestampsAbsB, timestampsEpsilon, firstDatenum, absDensityError, ...
 absWindError, noiseAffected, eclipseAffected, isMorningPass, ionThrusterActive] = readDensityFile(timestampsAeDatenum, timestampsAbsBDatenum, timestampsEpsilonDatenum);

[aeIntegrals] = computeAeIntegrals(ae, timestamps1min, timestamps10sFixed, absB, timestampsAbsB);

apAll = readApFile(timestampsDensityDatenum);

[F107A81Days, F107Yesterday] = giveF107ValsForMSIS(F107values, F107datenum, timestampsDensityDatenum);

[apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, ApDaily, ap, timestamps3h, timestamps3hFixed] =...
    giveApValsForMSIS(apAll, timestamps10sFixed, timestamps1minFixed, timestamps1min, timestampsDensityDatenum);

[densityNoBg, msisDensityVariableAlt, msisDensity270km, msisDensity270kmNoAp, densityIndex, densityIndex1min, densityIndexNoBg, averagedDensity, averagedDensityNoBg, density3h] = relateMsisToDensity(density, altitude, timestampsDensityDatenum, doy,...
    timestampsAeDatenum, timestamps1minFixed, timestamps3hFixed, solarTime, latitude, longitude, F107A81Days, F107Yesterday, ApDaily, apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h);
measuredDensity = density;

[morningTimestamps10s, eveningTimestamps10s, morningLatitude, eveningLatitude] = ...
    splitBySolarTime(timestamps10sFixed, latitude, solarTime);

[morningGrid, morningBins] = computeTimeCells(morningTimestamps10s, firstDatenum, morningLatitude);
[eveningGrid, eveningBins] = computeTimeCells(eveningTimestamps10s, firstDatenum, eveningLatitude);

msisSimulatedResidue = msisDensity270km - msisDensity270kmNoAp;
[morningFourierGrid, eveningFourierGrid] = computeFourierGrids(morningGrid, morningBins, eveningGrid, eveningBins, doy, timestamps10sFixed, msisSimulatedResidue, aeIntegrals(:,9));

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
save('goceVariables.mat', 'densityNoBg', '-append')
save('goceVariables.mat', 'msisDensityVariableAlt', '-append')
save('goceVariables.mat', 'msisDensity270km', '-append')
save('goceVariables.mat', 'msisDensity270kmNoAp', '-append')
save('goceVariables.mat', 'morningFourierGrid', '-append')
save('goceVariables.mat', 'eveningFourierGrid', '-append')
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

function ap = readApFile(timestampsDensityDatenum)
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


end

function [F107values, F107datenum] = readF107File
%
fprintf('%s\n', 'Began reading F107 file')

F107File = fopen('F107data');
if F107File == -1
    error('F10.7 File open unsuccesful! Check that you have a file "F107data" in your WORKING DIRECTORY.')
end

F10data = textscan(F107File, '%f %f %f %f', 'MultipleDelimsAsOne',1);
F107values = F10data{4};
F10year = F10data{1};
F10dayOfYear = F10data{2};
F107datenum = datenum(F10year, 1, 1) + F10dayOfYear - 1;

indicesToConserve = find(F107values < 990);
tF107ToInterpolate = F107datenum(indicesToConserve);
F107valuesToInterpolate = F107values(indicesToConserve);
F107values = interp1(tF107ToInterpolate, F107valuesToInterpolate, F107datenum, 'linear', 'extrap');

end

function [absB, timestampsAbsBDatenum, akasofuEpsilon, epsilonQualityFlag, timestampsEpsilonDatenum, vBz] = ...
    readSolarParameterFiles(timestampsAeDatenum)
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
theta = atan(By ./ Bz);

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
timestampsAbsB = round((tBToInterpolate-tBToInterpolate(1))*secondsInDay + find(timestampsAeDatenum<tBToInterpolate(1), 1, 'last') * 60);
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

function [density, longitude, latitude, altitude, solarTime, magneticLatitude, crwindEast, crwindNorth, crwindUp,...
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

magneticLatitude = convertToMagneticCoordinates(latitude, longitude, altitude);

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

function [F107A81Days, F107Yesterday] = giveF107ValsForMSIS(F107values, F107datenum, timestampsDensityDatenum)
%
fprintf('%s\n', 'Computing F10.7 values for msis')

F107beginDay = floor(timestampsDensityDatenum(1)) - 81;
F107endDay = floor(timestampsDensityDatenum(end));
F107values = F107values(F107datenum >= F107beginDay & F107datenum <= F107endDay);
F107datenum = (F107beginDay:F107endDay)';
F107smoothed = smooth(F107values, 81);
yesterdayDatenums = floor(timestampsDensityDatenum) - 1;
uniqueDatenums = unique(yesterdayDatenums, 'stable');
F107Yesterday = zeros(length(yesterdayDatenums), 1);
F107A81Days = zeros(length(yesterdayDatenums), 1);
k = 1;
for i = 1:length(uniqueDatenums)
    numOfOccurences = length(find(yesterdayDatenums == uniqueDatenums(i)));
    F107Yesterday(k:numOfOccurences + k - 1) = F107values(F107datenum == uniqueDatenums(i));
    F107A81Days(k:numOfOccurences + k - 1) = F107smoothed(find(F107datenum == uniqueDatenums(i)) - 40);
    k = k + numOfOccurences;
end

end

function [apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, ApDaily, apAllFixed, timestamps3h, timestamps3hFixed] = ...
    giveApValsForMSIS(apAll, timestamps10sFixed, timestamps1minFixed, timestamps1min, timestampsDensityDatenum)
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

timestamps3h(timestamps3h <= firstHour * threeHinSec) = nan(1);
apAllFixed = apAll(ismember(timestamps3h, timestamps1min));
timestamps3hFixed = timestamps3h(ismember(timestamps3h, timestamps1minFixed))';
timestamps3h = timestamps3h(ismember(timestamps3h, timestamps1min))';

end

function [correctedDensity, msisDensityVariableAlt, msisDensity270km, msisDensity270kmNoAp, densityIndex, densityIndex1min, densityIndexNoBg, averagedDensity, averagedDensityNoBg, density3h] = relateMsisToDensity(density, altitude, timestampsDatenum, doy,...
          timestampsAeDatenum, timestamps1minFixed, timestamps3hFixed, solarTime, latitude, longtitude, F107A, F107, ApDaily, apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h)
      
fprintf('%s\n', 'Computing normalized densities with msis. This may take even half an hour.')

msisDensity270km = nan(size(density));
msisDensityVariableAlt = nan(size(density));
msisDensity270kmNoAp = nan(size(density));

modelingIndices = 1:length(density);

altitudeInKm = altitude / 1000;

doy = ceil(doy);
doy(doy == 0) = 1;

% second of the day
secondsInDay = 24 * 60 * 60;
seconds = (timestampsDatenum - floor(timestampsDatenum)) * secondsInDay;

targetCount = round(length(modelingIndices) / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running MSIS, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

parfor i = modelingIndices
      [~,~,~,~,~,msisDensityVariableAlt(i),~,~,~,~,~]...
          =nrlmsise_mex(doy(i),seconds(i),altitudeInKm(i),latitude(i),longtitude(i),solarTime(i),F107A(i),F107(i),...
          ApDaily(i),apNow(i),ap3h(i),ap6h(i),ap9h(i),apAver12To33h(i),apAver36To57h(i));
      
      [~,~,~,~,~,msisDensity270km(i),~,~,~,~,~]...
          =nrlmsise_mex(doy(i),seconds(i),270,latitude(i),longtitude(i),solarTime(i),F107A(i),F107(i),...
          ApDaily(i),apNow(i),ap3h(i),ap6h(i),ap9h(i),apAver12To33h(i),apAver36To57h(i));
      
      [~,~,~,~,~,msisDensity270kmNoAp(i),~,~,~,~,~]...
          =nrlmsise_mex(doy(i),seconds(i),270,latitude(i),longtitude(i),solarTime(i),F107A(i),F107(i),3,3,3,3,3,3,3);
      
      if mod(i, 10000) == 0
         p.progress;
      end
end
p.stop;

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

function [morningTimestamps10s, eveningTimestamps10s, morningLatitude, eveningLatitude] = ...
    splitBySolarTime(timestamps10s, latitude, solarTime)
%

morningIndices = find(solarTime <= 12);
eveningIndices = find(solarTime > 12);

morningTimestamps10s = timestamps10s(morningIndices);
morningLatitude = latitude(morningIndices);

eveningTimestamps10s = timestamps10s(eveningIndices);
eveningLatitude = latitude(eveningIndices);

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

function [morningFourierGrid, eveningFourierGrid] = computeFourierGrids(morningGrid, morningBins, eveningGrid, eveningBins, doy, timestamps10sFixed, densityIndex, aeIntegral)
%

fprintf('%s\n', 'Computing fourier fits.')

targetCount = length(morningBins);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running Parametrizations, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
%aeLimit = aeIntegral < quantile(aeIntegral, 0.99);

binWidth = 10;
doy1day = (1:365)';

% parfor i = 1:length(morningBins)
%     timeThisLat = morningGrid{i};
%     thisLatIndices = ismember(timestamps10sFixed, timeThisLat);
%     densityThisLat = densityIndex(thisLatIndices);
%     
%     doyThisLat = ceil(doy(thisLatIndices));
%     
%     thisLatDens = zeros(length(doy1day),1);
%     for k = 1:length(doy1day)
%         indices = doyThisLat == doy1day(k);
%         thisLatDens(k,:) = mean(densityThisLat(indices));
%     end
%    
%     morningDensities(:,i) = thisLatDens; 
%     
%     p.progress;
% end
% 
% parfor i = 1:length(eveningBins)
%     timeThisLat = eveningGrid{i};
%     thisLatIndices = ismember(timestamps10sFixed, timeThisLat);
%     densityThisLat = densityIndex(thisLatIndices);
% 
%     doyThisLat = ceil(doy(thisLatIndices));
%     
%     thisLatDens = zeros(length(doy1day),1);
%     for k = 1:length(doy1day)
%         indices = doyThisLat == doy1day(k);
%         thisLatDens(k,:) = mean(densityThisLat(indices));
%     end
% 
%     eveningDensities(:,i) = thisLatDens;
%     
%     p.progress;
% end

latitude = morningBins;
morningDensities = zeros(length(doy1day), length(latitude));
eveningDensities = zeros(length(doy1day), length(latitude));
morningSolarTime = 7.0;
eveningSolarTime = 19.0;
seconds = 0 : 3 * 60 * 60 : 24 * 60 * 60 - 1;
morningMsis270km = zeros(length(seconds), 1);
morningMsisNoAp = zeros(length(seconds), 1);
eveningMsis270km = zeros(length(seconds), 1);
eveningMsisNoAp = zeros(length(seconds), 1);

F107A = 100;
F107 = 100;
ApDaily = 39;
apNow = 56;
ap3h = 94;
ap6h = 80;
ap9h = 39;
apAver12To33h = 6;
apAver36To57h = 6;
for i = 1:length(latitude)
    for j = 1:length(doy1day)
        for k = 1:length(seconds)
            morningLongitude = 180 * (morningSolarTime - (seconds(k)/3600)) / 12;
            eveningLongitude = 180 * (eveningSolarTime - (seconds(k)/3600)) / 12;
            morningLongitude(morningLongitude < -180) = morningLongitude + 360;
            eveningLongitude(eveningLongitude > 180) = eveningLongitude - 360;
            
            [~,~,~,~,~,morningMsis270km(k),~,~,~,~,~]...
            =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),morningLongitude,morningSolarTime,F107A,F107,...
            ApDaily,apNow,ap3h,ap6h,ap9h,apAver12To33h,apAver36To57h);

            [~,~,~,~,~,morningMsisNoAp(k),~,~,~,~,~]...
            =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),morningLongitude,morningSolarTime,F107A,F107,3,3,3,3,3,3,3);
        
            [~,~,~,~,~,eveningMsis270km(k),~,~,~,~,~]...
            =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),eveningLongitude,eveningSolarTime,F107A,F107,...
            ApDaily,apNow,ap3h,ap6h,ap9h,apAver12To33h,apAver36To57h);

            [~,~,~,~,~,eveningMsisNoAp(k),~,~,~,~,~]...
            =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),eveningLongitude,eveningSolarTime,F107A,F107,3,3,3,3,3,3,3);            
        end
        
        morningResidue = mean(morningMsis270km - morningMsisNoAp);
        eveningResidue = mean(eveningMsis270km - eveningMsisNoAp);
        morningDensities(j,i) = morningResidue * 1e14;        
        eveningDensities(j,i) = eveningResidue * 1e14;
    end
    p.progress;    
end

p.stop;

% morningDensities = bsxfun(@minus, morningDensities, mean(morningDensities));
% eveningDensities = bsxfun(@minus, eveningDensities, mean(eveningDensities));
morningDensities = morningDensities - mean(morningDensities(:));
eveningDensities = eveningDensities - mean(eveningDensities(:));
morningDensities = morningDensities';
eveningDensities = eveningDensities';

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

end

function [magneticLatitude] = convertToMagneticCoordinates(latitude, longitude, altitude)
% [magneticLatitude] = convertToMagneticCoordinates(latitude, longitude, altitude)

ecefXYZ = geod2ecef(latitude, longitude, altitude)';
ecefToMagTransform = [0.33907 -0.91964 -0.19826; ...
                      0.93826  0.34594  0      ; ...
                      0.06859  0.18602  0.98015];
magXYZ = ecefToMagTransform * ecefXYZ;
r = sqrt(magXYZ(1,:).^2 + magXYZ(2,:).^2 + magXYZ(3,:).^2);
magneticLatitude = pi/2 - acos(magXYZ(3,:) ./ r);
magneticLatitude = magneticLatitude' * 180 / pi;

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


