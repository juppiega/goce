function goceDataAnalyze( varargin )
% goceDataAnalyze( AEFilename, apFilename, densityFilename, F107Filename, absBFilename, computeLatitudes, threshold )
%   Detailed explanation goes here

initialize()
%KAYTA STRUCTEJA!!!!
[ AEFilename, apFilename, densityFilename, F107Filename, absBFilename, computeLatitudes, threshold] ...
    = processInputArguments(varargin, nargin);

[ae, aeNoBg, ap, absB, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg,   ...
morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1min, ...
timestamps1minFixed, timestamps3h, timestamps4min, timestamps4minFixed, density4min, density3h, ...
morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = variables(AEFilename, apFilename, densityFilename, F107Filename, absBFilename, threshold);

for i = 1:cellArrayLength
    plotTimeseries(timestamps1min{i}, timestamps1minFixed{i}, timestamps4min{i}, timestamps4minFixed{i},...
        timestamps3h{i}, ae{i}, ap{i}, absB{i},...
        averagedDensityNoBg{i}, density3h{i}, density4min{i})

    plotAndCalculateCorrelation(timestamps1min{i}, timestamps1minFixed{i}, ae{i}, averagedDensityNoBg{i}, 'AE'); 
    plotAndCalculateCorrelation(timestamps3h{i}, timestamps3h{i}, ap{i}, density3h{i}, 'ap'); 
    plotAndCalculateCorrelation(timestamps4min{i}, timestamps4minFixed{i}, absB{i}, density4min{i}, 'IMF |B|'); 

    plotAndAnalyzeDensityByLatitude(ae{i}, ap{i}, absB{i}, timestamps1min{i}, timestamps1minFixed{i}, timestamps4min{i}, timestamps3h{i}, ...
        morningDensityNoBg{i}, morningTimestamps10s{i}, morningMagneticLatitude{i}, computeLatitudes, 'Morning (~06:00 LST)');
    plotAndAnalyzeDensityByLatitude(ae{i}, ap{i}, absB{i}, timestamps1min{i}, timestamps1minFixed{i}, timestamps4min{i}, timestamps3h{i}, ...
        eveningDensityNoBg{i}, eveningTimestamps10s{i}, eveningMagneticLatitude{i}, computeLatitudes, 'Evening (~18:00 LST)');

    plotAndAnalyzeChangesByOrbit(morningDensityNoBg{i}, morningMagneticLatitude{i}, averagedDensityNoBg{i},...
        timestamps1minFixed{i}, morningTimestamps10s{i}, 'Morning (~06:00 LST)')
    plotAndAnalyzeChangesByOrbit(eveningDensityNoBg{i}, eveningMagneticLatitude{i}, averagedDensityNoBg{i},...
        timestamps1minFixed{i}, eveningTimestamps10s{i}, 'Evening (~18:00 LST)')

    simpleColorMap(morningDensityNoBg{i}, morningMsisDensity{i}, morningMagneticLatitude{i}, morningTimestamps10s{i}, 'Morning');
    simpleColorMap(eveningDensityNoBg{i}, eveningMsisDensity{i}, eveningMagneticLatitude{i}, eveningTimestamps10s{i}, 'Evening');
end

end

function initialize()
% initialize()
clear;
format compact
if matlabpool('size') <= 0
    matlabpool open
end

end

function [ AEFilename, apFilename, densityFilename, F107Filename, absBFilename, computeLatitudes, threshold] ...
            = processInputArguments(inputArgs, numOfInputArgs)
% [ AEFilename, densityFilename, computeLatitudes, threshold] = processInputArguments(varargin, nargin)

if numOfInputArgs == 7
    AEFilename = inputArgs{1};
    apFilename = inputArgs{2};
    densityFilename = inputArgs{3};
    F107Filename = inputArgs{4};
    absBFilename = inputArgs{5};
    computeLatitudes = round(inputArgs{6});
    threshold = inputArgs{7};
    
else
    fprintf(2, '%s %d %s\n', 'Wrong number of input arguments: ', numOfInputArgs, '. See >>help goceDataAnalyze')
    error('goceDataAnalyze: Argin error')
end

end

function [ae, aeNoBg, ap, absB, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1min, ...
 timestamps1minFixed, timestamps3h, timestamps4min, timestamps4minFixed, density4min, density3h, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = variables(AEFilename, apFilename, densityFilename, F10Filename, absBFilename, threshold)
% [ae, averagedDensity, timestamps] = readAEDensityAndTimestamps( AEFilename, densityFilename )

[ae, timestampsAeDatenum] = readAeFile(AEFilename);
intervalsOfInterest = findInterestingIntervals(ae, timestampsAeDatenum, threshold);

apMonth = readApFile(apFilename, timestampsAeDatenum);

[F107values, F107datenum] = readF107File(F10Filename);

[absB, timestampsAbsBDatenum] = readAbsBFile(absBFilename);

[density, averagedDensity, averagedDensityNoBg, density4min, density3h, longitude, latitude, altitude, solarTime, magneticLatitude, ...
    timestamps10sFixed, timestamps1min, timestamps1minFixed, timestamps4min, timestamps4minFixed, timestampsDensityDatenum, ae, aeNoBg] = ...
    readDensityFile(densityFilename, timestampsAeDatenum, timestampsAbsBDatenum, ae);

[F107A81Days, F107Yesterday] = giveF107ValsForMSIS(F107values, F107datenum, timestampsDensityDatenum);

[apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, ApDaily, ap, timestamps3h] =...
    giveApValsForMSIS(apMonth, timestamps10sFixed, timestamps1minFixed);

[densityNoBg, msisDensity] = relateMsisToDensity(density, altitude, timestampsDensityDatenum, timestamps10sFixed,...
    timestamps1min, latitude, longitude, F107A81Days, F107Yesterday, ae, ApDaily, apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h);

[morningTimestamps10s, morningMagneticLatitude, morningDensityNoBg, morningMsisDensity, eveningTimestamps10s,  ...
    eveningMagneticLatitude, eveningDensityNoBg, eveningMsisDensity] = ...
    splitBySolarTime(timestamps10sFixed, magneticLatitude, densityNoBg, msisDensity, solarTime);

[ae, aeNoBg, ap, absB, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestamps4min, ...
 timestamps1minFixed, timestamps4minFixed, timestamps3h, density3h, density4min, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = sliceToInterestingIntervals(ae, aeNoBg, ap, absB, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1minFixed, timestamps4minFixed, timestamps3h, ...
 density3h, density4min, morningMagneticLatitude, eveningMagneticLatitude, timestamps10sFixed, timestamps1min, timestamps4min, intervalsOfInterest);

if fclose('all') ~= 0
    display('File close unsuccesful - check that it is not reserved by another editor.')
end

end

function [ae, timestampsAeDatenum] = readAeFile(AEFilename)
%
aeFile = fopen(AEFilename);
if aeFile == -1
    error('aeFile open unsuccesful')
end

aeData = textscan(aeFile, '%s %s %f %f %f %f %f', 'MultipleDelimsAsOne',1, 'HeaderLines',15);

ae = aeData{4};
timestampAeString = strcat(aeData{1}, aeData{2});
timestampsAeDatenum = datenum(timestampAeString, 'yyyy-mm-ddHH:MM:SS.FFF');

end

function ap = readApFile(apFilename, timestampsAeDatenum)
%
apFile = fopen(apFilename);
if apFile == -1
    error('apFile open unsuccesful')
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
apRows = find(timestampsAp >= timestampsAeDatenum(1) - 3 & timestampsAp <= timestampsAeDatenum(end));
ap = reshape(apValues(apRows, :)', numel(apValues(apRows, :)), 1);

end

function [F107values, F107datenum] = readF107File(F10Filename)
%
F10File = fopen(F10Filename);
if F10File == -1
    error('F107File open unsuccesful')
end

F10data = textscan(F10File, '%f %f %f %f', 'MultipleDelimsAsOne',1);
F107values = F10data{4};
F10year = F10data{1};
F10dayOfYear = F10data{2};
F107datenum = datenum(F10year, 1, 1) + F10dayOfYear - 1;

end

function [absB, timestampsAbsBDatenum] = readAbsBFile(absBFilename)
%

absBFile = fopen(absBFilename);
if absBFile == -1
    error('absBFile open unsuccesful')
end
textscan(absBFile, '%s %s %f', 'commentStyle','#', 'MultipleDelimsAsOne',1);
fgets(absBFile);
absBData = textscan(absBFile, '%s %s %f', 'commentStyle','#', 'MultipleDelimsAsOne',1);
timestampsAbsBString = strcat(absBData{1}, absBData{2});
timestampsAbsBDatenum = datenum(timestampsAbsBString, 'dd-mm-yyyyHH:MM:SS.FFF');
absB = absBData{3};
nanIndices = find(absB < 0);
timestampsAbsBDatenum(nanIndices) = [];
absB(nanIndices) = [];

end

function [density, averagedDensity, averagedDensityNoBg, density4min, density3h, longitude, latitude, altitude, solarTime, magneticLatitude,...
    timestamps10sFixed, timestamps1min, timestamps1minFixed, timestamps4min, timestamps4minFixed, timestampsDensityDatenum, ae, aeNoBg] = ...
    readDensityFile(densityFilename, timestampsAeDatenum, timestampsAbsBDatenum, ae)
%

densityFile = fopen(densityFilename);

if densityFile == -1
    error('densityFile open unsuccesful')
end

densityData = textscan(densityFile, '%s %s %s %f %f %f %f %f %f %f %f %f', 'commentStyle','#');

density = densityData{9} * power(10, 11);
longitude = densityData{5};
latitude = densityData{6};
altitude = densityData{4};
solarTime = densityData{7};
magneticLatitude = convertToMagneticCoordinates(latitude, longitude, altitude);

averagedDensity = smooth(density, 7);

secondsInDay = 60 * 60 * 24;
timestampDensityString = strcat(densityData{1}, densityData{2});
timestampsDensityDatenum = datenum(timestampDensityString, 'yyyy-mm-ddHH:MM:SS.FFF');

%ae = ae(ismember(timestampsAeDatenum, timestampsDensityDatenum));
aeNoBg = removePeriodicBackground(ae, 125, 1, 0);

averagedDensity = averagedDensity(ismember(timestampsDensityDatenum, timestampsAeDatenum));
averagedDensityNoBg = removePeriodicBackground(averagedDensity, 125, 1, 1);
aeNoBg = normalize(aeNoBg, ae);
averagedDensityNoBg = normalize(averagedDensityNoBg, averagedDensity);

timestamps10sFixed = timestampsDensityDatenum * secondsInDay;
timestamps1minFixed = timestamps10sFixed(ismember(timestampsDensityDatenum, timestampsAeDatenum));
monthFirstDatenum = datenum(datestr(timestampsDensityDatenum(1), 'yyyy-mm-dd'));
timestamps10sFixed = timestamps10sFixed - monthFirstDatenum * secondsInDay;
timestamps1minFixed = timestamps1minFixed - monthFirstDatenum * secondsInDay;
timestamps10sFixed = round(timestamps10sFixed);
timestamps1minFixed = round(timestamps1minFixed);
timestamps1min = round((timestampsAeDatenum - monthFirstDatenum) * secondsInDay);
timestamps4min = round((timestampsAbsBDatenum - monthFirstDatenum) * secondsInDay);
figure;
plot(timestamps1min / 86400 + 1, ae);
title('AE For whole month')
ylabel('AE')
xlabel('Calendar days from the beginning of the month');
grid on
threeHinSec = 3 * 60 * 60;
timestamps3h = threeHinSec : threeHinSec : 31 * secondsInDay;
density3h = smooth(averagedDensity, 179);
density3h = density3h(find(ismember(timestamps1minFixed, timestamps3h)) - 90);

density4min = smooth(averagedDensityNoBg, 5);
density4min = density4min(ismember(timestamps1minFixed, timestamps4min));
timestamps4minFixed = timestamps4min(ismember(timestamps4min, timestamps1minFixed)); 

end

function [F107A81Days, F107Yesterday] = giveF107ValsForMSIS(F107values, F107datenum, timestampsDensityDatenum)
%

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

function [apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, ApDaily, apMonthFixed, timestamps3h] = ...
    giveApValsForMSIS(apMonth, timestamps10sFixed, timestamps1min)
% 

threeHinSec = 3 * 60 * 60;
secondsInDay = 24 * 60 * 60;
timestamps3h = -3 * 8 * threeHinSec + threeHinSec : threeHinSec : 31 * secondsInDay;
timestamps10s = min(timestamps10sFixed) : 10 : max(timestamps10sFixed);
nearestThreeHStamp = threeHinSec * floor(timestamps10s / threeHinSec) + threeHinSec;
ap10s = reshape(repmat(apMonth', 180 * 6, 1),[],1);
apTime10s = reshape(repmat(timestamps3h, 180 * 6, 1),[],1);
ap10s(1) = [];
apTime10s(1) = [];
apNow = ap10s(ismember(apTime10s, nearestThreeHStamp));
ap3h = ap10s(ismember(apTime10s, nearestThreeHStamp - threeHinSec));
ap6h = ap10s(ismember(apTime10s, nearestThreeHStamp - 2 * threeHinSec));
ap9h = ap10s(ismember(apTime10s, nearestThreeHStamp - 3 * threeHinSec));
ap10sSmoothed = reshape(repmat(smooth(apMonth, 9)', 180 * 6, 1),[],1);
ap10sSmoothed(1) = [];
apAver12To33h = ap10sSmoothed(ismember(apTime10s, nearestThreeHStamp - 8 * threeHinSec));
apAver36To57h = ap10sSmoothed(ismember(apTime10s, nearestThreeHStamp - 16 * threeHinSec));
ApDaily = ap10sSmoothed(ismember(apTime10s, nearestThreeHStamp - 4 * threeHinSec));
finalIndices = ismember(timestamps10s, timestamps10sFixed);
apNow = apNow(finalIndices);
ap3h = ap3h(finalIndices);
ap6h = ap6h(finalIndices);
ap9h = ap9h(finalIndices);
apAver12To33h = apAver12To33h(finalIndices);
apAver36To57h = apAver36To57h(finalIndices);
ApDaily = ApDaily(finalIndices);

timestamps3h(timestamps3h <= 0) = -1;
apMonthFixed = apMonth(ismember(timestamps3h, timestamps1min));
timestamps3h = timestamps3h(ismember(timestamps3h, timestamps1min))';

end

function [morningTimestamps10s, morningMagneticLatitude, morningDensityNoBg, morningMsisDensity, eveningTimestamps10s,  ...
    eveningMagneticLatitude, eveningDensityNoBg, eveningMsisDensity] = ...
    splitBySolarTime(timestamps10s, magneticLatitude, densityNoBg, msisDensity, solarTime)
%

morningIndices = find(solarTime <= 12);
eveningIndices = find(solarTime > 12);

morningTimestamps10s = timestamps10s(morningIndices);
morningMagneticLatitude = magneticLatitude(morningIndices);
morningDensityNoBg = densityNoBg(morningIndices);
morningMsisDensity = msisDensity(morningIndices);

eveningTimestamps10s = timestamps10s(eveningIndices);
eveningMagneticLatitude = magneticLatitude(eveningIndices);
eveningDensityNoBg = densityNoBg(eveningIndices);
eveningMsisDensity = msisDensity(eveningIndices);

end

function intervalsOfInterest = findInterestingIntervals(ae, timestampsAeDatenum, threshold)
%

intervalsOfInterest = zeros(1, 2);
medianCrossings = findCrossings(ae, 'smooth', 'mean');
calmDays = 2;
if length(medianCrossings) > 1
    for i = 1:length(medianCrossings) - 1
        if ~isempty(find(ae(medianCrossings(i):medianCrossings(i + 1)) >= threshold, 1))
            [intervalBegin, intervalEnd] = giveInterval(medianCrossings(i), medianCrossings(i + 1), timestampsAeDatenum, calmDays);
            intervalsOfInterest = vertcat(intervalsOfInterest, [intervalBegin intervalEnd]);
        end
    end
end

intervalsOfInterest(1,:) = [];
if isempty(intervalsOfInterest)
    fprintf(2, '%s %d %s\n', 'There were no AE values above threshold: ', threshold, '. Either the threshold is too high,')
    fprintf(2, '%s %d %s\n', 'or length(medianCrossings): ', length(medianCrossings), ' == 1, which should not happen.')
    error('No AE values above threshold found!')
end

if length(intervalsOfInterest(:,1)) > 1
    newIntervals = intervalsOfInterest;
    indicesInDay = 24 * 60;
    numOfGoodIntervals = 0;
    while numOfGoodIntervals ~= length(intervalsOfInterest(:,1)) - 1
        for i = 1:length(intervalsOfInterest(:,1)) - 1
            if intervalsOfInterest(i, 2) - intervalsOfInterest(i + 1,1) > calmDays * indicesInDay
                newIntervals(i,:) = [intervalsOfInterest(i,1) intervalsOfInterest(i + 1,2)];
                newIntervals(i + 1,:) = [];
                break;
            else
                numOfGoodIntervals = numOfGoodIntervals + 1;
            end
        end
        intervalsOfInterest = newIntervals;
    end
end

end

function [intervalBegin, intervalEnd] = giveInterval(peakBeginIndex, peakEndIndex, timestampsAeDatenum, calmDays)
%

firstAeDay = timestampsAeDatenum(1);
lastAeDay = timestampsAeDatenum(end);
peakBegin = timestampsAeDatenum(peakBeginIndex);
peakEnd = timestampsAeDatenum(peakEndIndex);

if floor(peakBegin - calmDays) < firstAeDay
    intervalBegin = firstAeDay;
else
    intervalBegin = floor(peakBegin - calmDays);
end

if ceil(peakEnd + calmDays) > lastAeDay
    intervalEnd = floor(lastAeDay);
else
    intervalEnd = ceil(peakEnd + calmDays);
end

intervalBegin = find(timestampsAeDatenum == intervalBegin);
intervalEnd = find(timestampsAeDatenum == intervalEnd);

end

function [ae, aeNoBg, ap, absB, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1minOut, timestamps4min, ...
 timestamps1minFixed, timestamps4minFixed, timestamps3h, density3h, density4min, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = sliceToInterestingIntervals(ae, aeNoBg, ap, absB, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1minFixed, timestamps4minFixed, timestamps3h, ...
 density3h, density4min, morningMagneticLatitude, eveningMagneticLatitude, timestamps10sFixed, timestamps1min, timestamps4min, intervalsOfInterest)
%

aeTemp = ae; aeNoBgTemp = aeNoBg; apTemp = ap; averagedDensityTemp = averagedDensity;...
averagedDensityNoBgTemp = averagedDensityNoBg; morningDensityNoBgTemp = morningDensityNoBg; ...
eveningDensityNoBgTemp = eveningDensityNoBg; morningMsisDensityTemp = morningMsisDensity; ...
eveningMsisDensityTemp = eveningMsisDensity; morningTimestamps10sTemp = morningTimestamps10s; ...
eveningTimestamps10sTemp = eveningTimestamps10s; timestamps1minFixedTemp = timestamps1minFixed; ...
timestamps3hTemp = timestamps3h; density3hTemp = density3h; morningMagneticLatitudeTemp = morningMagneticLatitude;...
eveningMagneticLatitudeTemp = eveningMagneticLatitude; timestamps4minTemp = timestamps4min; density4minTemp = density4min; ...
timestamps4minFixedTemp = timestamps4minFixed; absBTemp = absB;

ae = {}; aeNoBg = {}; ap = {}; averagedDensity = {}; averagedDensityNoBg = {}; morningDensityNoBg = {};
eveningDensityNoBg = {}; morningMsisDensity = {}; eveningMsisDensity = {}; morningTimestamps10s = {};
eveningTimestamps10s = {}; timestamps1minFixed = {}; timestamps3h = {}; density3h = {}; morningMagneticLatitude = {};
eveningMagneticLatitude = {}; timestamps4min = {}; density4min = {}; timestamps4minFixed = {}; absB = {};

cellArrayLength = length(intervalsOfInterest(:,1));

for i = 1:cellArrayLength
    beginIndex = intervalsOfInterest(i,1);
    endIndex = intervalsOfInterest(i,2);
    timestamps1minOut{i} = timestamps1min(beginIndex:endIndex);
    ae{i} = aeTemp(beginIndex:endIndex);
    aeNoBg{i} = aeNoBgTemp(beginIndex:endIndex);
    
    threeHinSec = 3 * 60 * 60;
    [~, beginIndex3h] = min(abs(timestamps3hTemp - threeHinSec * round(timestamps1min(beginIndex) / threeHinSec)));
    [~, endIndex3h] = min(abs(timestamps3hTemp - threeHinSec * round(timestamps1min(endIndex) / threeHinSec)));
    ap{i} = apTemp(beginIndex3h:endIndex3h);
    timestamps3h{i} = timestamps3hTemp(beginIndex3h:endIndex3h);
    density3h{i} = density3hTemp(beginIndex3h:endIndex3h);
    
    fourMinInSec = 4 * 60;
    [~, beginIndex4min] = min(abs(timestamps4minTemp - fourMinInSec * round(timestamps1min(beginIndex) / fourMinInSec)));
    [~, endIndex4min] = min(abs(timestamps4minTemp - fourMinInSec * round(timestamps1min(endIndex) / fourMinInSec)));
    absB{i} = absBTemp(beginIndex4min:endIndex4min);
    timestamps4min{i} = timestamps4minTemp(beginIndex4min:endIndex4min);
    
    [~, beginIndex4minFixed] = min(abs(timestamps4minFixedTemp - fourMinInSec * round(timestamps1min(beginIndex) / fourMinInSec)));
    [~, endIndex4minFixed] = min(abs(timestamps4minFixedTemp - fourMinInSec * round(timestamps1min(endIndex) / fourMinInSec)));
    density4min{i} = density4minTemp(beginIndex4minFixed:endIndex4minFixed);
    timestamps4minFixed{i} = timestamps4minFixedTemp(beginIndex4minFixed:endIndex4minFixed);
    
    [~, beginIndex1min] = min(abs(timestamps1minFixedTemp - timestamps1min(beginIndex)));
    [~, endIndex1min] = min(abs(timestamps1minFixedTemp - timestamps1min(endIndex)));
    averagedDensity{i} = averagedDensityTemp(beginIndex1min:endIndex1min);
    averagedDensityNoBg{i} = averagedDensityNoBgTemp(beginIndex1min:endIndex1min);
    timestamps1minFixed{i} = timestamps1minFixedTemp(beginIndex1min:endIndex1min);
    
    [~, beginIndexMorning10s] = min(abs(morningTimestamps10sTemp - timestamps1min(beginIndex)));
    [~, endIndexMorning10s] = min(abs(morningTimestamps10sTemp - timestamps1min(endIndex)));   
    [~, beginIndexEvening10s] = min(abs(eveningTimestamps10sTemp - timestamps1min(beginIndex)));
    [~, endIndexEvening10s] = min(abs(eveningTimestamps10sTemp - timestamps1min(endIndex)));
    morningMsisDensity{i} = morningMsisDensityTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningMsisDensity{i} = eveningMsisDensityTemp(beginIndexEvening10s:endIndexEvening10s);
    morningDensityNoBg{i} = morningDensityNoBgTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningDensityNoBg{i} = eveningDensityNoBgTemp(beginIndexEvening10s:endIndexEvening10s);
    morningTimestamps10s{i} = morningTimestamps10sTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningTimestamps10s{i} = eveningTimestamps10sTemp(beginIndexEvening10s:endIndexEvening10s);
    morningMagneticLatitude{i} = morningMagneticLatitudeTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningMagneticLatitude{i} = eveningMagneticLatitudeTemp(beginIndexEvening10s:endIndexEvening10s);
end

end

function valueNoBg = normalize(valueNoBg, value)
% [aeNoBg, averagedDensityNoBg] = normalize(aeNoBg, averagedDensityNoBg)

valueNoBg = mean(value) / mean(valueNoBg) * valueNoBg;
end

function [values] = removePeriodicBackground( values, removeShorterPeriods, samplingFreq, plotOrNot )
% [averagedDensity] = removePeriodicBackground( values, removeShorterPeriods, samplingFreq, plotOrNot )

numOfSamples = length(values);
valueFFT = fft(values);
Freq = (0:numOfSamples - 1) * samplingFreq / numOfSamples;
period = 1 ./ Freq;

if plotOrNot > 0
    figure;
    plotIndices = find(period < 500);
    plot(period(plotIndices), abs(valueFFT(plotIndices)))
    title('Density FFT')
    xlabel('T / min')
end

indicesToRemove = period < removeShorterPeriods;
valueFFT(indicesToRemove) = 0;

values = abs(ifft(valueFFT));

end

function [valueNearPeak, peakBegin, peakEnd] = limitToNearPeak(value, smoothOption, meanOrMedian)
% [aeNearPeak] = limitToNearPeak(ae)

% figure;
% subplot(2,1,1)
% plot(1:length(value), value);
[medianCrossings, subtractedVal] = findCrossings(value, smoothOption, meanOrMedian);

integralVals = zeros(length(medianCrossings), 1);
if length(medianCrossings) == 1
   if trapz(subtractedVal(medianCrossings:end)) >= trapz(subtractedVal(1:medianCrossings))
      peakBegin = medianCrossings + 1;
      peakEnd = length(value);
      integralVals = trapz(subtractedVal(medianCrossings:end));
   else
       peakBegin = 1;
       peakEnd = medianCrossings;
       integralVals = trapz(subtractedVal(1:medianCrossings));
   end
else
    for i = 1:length(medianCrossings) - 1
        integralVals(i) = trapz(subtractedVal(medianCrossings(i):medianCrossings(i+1)));
    end
    maxIndex = find(integralVals == max(integralVals));
    peakBegin = medianCrossings(maxIndex) + 1;
    peakEnd = medianCrossings(maxIndex + 1);
end

value([(1:peakBegin - 1) (peakEnd + 1:length(value))]) = 0;
valueNearPeak = value;
% subplot(2,1,2)
% plot(1:length(value), value);
end

function [crossings, subtractedVal] = findCrossings(value, smoothOption, meanOrMedian)
%

%fftSmoothedValue = removePeriodicBackground(value, 432, 1, 0);
%subtractedVal = fftSmoothedValue - mean(fftSmoothedValue);
if strcmpi(meanOrMedian, 'mean')
    if strcmpi(smoothOption, 'smooth')
        fftSmoothedValue = removePeriodicBackground(value, 432, 1, 0);
        subtractedVal = fftSmoothedValue - mean(fftSmoothedValue);
    else
        subtractedVal = value - mean(value);
    end
else
    if strcmpi(smoothOption, 'smooth')
        fftSmoothedValue = removePeriodicBackground(value, 432, 1, 0);
        subtractedVal = fftSmoothedValue - median(fftSmoothedValue);
    else
        subtractedVal = value - median(value);
    end
end

crossings = find(subtractedVal(1:end-1) .* subtractedVal(2:end) < 0);

end

function plotTimeseries(timestamps1min, timestamps1minFixed, timestamps4min, timestamps4minFixed, timestamps3h, ae, ap, absB, ...
    averagedDensityNoBg, density3h, density4min)
% plotTimeseries(timestamps, ae, averagedDensity)
global timeseriesFigHandle
timeseriesFigHandle = figure;
subplot(2,2,1)
secondsInDay = 60 * 60 * 24;
timestampsInDays1min = timestamps1min / secondsInDay + 1;
timestampsInDays4min = timestamps4min / secondsInDay + 1;
timestampsInDays4minFixed = timestamps4minFixed / secondsInDay + 1;
[hAx,~,~] = plotyy(timestampsInDays4min, absB, timestampsInDays4minFixed, density4min);
title('Timeseries of IMF |B| and density')
xlabel('t / days')
ylabel(hAx(1), 'IMF |B| / nT')
ylabel(hAx(2), 'Density')
set(hAx, 'XLim', [min(timestampsInDays1min) max(timestampsInDays1min)]);
grid on

subplot(2,2,3)
timestampsInDays3h = timestamps3h / secondsInDay + 1;
[hAx,~,~] = plotyy(timestampsInDays3h, ap, timestampsInDays3h, density3h);
title('Timeseries of ap and density')
xlabel('t / days')
ylabel(hAx(1), 'ap')
ylabel(hAx(2), 'previous 3h average density')
set(hAx, 'XLim', [min(timestampsInDays3h) max(timestampsInDays3h)]);
grid on

timestampsInDays1minFixed = timestamps1minFixed / secondsInDay + 1;
subplot(2,2,4)
[hAx,~,~] = plotyy(timestampsInDays1min, ae, timestampsInDays1minFixed, averagedDensityNoBg);
title('Timeseries of AE and density')
xlabel('t / days')
ylabel(hAx(1), 'AE')
ylabel(hAx(2), 'Density')
set(hAx, 'XLim', [min(timestampsInDays1min) max(timestampsInDays1min)]);
grid on

end

function [densityIndexTimelag] = plotAndGiveMaxCrossCorrelation(density, geomIndex, indexName)
% plotCrossCorrelation(averagedDensityNoBg, ae)

if strcmpi(indexName, 'ae'); maxLag = 24 * 60; else maxLag = 8; end
[correlations, lags] = xcorr(density, geomIndex, maxLag, 'coeff');
correlations = correlations(lags >= 0);
lags = lags(lags >= 0);
if strcmpi(indexName, 'ae'); lagsInHours = lags/60; else lagsInHours = lags * 3; end

figure;
plot(lagsInHours, correlations);
title([indexName, ' and Density xcorr'])
xlabel('lags / h')
densityIndexTimelag = lagsInHours(correlations == max(correlations));
fprintf('\n%s %f %s\n\n', [indexName, ' - Density timelag:'], densityIndexTimelag, 'h')
ylimits = get(gca, 'ylim');
line([densityIndexTimelag densityIndexTimelag], [ylimits(1) max(correlations)], 'LineStyle', '--');
textYLocation = mean([ylimits(1) max(correlations)]);
text(densityIndexTimelag, textYLocation, ['\leftarrow Lag = ', num2str(densityIndexTimelag), ' h'], 'FontSize', 14)
if strcmpi(indexName, 'ae'); densityIndexTimelag = densityIndexTimelag * 60; else densityIndexTimelag = densityIndexTimelag / 3;end

end

function [ae, averagedDensity] = moveDataseries(ae, averagedDensity, timelag)
% [ae, averagedDensity] = moveDataseries(ae, averagedDensity, timelag)
timelag = round(timelag);
if timelag > 0
   averagedDensity = averagedDensity((timelag + 1) : end);
   ae = ae(1 : (length(ae) - timelag));
elseif timelag < 0
   averagedDensity = averagedDensity(1 : (length(averagedDensity) - timelag));
   ae = ae((timelag + 1) : end);
end

end

function plotAndCalculateCorrelation(timestamps, timestampsFixed, geomIndex, density, indexName)
% [r, r2] = plotAndCalculateCorrelation(ae, averagedDensity, timelag)

global timeseriesFigHandle
if exist('timeseriesFigHandle', 'var')
    timeseriesHandle = timeseriesFigHandle;
end
%[densityNearPeak] = limitToNearPeak(density);
%[geomIndexNearPeak] = limitToNearPeak(geomIndex);
%figure;
%plotyy(1:length(geomIndexNearPeak), geomIndexNearPeak, 1:length(densityNearPeak), densityNearPeak)
[timelag, bestIntegral] = compareDensityToGeomIndexIntegral(density, geomIndex, timestamps, timestampsFixed, indexName);
% if strcmpi(indexName, 'ae')
%     [timelag, bestIntegral] = compareDensityToGeomIndexIntegral(density, geomIndex, indexName);
% elseif ~isempty(strfind(upper(indexName), '|B|'));
%     [timelag, bestIntegral] = compareDensityToGeomIndexIntegral(density, geomIndex, indexName);
% else
%     [timelag, bestIntegral] = compareDensityToGeomIndexIntegral(density, geomIndex, indexName);
% end
%[timelag] = plotAndGiveMaxCrossCorrelation(densityNearPeak, geomIndexNearPeak, indexName);
%[geomIndex, density] = moveDataseries(geomIndex, density, timelag);
geomIndexFixed = geomIndex(ismember(timestamps, timestampsFixed));
densityFixed = density(ismember(timestampsFixed, timestamps));
timestampsFixed = timestampsFixed(ismember(timestampsFixed, timestamps));
densityIndices = ismember(timestampsFixed, timestamps(timelag + 1:end));
    
if strcmpi(indexName, 'ae')
    figure(timeseriesHandle);
    subplot(2,2,2)
    secondsInDay = 24 * 60 * 60;
    timestampsInDays = timestamps(timelag + 1:end) / secondsInDay + 1;
    timestampsInDaysFixed = timestampsFixed(densityIndices) / secondsInDay + 1;
    [hAx,~,~] = plotyy(timestampsInDays, bestIntegral, timestampsInDaysFixed, density(densityIndices));
    xlabel('t / days')
    ylabel(hAx(1), 'AE Integral')
    ylabel(hAx(2), 'Density')
    title('Raw AE integral vs FFT smoothed density');
    set(hAx, 'XLim', [min(timestampsInDaysFixed) max(timestampsInDaysFixed)]);    
    grid on;
end

bestIntegral = bestIntegral(ismember(timestamps(timelag + 1:end), timestampsFixed));
plotCorrelation(bestIntegral, density(densityIndices), [indexName, ' Integral'], 'Density');

geomIndexNoBg = removePeriodicBackground(geomIndexFixed, 125, 1, 0);
geomIndexNoBg = normalize(geomIndexNoBg, geomIndexFixed);
if strcmpi(indexName, 'ae') || ~isempty(strfind(upper(indexName), '|B|'))
    plotCorrelation(geomIndexNoBg, densityFixed, indexName, 'Density at 270 km');
else
    plotCorrelation(geomIndexFixed, densityFixed, indexName, 'Density at 270 km');
end

end

function plotCorrelation(xvals, yvals, xvalName, yvalName)
%

r = corr(xvals, yvals);
r2 = r * r;
fprintf('%s %f\n', 'Pearson correlation: ', r)
fprintf('%s %d %s\n', 'Thus, ', round(r2 * 100), ['% of variation in ', yvalName, ' can be explained by changes in ', xvalName])

% xScaled = (xvals - mean(xvals)) ./ std(xvals);
% yScaled = (yvals - mean(yvals)) ./ std(yvals);
figure;
plot(xvals, yvals, '.')
p = polyfit(xvals, yvals, 1);
m = p(1);
b = p(2);
ylimits = get(gca, 'ylim');
xlimits = get(gca, 'xlim');
regLineYAxisCross = m * xlimits(1) + b;
regLineXAxisCross = (ylimits(2) - b) / m;
line([xlimits(1) regLineXAxisCross], [regLineYAxisCross ylimits(2)]);
title([yvalName, ' vs ', xvalName ,' correlation'])
xlabel(xvalName)
ylabel(yvalName)
textYLocation = ylimits(2) - 0.05 * (ylimits(2) - ylimits(1));
textXLocation = xlimits(2) - 0.05 * (xlimits(2) - xlimits(1));
str1 = [' r =  ', num2str(r, '%07.4f')];
str2 = ['r^2 = ', num2str(r2, '%07.4f')];
textString = [str1 ; str2];
text(textXLocation, textYLocation, textString, 'FontSize', 12, ...
    'VerticalAlignment','top', 'HorizontalAlignment','right', 'EdgeColor', [0 0 0]);

end

function [bestLag, geomIndexBestInt] = compareDensityToGeomIndexIntegral(density, ...
    geomIndex, timestamps, timestampsFixed, indexName)
%

maxDays = 3;
if strcmpi(indexName, 'ae'); maxLag = maxDays * 24 * 60; 
elseif ~isempty(strfind(upper(indexName), '|B|')); maxLag = (maxDays * 24 * 60) / 4;
else maxLag = maxDays * 8; end

cumulativeGeomIndex = cumtrapz(geomIndex);
correlations = zeros(1, maxLag + 1);
lags = zeros(1, maxLag + 1);
correlationType = 'Spearman';
parfor lag = 1:maxLag
    geomIndexInt = cumulativeGeomIndex(lag + 1 : end) - cumulativeGeomIndex(1 : end - lag);
    timestampsInt = timestamps(lag + 1 : end);
    densityIndices = ismember(timestampsFixed, timestampsInt);
    integralIndices = ismember(timestampsInt, timestampsFixed);
    correlations(lag + 1) = corr(density(densityIndices), geomIndexInt(integralIndices), 'type', correlationType);
    lags(lag + 1) = lag;
end
correlations(1) = corr(density(ismember(timestampsFixed, timestamps)), geomIndex(ismember(timestamps, timestampsFixed)),...
    'type', correlationType);
lags(1) = 0;
if strcmpi(indexName, 'ae'); lagsInHours = lags / 60; 
elseif ~isempty(strfind(upper(indexName), '|B|')); lagsInHours = lags * 4 / 60;
else lagsInHours = lags * 3; end

bestLag = lags(correlations == max(correlations));
geomIndexBestInt = cumulativeGeomIndex(bestLag + 1 : end) - cumulativeGeomIndex(1 : end - bestLag);

figure;
plot(lagsInHours, correlations);
title([indexName, ' integral optimal window length'])
ylabel([correlationType, ' correlation']);
xlabel('lags / h')
integralWindowSize = lagsInHours(correlations == max(correlations));
fprintf('\n%s %f %s\n\n', [indexName, ' integral - Density timelag:'], integralWindowSize, 'h')
ylimits = get(gca, 'ylim');
line([integralWindowSize integralWindowSize], [ylimits(1) max(correlations)], 'LineStyle', '--');
textYLocation = mean([ylimits(1) max(correlations)]);
text(integralWindowSize, textYLocation, ['\leftarrow Lag = ', num2str(integralWindowSize), ' h'], 'FontSize', 14)
% if strcmpi(indexName, 'ae'); integralWindowSize = integralWindowSize * 60;
% elseif ~isempty(strfind(upper(indexName), '|B|')); integralWindowSize = integralWindowSize * 60 / 4;
% else integralWindowSize = integralWindowSize / 3;end

end

function plotDensityLatitudeTimeSurf(correctedDensity, magneticLatitude, timestamps)
% plotDensityLatitudeTimeSurf(averagedDensity, averagedLatitude, timestamps)

figure;
secondsInDay = 60 * 60 * 24;
timestampsInDays = timestamps / secondsInDay; % EI nyt ole kylla kalenteripaivia!
[tInterp, latitudeInterp] = meshgrid(0:120/secondsInDay:timestampsInDays(length(timestampsInDays)), ...
    min(magneticLatitude):1:max(magneticLatitude));
densityInterp = griddata(timestampsInDays, magneticLatitude, correctedDensity, tInterp, latitudeInterp, 'cubic');
surf(tInterp, latitudeInterp, densityInterp, 'EdgeColor', 'None')
caxis([min(correctedDensity) max(correctedDensity)])
colormap jet(500)
colorbar
view(2);
xlabel('t / days')
ylabel('Lat (°)')
zlabel('Density')

% hold on;
% maxZ = ones(size(timestamps)) * max(correctedDensity);
% plot3(timestampsInDays, averagedLatitude, maxZ, 'r--');
% hold off;
end

function plotAndAnalyzeDensityByLatitude(ae, ap, absB, timestamps1min, timestamps1minFixed, timestamps4min, timestamps3h, correctedDensity, ...
    timestamps10s, magneticLatitude, computeLatitudes, timeOfDay)
% plotCorrectedDensityLatitudes(ae, timestamps1min, correctedDensity, timestamps10s, latitude, timestampsDatenum, computeLatitudes);
    
[latBeginIndex, latEndIndex] = limitLatitudeToIntegerMultipleOfOrbitalPeriod(magneticLatitude);
limitedLatitude = magneticLatitude(latBeginIndex:latEndIndex);
limitedTimestamps = timestamps10s(latBeginIndex:latEndIndex);
[~, exactOrbitIndices] = splitIntoOrbits(limitedLatitude);
limitedLatitude = limitedLatitude(exactOrbitIndices);
limitedTimestamps = limitedTimestamps(exactOrbitIndices);
[orbits, ~] = splitIntoOrbits(limitedLatitude);

[minAllowedLatitude, maxAllowedLatitude] = findInterpolationLimits(limitedLatitude);
orbitsToDelete = find(abs(limitedLatitude(orbits(:,2)) - limitedLatitude(orbits(:,1))) < ...
    (maxAllowedLatitude - minAllowedLatitude));
newIndices = 1:length(limitedLatitude);
for i = 1:length(orbitsToDelete)
    newIndices = setdiff(newIndices, (orbits(orbitsToDelete(i), 1) : orbits(orbitsToDelete(i), 2)));
end
limitedTimestamps = limitedTimestamps(newIndices);
limitedLatitude = limitedLatitude(newIndices);

F = scatteredInterpolant(timestamps10s, magneticLatitude, correctedDensity);

tooSmallComputeLatIndices = find(computeLatitudes < minAllowedLatitude);
if ~isempty(tooSmallComputeLatIndices)
    fprintf(2,'%s %.1f %s %d%s\n', 'Warning: latitude(s)', computeLatitudes(tooSmallComputeLatIndices), ...
        ['is/are smaller than the minimum allowed by dataset (', timeOfDay, '): '], minAllowedLatitude,...
           '. Using min latitude instead.')
       computeLatitudes(tooSmallComputeLatIndices) = minAllowedLatitude;
end

tooBigComputeLatIndices = find(computeLatitudes > maxAllowedLatitude);
if ~isempty(tooBigComputeLatIndices)
    fprintf(2,'%s %.1f %s %d%s\n', 'Warning: latitude(s)', computeLatitudes(tooBigComputeLatIndices), ...
        ['is/are greater than the maximum allowed by dataset (', timeOfDay, '): '], maxAllowedLatitude,...
           '. Using max latitude instead.')
       computeLatitudes(tooBigComputeLatIndices) = maxAllowedLatitude;
end
computeLatitudes = unique(computeLatitudes, 'stable');

oneDegreeStep = minAllowedLatitude:1:maxAllowedLatitude;

plotDensityLength = length(latitudeCrossingTimes(limitedLatitude, limitedTimestamps, minAllowedLatitude));
plotDensityByLatitude = zeros(plotDensityLength, 1);
plotDensityCrossingTimes = zeros(plotDensityLength, 1);
parfor i = 1:length(oneDegreeStep)
    %oneDegreeStep(i) = checkOutOfBounds(minAllowedLatitude, maxAllowedLatitude, oneDegreeStep(i));
    crossingTimes(:,i) = latitudeCrossingTimes(limitedLatitude, limitedTimestamps, oneDegreeStep(i));

    latitudeInterp = ones(length(crossingTimes(:,i)), 1) * oneDegreeStep(i);    
    densityByLatitude(:,i) = F(crossingTimes(:,i), latitudeInterp);

    if ~isempty(find(computeLatitudes == oneDegreeStep(i), 1))
        plotDensityByLatitude = horzcat(plotDensityByLatitude, densityByLatitude(:,i));
        plotDensityCrossingTimes = horzcat(plotDensityCrossingTimes, crossingTimes(:,i));
    end
end
plotDensityByLatitude = plotDensityByLatitude(:, 2:end);
plotDensityCrossingTimes = plotDensityCrossingTimes(:, 2:end);
limitedTimestamps = limitedTimestamps(ismember(limitedTimestamps, timestamps1minFixed));

plotDensity(limitedTimestamps, computeLatitudes, plotDensityByLatitude, plotDensityCrossingTimes, ...
    densityByLatitude, minAllowedLatitude, maxAllowedLatitude, crossingTimes, timestamps1min, timestamps4min, timestamps3h, ae, ap, absB, timeOfDay, 0);
plotDensity(limitedTimestamps, computeLatitudes, plotDensityByLatitude, plotDensityCrossingTimes, ...
    densityByLatitude, minAllowedLatitude, maxAllowedLatitude,crossingTimes, timestamps1min, timestamps4min, timestamps3h, ae, ap, absB, timeOfDay, 1);


end

function [orbits, exactOrbitIndices] = splitIntoOrbits(latitude)
%

if satelliteIsGoingSouth(latitude)
    orbitEndIndices = [find(latitude(1:end-1) < latitude(2:end)); length(latitude)];
else
    orbitEndIndices = [find(latitude(1:end-1) > latitude(2:end)); length(latitude)];
end
orbitBeginIndices = [1; (orbitEndIndices(1:end-1) + 1)];

orbits = [orbitBeginIndices orbitEndIndices];

orbitsToDelete = orbits(:,2) == orbits(:,1);
orbits = orbits(~orbitsToDelete, :);

exactOrbitIndices = [];
for i = 1:length(orbits(:,1))
    exactOrbitIndices = [exactOrbitIndices (orbits(i,1) : orbits(i,2))];
end

end

function [minAllowedLatitude, maxAllowedLatitude] = findInterpolationLimits(magneticLatitude)
%
maxAllowedLatitude = 90;
numOfSameNumOfCrossings = 0;
previousNumOfCrossings = -1;
while 1
   subtractedLatitude = magneticLatitude - maxAllowedLatitude;
   numOfCrossings = length(find(subtractedLatitude(1:end-1) .* subtractedLatitude(2:end) < 0));
   if numOfCrossings == previousNumOfCrossings
      numOfSameNumOfCrossings = numOfSameNumOfCrossings + 1; 
   else
      numOfSameNumOfCrossings = 0;
   end
   
   if numOfSameNumOfCrossings > 2
       break;
   end
   
   previousNumOfCrossings = numOfCrossings;
   maxAllowedLatitude = maxAllowedLatitude - 1;
end

minAllowedLatitude = -90;
numOfSameNumOfCrossings = 0;
previousNumOfCrossings = -1;
while 1
   subtractedLatitude = magneticLatitude - minAllowedLatitude;
   numOfCrossings = length(find(subtractedLatitude(1:end-1) .* subtractedLatitude(2:end) < 0));
   if numOfCrossings == previousNumOfCrossings
      numOfSameNumOfCrossings = numOfSameNumOfCrossings + 1; 
   else
      numOfSameNumOfCrossings = 0;
   end
   
   if numOfSameNumOfCrossings > 2
       break;
   end
   
   previousNumOfCrossings = numOfCrossings;
   minAllowedLatitude = minAllowedLatitude + 1;
end

end

function [crossingTimes] = latitudeCrossingTimes(limitedLatitude, limitedTimestamps, computeLatitude)
% latitudeCrossingTimes(limitedLatitude, limitedTimestamps, computeLatitude)

subtractedLatitude = limitedLatitude - computeLatitude;
crossingIndices = find(subtractedLatitude(1:end-1) .* subtractedLatitude(2:end) < 0);
crossingIndices = crossingIndices(abs(limitedLatitude(crossingIndices) - limitedLatitude(crossingIndices + 1)) < 20);
[crossingIndices] = checkConsistency(crossingIndices);
tBeforeCrossing = limitedTimestamps(crossingIndices);
tAfterCrossing = limitedTimestamps(crossingIndices + 1);
subLatBeforeCrossing = abs(subtractedLatitude(crossingIndices));
subLatAfterCrossing = abs(subtractedLatitude(crossingIndices + 1));
crossingTimes = (subLatAfterCrossing .* tBeforeCrossing + subLatBeforeCrossing .* tAfterCrossing) ...
                                    ./ (subLatAfterCrossing + subLatBeforeCrossing);
end

function plotDensity(limitedTimestamps, computeLatitudes, plotDensityByLatitude, plotDensityCrossingTimes, densityByLatitude, ...
    minAllowedLatitude, maxAllowedLatitude, crossingTimes, timestamps1min, timestamps4min, timestamps3h, ae, ap, absB, timeOfDay, fftSmoothing)
% plotDensity(limitedTimestamps, computeLatitudes, plotDensityByLatitude, plotDensityCrossingTimes, densityByLatitude, timestamps1min, timestampsDatenum, ae, fftSmoothing)

plotTitle = ['270 km density for different ', lower(timeOfDay), ' latitudes'];
if fftSmoothing ~= 0
    secondsInDay = 60 * 60 * 24;
    parfor i = 1:length(computeLatitudes)
        samplingFreq = secondsInDay / mean(diff(plotDensityCrossingTimes(:,i)));
        plotDensityByLatitude(:,i) = removePeriodicBackground(plotDensityByLatitude(:,i), 0.3, samplingFreq, 0);
        %densityByLatitude(:,i) = removePeriodicBackground(densityByLatitude(:,i), 0.25, samplingFreq, 0);
    end
    plotTitle = ['FFT filtered 270 km density for different ', lower(timeOfDay), ' latitudes'];
end

parfor i = 1:length(computeLatitudes)
    splinedDensity(:,i) = spline(plotDensityCrossingTimes(:,i), plotDensityByLatitude(:,i), limitedTimestamps);
end

if fftSmoothing ~= 0
    timestampsInDays = limitedTimestamps / (60 * 60 * 24) + 1;
    %timestampsInDays = crossingTimes / (60 * 60 * 24) + 1;
    figure;
    subplot(2,1,1)
    t = repmat(timestampsInDays,1,length(computeLatitudes));
    plot(t,splinedDensity)
    %plot(timestampsInDays, densityByLatitude);
    legends = strsplit(num2str(computeLatitudes));
    legend(legends)
    xlabel('t / days');
    ylabel('Density relative to MSIS model');
    title(plotTitle)
    grid on
        
    secondsInDay = 60 * 60 * 24;
    timestampsInDays1min = timestamps1min / secondsInDay + 1;
    timestampsInDays4min = timestamps4min / secondsInDay + 1;
    timestampsInDays3h = timestamps3h / secondsInDay + 1;

    subplot(2,1,2)
    [hAx,~,~] = plotyy(timestampsInDays1min, ae, timestampsInDays4min, 10 * absB);
    hold(hAx(1), 'on'); hold(hAx(2), 'on');
    plot(hAx(2), timestampsInDays3h, ap, 'r');
    hold off;
    xlabel('t / days')
    ylabel(hAx(1), 'AE')
    ylabel(hAx(2), 'ap / IMF |B| [10^{-10} T]')
end
%axis([min(timestampsInDays) max(timestampsInDays) min(ae) max(ae)]);

if fftSmoothing == 0
    ae = ae(ismember(timestamps1min, limitedTimestamps));
    writeAndPlotPeakAnalysis(limitedTimestamps, ae, densityByLatitude, crossingTimes, minAllowedLatitude, maxAllowedLatitude, timeOfDay);
end
end

function [crossingIndices] = checkConsistency(crossingIndices)
% [crossingIndices] = checkConsistency(crossingIndices)

diffOfIndices = diff(crossingIndices);
inconsistentIndices = find(diffOfIndices > 1.5 * mean(diffOfIndices));
incIndicesApproxVals = round((crossingIndices(inconsistentIndices) + crossingIndices(inconsistentIndices + 1)) ./ 2);
for i = 1:length(inconsistentIndices)
   crossingIndices = [crossingIndices(1:inconsistentIndices(i)); incIndicesApproxVals(i); ...
                        crossingIndices(inconsistentIndices(i) + 1:end)]; 
   inconsistentIndices = inconsistentIndices + 1;
end
end

function writeAndPlotPeakAnalysis(limitedTimestamps, ae, densityByLatitude, crossingTimes, minAllowedLatitude, maxAllowedLatitude, timeOfDay)
% writeAndPlotPeakAnalysis(timestamps, ae, splinedDensity, computeLatitutes)

global aeDensityLagByLatHandle
%global densityPeakMaxTimesByLatHandle

ae = removePeriodicBackground(ae, 432, 1, 0); % 432 min =^ 0.3 days
% ae = smooth(ae, 431);
[ae, ~, ~] = limitToNearPeak(ae, 'noSmooth', 'mean');
for i = 1:length(densityByLatitude(1,:))
   samplingFreq = 86400 / mean(diff(crossingTimes(:,i)));
   densityByLatitude(:,i) = removePeriodicBackground(densityByLatitude(:,i), 0.3, samplingFreq, 0);
end

oneDegreeStep = minAllowedLatitude:maxAllowedLatitude;
parfor i = 1:length(oneDegreeStep)
    splinedDensity(:,i) = spline(crossingTimes(:,i), densityByLatitude(:,i), limitedTimestamps);
end

step = 8;
northernInterval = -step/2:step:maxAllowedLatitude;
northernInterval(ismember(northernInterval, maxAllowedLatitude)) = [];
allLatitudesAeLag = zeros(maxAllowedLatitude - minAllowedLatitude + 1, 1);
catTimelag = 0;
for i = northernInterval
    densityInterval = i - minAllowedLatitude + 1 : i - minAllowedLatitude + 1 + step;
    if ~isempty(find(densityInterval > maxAllowedLatitude - minAllowedLatitude + 1, 1))
        densityInterval = i - minAllowedLatitude + 1 : maxAllowedLatitude - minAllowedLatitude + 1;
        catTimelag = 1;
    end
    
    maxlag = 24 * 60;
    lagsAllIntervalLats = zeros(length(densityInterval), 1);
    maxIndicesAllIntervalLats = zeros(length(densityInterval), 1);
    for k = densityInterval
       [correlations, lags] = xcorr(splinedDensity(:,k), ae, maxlag, 'coeff');
       correlations = correlations(lags >= 0);
       lags = lags(lags >= 0); 
       allLatitudesAeLag(k) = lags(correlations == max(correlations));
       lagsAllIntervalLats(k - min(densityInterval) + 1) = lags(correlations == max(correlations));
       [~, maxIndicesAllIntervalLats(k - min(densityInterval) + 1)] = max(splinedDensity(:,k));
    end
    
    tolerance = 60 * 6;
    maxIndicesAllIntervalLats = fixPossibleOutliers(maxIndicesAllIntervalLats, splinedDensity, tolerance);
    
    if catTimelag > 0
        timelag = [timelag mean(lagsAllIntervalLats)];
        errInLag = [errInLag (std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats)))];  
        maxTimes = limitedTimestamps(maxIndicesAllIntervalLats);
        maxTimesAverage = [maxTimesAverage mean(maxTimes)];
        errInMax = [errInMax (std(maxTimes) / sqrt(length(maxTimes)))];
    else
        timelag(length(i:-step:minAllowedLatitude) + 1) = mean(lagsAllIntervalLats);
        errInLag(length(i:-step:minAllowedLatitude) + 1) = (std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats)));
        maxTimes = limitedTimestamps(maxIndicesAllIntervalLats);
        maxTimesAverage(length(i:-step:minAllowedLatitude) + 1) = mean(maxTimes);
        errInMax(length(i:-step:minAllowedLatitude) + 1) = std(maxTimes) / sqrt(length(maxTimes));
    end
end

southernInterval = -step/2:-step:minAllowedLatitude;
southernInterval(ismember(southernInterval, minAllowedLatitude)) = [];
for i = southernInterval
    densityInterval = i - minAllowedLatitude + 1 - step : i - minAllowedLatitude + 1;
    if ~isempty(find(densityInterval < 1, 1))
        densityInterval = 1 : i - minAllowedLatitude + 1;
    end
    
    maxlag = 24 * 60;
    lagsAllIntervalLats = zeros(length(densityInterval), 1);
    maxIndicesAllIntervalLats = zeros(length(densityInterval), 1);
    for k = densityInterval
       [correlations, lags] = xcorr(splinedDensity(:,k), ae, maxlag, 'coeff');
       correlations = correlations(lags >= 0);
       lags = lags(lags >= 0); 
       allLatitudesAeLag(k) = lags(correlations == max(correlations));
       lagsAllIntervalLats(k - min(densityInterval) + 1) = lags(correlations == max(correlations));
       [~, maxIndicesAllIntervalLats(k - min(densityInterval) + 1)] = max(splinedDensity(:,k));      
    end
    
    maxIndicesAllIntervalLats = fixPossibleOutliers(maxIndicesAllIntervalLats, splinedDensity, tolerance);
    
    timelag(length(i:-step:minAllowedLatitude)) = mean(lagsAllIntervalLats);
    errInLag(length(i:-step:minAllowedLatitude)) = std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats));
    maxTimes = limitedTimestamps(maxIndicesAllIntervalLats);
    maxTimesAverage(length(i:-step:minAllowedLatitude)) = mean(maxTimes);
    errInMax(length(i:-step:minAllowedLatitude)) = std(maxTimes) / sqrt(length(maxTimes));
end

densityAeTimelag = mean(allLatitudesAeLag) / 60;
densityAeError = std(allLatitudesAeLag) / 60;
densityAeError = ceil(densityAeError * 10) / 10;
fprintf('%s %s %s %.1f %s %.1f % s\n', 'AE - Density ', timeOfDay, ' timelag: ', densityAeTimelag, '±',densityAeError, 'h')

intervalLatitudes = (southernInterval(end) : step : northernInterval(end)) + step/2;
intervalLatitudes(end) = [];
intervalLatitudes = [ round(mean([southernInterval(end) minAllowedLatitude]))...
                      intervalLatitudes                                      ...
                      round(mean([northernInterval(end) maxAllowedLatitude])) ]; 

if ~isempty(strfind(lower(timeOfDay), 'morning')); aeDensityLagByLatHandle = figure; subplotNum = 1; else subplotNum = 2; end
figure(aeDensityLagByLatHandle);
subplot(1,2,subplotNum)
timestampsInHours = timelag / 60;
errInHours = errInLag / 60;
herrorbar(timestampsInHours, intervalLatitudes, errInHours, '-s')
xlabel('t / h');
ylabel('Geomagnetic latitude');
title(timeOfDay)
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Density-AE lag for different latitudes', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

% if ~isempty(strfind(lower(timeOfDay), 'morning')); densityPeakMaxTimesByLatHandle = figure; subplotNum = 1; else subplotNum = 2; end
% figure(densityPeakMaxTimesByLatHandle);
% subplot(1,2,subplotNum)
% timestampsInHours = (maxTimesAverage - limitedTimestamps(ae == max(ae))) / 3600;
% errInMaxInHours = errInMax / 3600;
% herrorbar(timestampsInHours, intervalLatitudes, errInMaxInHours, '-s');
% xlabel('t / h');
% ylabel('Geomagnetic latitude');
% title(timeOfDay)
% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', 'Density maximum times for different latitudes', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center')

end

function maxIndicesAllIntervalLats = fixPossibleOutliers(maxIndicesAllIntervalLats, splinedDensity, tolerance)
%

outliers = find(maxIndicesAllIntervalLats > median(maxIndicesAllIntervalLats) + tolerance |...
                maxIndicesAllIntervalLats < median(maxIndicesAllIntervalLats) - tolerance);
newMean = mean(maxIndicesAllIntervalLats(~ismember(1:length(maxIndicesAllIntervalLats), outliers)));
if outliers > 0
   if round(newMean) - tolerance > 0; beginIndex = round(newMean) - tolerance; else beginIndex = 1; end
   if round(newMean) + tolerance < length(splinedDensity(:,1)); endIndex = round(newMean) + tolerance; 
       else endIndex = length(splinedDensity(:,1)); end
   lookInInterval = beginIndex : endIndex;
   [~, maxIndicesAllIntervalLats(outliers)] = max(splinedDensity(lookInInterval, outliers));
   maxIndicesAllIntervalLats(outliers) = maxIndicesAllIntervalLats(outliers) + min(lookInInterval);
end

end

function [latBeginIndex, latEndIndex] = limitLatitudeToIntegerMultipleOfOrbitalPeriod(latitude)
% [latBeginIndex, latEndIndex] = limitLatitudeToIntegerMultipleOfOrbitalPeriod(latitude)

orbitalPeriod = 45 * 6; % of 10 second ticks
goingSouth = satelliteIsGoingSouth(latitude);
if goingSouth == 1
    latBeginIndex = find(latitude(1:orbitalPeriod) == max(latitude(1:orbitalPeriod)));
    latEndIndex = find(latitude(end-orbitalPeriod : end) == min(latitude(end-orbitalPeriod : end)));
else
    latBeginIndex = find(latitude(1:orbitalPeriod) == min(latitude(1:orbitalPeriod)));
    latEndIndex = find(latitude(end-orbitalPeriod : end) == max(latitude(end-orbitalPeriod : end)));
end
latEndIndex = latEndIndex + (length(latitude) - orbitalPeriod - 1);
end

function simpleColorMap(correctedDensity, msisDensity, latitude, timestamps10s, timeOfDay)
% simpleDataPlot(correctedDensity, latitude, timestamps10s)

figure;
timestampsInDays = timestamps10s / (60 * 60 * 24) + 1;
%timestampsInDays = timestampsInDays - timestampsInDays(1);
scatter(timestampsInDays, latitude, 60, correctedDensity, '.')
xlabel('t / d')
ylabel('latitude')
title(['GOCE ', timeOfDay, ' density at 270 km'])
colorbar
grid on

figure;
scatter(timestampsInDays, latitude, 60, msisDensity, '.')
xlabel('t / d')
ylabel('latitude')
title(['NRLMSISE00 ', timeOfDay, ' density at 270 km'])
colorbar
grid on

end

function [correctedDensity, msisDensity270km] = relateMsisToDensity(density, altitude, timestampsDatenum, timestamps10s, ...
    timestamps1min, latitude,longtitude, F107A, F107, ae, ApDaily, apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h)
% [correctedDensity] = subtractMsisFromDensity(density, altitude,
% timestampsDatenum, latitude, longtitude, F107A, F107, ap)

correctedDensity = zeros(size(density));
msisDensity270km = zeros(size(density));
msisDensityVariableAlt = zeros(size(density));

modelingIndices = 1:length(density);

altitudeKm = altitude / 1000;

parfor i = modelingIndices
    modelVariableAlt = run_nrlmsise00(altitudeKm(i), timestampsDatenum(i), latitude(i), longtitude(i), F107A(i), F107(i),...
        ApDaily(i), apNow(i), ap3h(i), ap6h(i), ap9h(i), apAver12To33h(i), apAver36To57h(i));
    modelConst270km = run_nrlmsise00(270, timestampsDatenum(i), latitude(i), longtitude(i), F107A(i), F107(i),...
        ApDaily(i), apNow(i), ap3h(i), ap6h(i), ap9h(i), apAver12To33h(i), apAver36To57h(i));
    correctedDensity(i) = density(i) * modelConst270km(:,7) / modelVariableAlt(:,7);
    msisDensity270km(i) = 1000 * modelConst270km(:,7);
    msisDensityVariableAlt(i) = 1000 * modelVariableAlt(:,7);
    %correctedDensity(i) = density(i) / (1000 * modelVariableAlt(:,7));
end

msisDensityVariableAlt = msisDensityVariableAlt * power(10, 11);
ratio = density ./ msisDensityVariableAlt;
ratioTrend = removePeriodicBackground(ratio, 125, 6, 0);
figure;
timestampsInDays10s = timestamps10s / 86400 + 1;
timestampsInDays1min = timestamps1min / 86400 + 1;
[hAx,hLine1,hLine2] = plotyy(timestampsInDays10s, ratio, timestampsInDays1min, ae);
set(hLine1, 'LineStyle', 'none', 'Marker', '.')
% For some reason Matlab can't handle axis limits properly, so they need to
% be set manually.
ratioYmin = 0.1 * floor(min(ratio) / 0.1);
ratioYmax = 0.1 * ceil(max(ratio) / 0.1);
set(hAx(1), 'YLim', [ratioYmin ratioYmax], 'YTick', ratioYmin:0.1:ratioYmax);
set(hAx, 'XLim', [min(timestampsInDays10s) max(timestampsInDays10s)], 'XMinorTick', 'on'); 
aeYmax = 500 * ceil(max(ae) / 500);
set(hAx(2),'YLim', [0 aeYmax], 'YTick', 0:250:aeYmax); 
hold on;
plot(timestampsInDays10s, ratioTrend, 'r-');
hold off;
%plotyy(timestampsInDays, correctedDensity, timestampsInDays, msisDensity270km);
xlabel('t / day of month')
ylabel(hAx(1), 'GOCE / NRLMSISE00 density')
ylabel(hAx(2), 'AE')
title('Ratio of GOCE NON-normalized densities to NRLMSISE00 prediction')

meanRatio = mean(ratio);
errRatio = std(ratio);
meanReciprocal = mean(1 ./ ratio);
errReciprocal = std(1 ./ ratio);
textString1 = ['Mean GOCE / NRLMSISE00 density: ', num2str(meanRatio, '%.2f'), ' ± ', num2str(errRatio, '%.2f')];
textString2 = ['Mean NRLMSISE00 / GOCE density: ', num2str(meanReciprocal, '%.2f'), ' ± ', num2str(errReciprocal, '%.2f')];
     
ylimits = get(gca, 'ylim');
xlimits = get(gca, 'xlim');
textYLocation = ylimits(2) - 0.05 * (ylimits(2) - ylimits(1));
textXLocation = xlimits(2) - 0.05 * (xlimits(2) - xlimits(1));
text(textXLocation, textYLocation, textString1, 'FontSize', 9, ...
    'VerticalAlignment','top', 'HorizontalAlignment','right');

textYLocation = ylimits(2) - 0.1 * (ylimits(2) - ylimits(1));
textXLocation = xlimits(2) - 0.05 * (xlimits(2) - xlimits(1));
text(textXLocation, textYLocation, textString2, 'FontSize', 9, ...
    'VerticalAlignment','top', 'HorizontalAlignment','right');

plotCorrelation(msisDensityVariableAlt, density, 'NRLMSISE00 density', 'GOCE measured density')

msisDensity270km = msisDensity270km * power(10, 11);

end

function plotAndAnalyzeChangesByOrbit(densityNoBg, magneticLatitude, averagedDensityNoBg, timestamps1minFixed, ...
    timestamps10s, timeOfDay)
% plotAndAnalyzeChangesByOrbit(densityNoBg, magneticLatitude, averagedDensityNoBg)
global densityByOrbitFigHandle
global residueFigHandle
global densityByOrbitAxesHandle
global residueAxisHandle

averagedDensityNoBg = averagedDensityNoBg(ismember(timestamps1minFixed, timestamps10s));
timestamps1minFixed = timestamps1minFixed(ismember(timestamps1minFixed, timestamps10s));
[~, peakBeginIndex, peakEndIndex] = limitToNearPeak(averagedDensityNoBg, 'noSmooth', 'mean');
peakBeginIndex = find(timestamps10s == timestamps1minFixed(peakBeginIndex));
peakEndIndex = find(timestamps10s == timestamps1minFixed(peakEndIndex));
orbits = splitIntoOrbits(magneticLatitude);

if ~isempty(strfind(lower(timeOfDay), 'morning')); densityByOrbitFigHandle = figure; subplotNum = 1; else subplotNum = 2; end
figure(densityByOrbitFigHandle);
densityByOrbitAxesHandle(subplotNum) = subplot(2,1,subplotNum);
%axis([-90 90 min(densityNoBg) max(densityNoBg)])
hold all;
linehandles = [];
relativeResidues = nan(size(magneticLatitude));
TADplot = nan(size(magneticLatitude));
calmOrbits = 5;
maxNumOfColorOrbits = 10;
[beginOrbit, endOrbit] = findBeginAndEndOrbits(orbits, peakBeginIndex, peakEndIndex, calmOrbits);
if endOrbit - beginOrbit - 2 * calmOrbits > maxNumOfColorOrbits
    endColorOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits - 1;
    loopOrbits = [(beginOrbit : endColorOrbit) (endOrbit - calmOrbits + 1 : endOrbit)]; 
else
    loopOrbits = beginOrbit:endOrbit;
    endColorOrbit = endOrbit;
end

surfPlotIndices = nan(size(magneticLatitude));
for i = 1:length(loopOrbits)
    indices = orbits(loopOrbits(i),1) : orbits(loopOrbits(i),2);
    smoothedDensity150s = smooth(densityNoBg(indices), 15);
    relativeResidues(indices) = (densityNoBg(indices) - smoothedDensity150s) ./ smoothedDensity150s;
    h = plot(magneticLatitude(indices), smoothedDensity150s, 'LineWidth', 2);
    linehandles = [linehandles h];
    
    if loopOrbits(i) <= endColorOrbit
        surfPlotIndices(indices) = indices; 
        smoothedDensity10400km = smooth(densityNoBg(indices), 133);
        smoothedDensity2600km = smooth(densityNoBg(indices), 33);
        TADplot(indices) = (smoothedDensity2600km - smoothedDensity10400km) ./ smoothedDensity10400km;
    end
    
end

%hold off;
%(calmOrbits + 1 : length(linehandles) - calmOrbits);
stormOrbits = (calmOrbits + 1 : length(linehandles) - calmOrbits);
[~,cmap] = cmapline('colormap',jet,'filled', 'lines', linehandles(stormOrbits)');
colormap(cmap)
colorbar

set(linehandles(1:calmOrbits), 'LineStyle', '--')
set(linehandles(1:calmOrbits), 'Color', 'k')
endLineHandles = length(linehandles) - calmOrbits + 1 : length(linehandles);
set(linehandles(endLineHandles), 'LineStyle', '-')
set(linehandles(endLineHandles), 'Color', 'k')
hold off;
xlabel('IGRF geomagnetic latitude')
ylabel('Density at 270 km')
title(timeOfDay)
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '150-s smoothed density along orbit', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
if ~isempty(strfind(lower(timeOfDay), 'evening'))
   linkaxes(densityByOrbitAxesHandle) 
end

magneticLatitudeResiduePlot = magneticLatitude(~isnan(relativeResidues));
relativeResidues = relativeResidues(~isnan(relativeResidues));

confIntervalStep = 10;
confIntervalX = [];
confIntervalMean = [];
confIntervalLower = [];
confIntervalUpper = [];
for i = -90:confIntervalStep:90 - confIntervalStep
    indicesInInterval = find(magneticLatitudeResiduePlot < i + 5 & magneticLatitudeResiduePlot >= i);
    if ~isempty(indicesInInterval)
        confIntervalX = [confIntervalX mean([i (i + confIntervalStep)])];
        confIntervalMean = [confIntervalMean mean(relativeResidues(indicesInInterval))];
        [~, ~, ci, ~] = ttest(relativeResidues(indicesInInterval));        
        confIntervalLower = [confIntervalLower ci(1)];
        confIntervalUpper = [confIntervalUpper ci(2)];
    end
end

if ~isempty(strfind(lower(timeOfDay), 'morning')); 
    residueFigHandle = figure('Color','w');
    residueAxisHandle = [];
    confIntervalColor = 'b'; 
    axisNum = 1;
else
    confIntervalColor = 'r';
    axisNum = 2;
end
figure(residueFigHandle);
residueAxisHandle = [residueAxisHandle ciplot(confIntervalLower, confIntervalUpper, confIntervalX, confIntervalColor)];
%if ~isempty(strfind(lower(timeOfDay), 'evening'))
xlim([min(confIntervalX) max(confIntervalX)]);
hold all;
plot(confIntervalX, confIntervalMean, 'k--');
hold off;
%end
xlabel('IGRF geomagnetic latitude')
ylabel('Relative residues 95 % confidence')
title('100-100 km variations')
if axisNum == 1
    hold on
else
    hold off
    legend(residueAxisHandle, '<-Morning', 'Evening->');
end

TADplot = TADplot(~isnan(TADplot));
indices = surfPlotIndices(~isnan(surfPlotIndices));
figure;
timestampsInDays = timestamps10s(indices) / 86400 + 1;
scatter(timestampsInDays, magneticLatitude(indices), 60, TADplot, '.')
xlabel('t / days')
ylabel('IGRF Magnetic Latitude')
title(['1300-5200km changes (TADs) [(2600 km smooth - 10400 km smooth) / 10400 km smooth] ', timeOfDay])
colorbar

end

function [orbitBegin, orbitEnd] = findBeginAndEndOrbits(orbits, peakBeginIndex, peakEndIndex, calmOrbits)
%

orbitBegin = find(orbits(:,1) <= peakBeginIndex, 1, 'last');
if orbitBegin > calmOrbits
    orbitBegin = orbitBegin - calmOrbits;
else
    orbitBegin = 1;
end

orbitEnd = find(orbits(:,2) >= peakEndIndex, 1, 'first');
if orbitEnd < length(orbits(:,2)) - calmOrbits
    orbitEnd = orbitEnd + calmOrbits;
else
    orbitEnd = length(orbits(:,2));
end

end

function goingSouth = satelliteIsGoingSouth(magneticLatitude)
%

orbitEndIndices = find(abs(diff(magneticLatitude)) > 45);
testIndex1 = round(mean([orbitEndIndices(3) orbitEndIndices(2)]));
testIndex2 = testIndex1 + 1;

if magneticLatitude(testIndex1) > magneticLatitude(testIndex2)
   goingSouth = 1; 
else
   goingSouth = 0;
end

end

function [magneticLatitude] = convertToMagneticCoordinates(latitude, longtitude, altitude)
% [magneticLatitude] = convertToMagneticCoordinates(latitude, longtitude, altitude)

ecefXYZ = geod2ecef(latitude, longtitude, altitude)';
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
% 
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

function hh = herrorbar(x, y, l, u, symbol)
%HERRORBAR Horizontal Error bar plot.
%   HERRORBAR(X,Y,L,R) plots the graph of vector X vs. vector Y with
%   horizontal error bars specified by the vectors L and R. L and R contain the
%   left and right error ranges for each point in X. Each error bar
%   is L(i) + R(i) long and is drawn a distance of L(i) to the right and R(i)
%   to the right the points in (X,Y). The vectors X,Y,L and R must all be
%   the same length. If X,Y,L and R are matrices then each column
%   produces a separate line.
%
%   HERRORBAR(X,Y,E) or HERRORBAR(Y,E) plots X with error bars [X-E X+E].
%   HERRORBAR(...,'LineSpec') uses the color and linestyle specified by
%   the string 'LineSpec'. See PLOT for possibilities.
%
%   H = HERRORBAR(...) returns a vector of line handles.
%
%   Example:
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      herrorbar(x,y,e)
%   draws symmetric horizontal error bars of unit standard deviation.
%
%   This code is based on ERRORBAR provided in MATLAB.   
%
%   See also ERRORBAR

% Copyright (c) 2009, Jos van der Geest
% All rights reserved.
% 
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

%   Jos van der Geest
%   email: jos@jasen.nl
%
%   File history:
%   August 2006 (Jos): I have taken back ownership. I like to thank Greg Aloe from
%   The MathWorks who originally introduced this piece of code to the
%   Matlab File Exchange. 
%   September 2003 (Greg Aloe): This code was originally provided by Jos
%   from the newsgroup comp.soft-sys.matlab:
%   http://newsreader.mathworks.com/WebX?50@118.fdnxaJz9btF^1@.eea3ff9
%   After unsuccessfully attempting to contact the orignal author, I
%   decided to take ownership so that others could benefit from finding it
%   on the MATLAB Central File Exchange.

if min(size(x))==1,
    npt = length(x);
    x = x(:);
    y = y(:);
    if nargin > 2,
        if ~isstr(l),
            l = l(:);
        end
        if nargin > 3
            if ~isstr(u)
                u = u(:);
            end
        end
    end
else
    [npt,n] = size(x);
end

if nargin == 3
    if ~isstr(l)
        u = l;
        symbol = '-';
    else
        symbol = l;
        l = y;
        u = y;
        y = x;
        [m,n] = size(y);
        x(:) = (1:npt)'*ones(1,n);;
    end
end

if nargin == 4
    if isstr(u),
        symbol = u;
        u = l;
    else
        symbol = '-';
    end
end

if nargin == 2
    l = y;
    u = y;
    y = x;
    [m,n] = size(y);
    x(:) = (1:npt)'*ones(1,n);;
    symbol = '-';
end

u = abs(u);
l = abs(l);

if isstr(x) | isstr(y) | isstr(u) | isstr(l)
    error('Arguments must be numeric.')
end

if ~isequal(size(x),size(y)) | ~isequal(size(x),size(l)) | ~isequal(size(x),size(u)),
    error('The sizes of X, Y, L and U must be the same.');
end

tee = (max(y(:))-min(y(:)))/100; % make tee .02 x-distance for error bars
% changed from errorbar.m
xl = x - l;
xr = x + u;
ytop = y + tee;
ybot = y - tee;
n = size(y,2);
% end change

% Plot graph and bars
hold_state = ishold;
cax = newplot;
next = lower(get(cax,'NextPlot'));

% build up nan-separated vector for bars
% changed from errorbar.m
xb = zeros(npt*9,n);
xb(1:9:end,:) = xl;
xb(2:9:end,:) = xl;
xb(3:9:end,:) = NaN;
xb(4:9:end,:) = xl;
xb(5:9:end,:) = xr;
xb(6:9:end,:) = NaN;
xb(7:9:end,:) = xr;
xb(8:9:end,:) = xr;
xb(9:9:end,:) = NaN;

yb = zeros(npt*9,n);
yb(1:9:end,:) = ytop;
yb(2:9:end,:) = ybot;
yb(3:9:end,:) = NaN;
yb(4:9:end,:) = y;
yb(5:9:end,:) = y;
yb(6:9:end,:) = NaN;
yb(7:9:end,:) = ytop;
yb(8:9:end,:) = ybot;
yb(9:9:end,:) = NaN;
% end change


[ls,col,mark,msg] = colstyle(symbol); if ~isempty(msg), error(msg); end
symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid

h = plot(xb,yb,esymbol); hold on
h = [h;plot(x,y,symbol)];

if ~hold_state, hold off; end

if nargout>0, hh = h; end

end


function varargout = cmapline(varargin)
% CMAPLINE - Apply a colormap to lines in plot
%
%   CMAPLINE finds all lines in an axis and specifies their
%   colors according to a colormap. Also accepts custom 
%   colormaps in the form of a n x 3 matrix. 
%    
% OPTIONS and SYNTAX
%
%   cmapline - with no inputs, cmapline finds all lines in the 
%   current axis and applies the colormap 'jet'.
% 
%   cmapline('ax',gca,'colormap','hot') - will find all lines
%   in the specified axis (in this case, the current axis)
%   and applies the colormap 'hot'.
%
%   cmapline('lines',handles) - applies colormap values to line
%   objects with specified handles.   
%
%   cmapline('filled') - will fill markers (if included in the 
%   line) with corresponding colormap colors.
%
%   lineh=cmapline - The optional output variable returns the
%   handles to the line objects.
%
%   [lineh, cmap]=cmapline - Two optional outputs returns both the 
%   the handles to the line objects and the applied colormap. 
%
% EXAMPLE 1 - color lines in two subplots according to different colormaps
%  
%   %generate some data
%   x=(0:0.3:2*pi);
%   m=10;
%   exdata=bsxfun(@plus,repmat(10.*sin(x),[m 1]),[1:m]');
%   
%   figure
%   subplot(121);
%   plot(x,exdata,'o-','linewidth',2)
%   cmapline('colormap','jet');
%   set(gca,'color','k')
%   title('jet colormap')
%
%   subplot(122);
%   plot(x,exdata,'o-','linewidth',2)
%   custommap=flipud(hot);
%   cmapline('colormap',custommap,'filled')
%   set(gca,'color','k')
%   title('reverse hot colormap, filled markers')  
%
% EXAMPLE 2 (uses data from example 1) - add a colorbar to your plot
%
%   figure
%   plot(x,exdata,'linewidth',2)
%   [lh,cmap]=cmapline('colormap','jet');
%   colormap(cmap)
%   colorbar
%
% SEE ALSO  colormap 

% Copyright (c) 2010, Andrew Stevens
% All rights reserved.
% 
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

% Andrew Stevens @ USGS, 8/15/2008
% astevens@usgs.gov

%default values
ax=gca;
cmap=@jet;
fillflag=0;
lh=[];

%parse inputs and do some error-checking
if nargin>0
    [m,n]=size(varargin);
    opts={'ax','lines','colormap','filled'};

    for i=1:n;
        indi=strcmpi(varargin{i},opts);
        ind=find(indi==1);
        if isempty(ind)~=1
            switch ind
                case 1
                    %make sure input is an axes handle, sort of
                    ax=varargin{i+1};
                    if ~ishandle(ax)
                        error(['Specified axes',...
                            ' must be a valid axis handle'])
                    end
                case 2
                    lh=varargin{i+1};
                    if ~all(ishandle(lh))
                        error('Invalid line handle')
                    else
                        lh=num2cell(lh);
                    end
                    
                case 3
                    cmap=varargin{i+1};
                    if isa(cmap,'function_handle')
                        cmap= func2str(cmap);
                    end
                    %check size of numeric colormap input
                    if isa(cmap,'numeric')
                        [m,n]=size(cmap);
                        if n~=3
                            error('Custom colormap must have 3 columns.')
                        end
                    end
                case 4
                    fillflag=1;

            end
        else
        end
    end
end

%find lines in axes
if isempty(lh)
    lh=num2cell(findobj(ax,'type','line'));
end

numlines=numel(lh);
if isempty(lh)
    fprintf('No lines present in specified axes.\n')
end

if isa(cmap,'numeric')
    %if needed, interpolate colormap to number of lines
    if numlines~=m
        int=m/numlines;
        ivec=1:m;
        ovec=1:int:1+(numlines-1)*int;

        cmap=num2cell(cmap,1);
        cmap=cellfun(@(x)(interp1(ivec,x,ovec,...
            'linear','extrap')'),cmap,'uni',0);
        colrs=num2cell(cell2mat(cmap),2);
    else
        colrs=num2cell(cmap,2);
    end
else
    %if standard colormap is supplied
    colrs=num2cell(feval(cmap,numlines),2);
end

%apply colors to lines
cellfun(@(x,y)(set(x,'color',y)),lh,colrs);

if strcmpi(get(lh{1},'marker'),'none')~=1 && ...
        fillflag==1;
    cellfun(@(x,y)(set(x,'markerfacecolor',y)),...
        lh,colrs);
end

%output 
if nargout>0
    varargout{1}=cell2mat(lh);
end
if nargout>1
    varargout{2}=cell2mat(colrs);
end

end


function h = ciplot(lower,upper,x,colour)
     
% ciplot(lower,upper)       
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.

% Raymond Reynolds 24/11/06

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end

if nargin<4
    colour='b';
end

if nargin<3
    x=1:length(lower);
end

% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end

h = fill([x fliplr(x)],[upper fliplr(lower)],colour);
%camlight; lighting gouraud;
alpha(0.5)

end




