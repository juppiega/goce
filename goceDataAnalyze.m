function goceDataAnalyze( varargin )
% goceDataAnalyze( threshold )
%   Detailed explanation goes here

initialize()
results = {};
%KAYTA STRUCTEJA!!!!
[threshold, plotDates] = processInputArguments(varargin, nargin);

[ae, ap, absB, akasofuEpsilon, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = variables(threshold, results);

plotFigures = 0;
for i = 1:cellArrayLength
    if userRequestedThisStormToPlot(plotDates, timestampsDatenum{i})
        plotFigures = 1;
    end
    results = writeTimeIntervalAndMaxAeToResultArray(timestampsDatenum{i}, ae{i}, results);
    progress(i, cellArrayLength);
    
%     if plotFigures ~= 0
%         plotTimeseries(timestamps1min{i}, timestamps1minFixed{i}, timestampsAbsB{i},...
%             timestamps3h{i}, timestamps3hFixed{i}, ae{i}, ap{i}, absB{i},averagedDensityNoBg{i}, density3h{i})
%     end
% 
%     results = plotAndCalculateCorrelation(timestamps1min{i}, timestamps1minFixed{i}, ae{i}, averagedDensityNoBg{i}, 'AE', plotFigures, results); 
%     results = plotAndCalculateCorrelation(timestamps3h{i}, timestamps3hFixed{i}, ap{i}, density3h{i}, 'ap', plotFigures, results); 
%     results = plotAndCalculateCorrelation(timestampsAbsB{i}, timestamps1minFixed{i}, absB{i}, averagedDensityNoBg{i}, 'IMF |B|', plotFigures, results); 
%     results = plotAndCalculateCorrelation(timestampsEpsilon{i}, timestamps1minFixed{i}, akasofuEpsilon{i}, averagedDensityNoBg{i}, 'Akasofu Epsilon', plotFigures, results);
%     
<<<<<<< HEAD
%     results = plotAndAnalyzeDensityByLatitude(ae{i}, timestamps1min{i}, timestamps1minFixed{i}, ...
%         morningDensityNoBg{i}, morningTimestamps10s{i}, morningMagneticLatitude{i}, 'Morning', plotFigures, results);
%     results = plotAndAnalyzeDensityByLatitude(ae{i}, timestamps1min{i}, timestamps1minFixed{i}, ...
%         eveningDensityNoBg{i}, eveningTimestamps10s{i}, eveningMagneticLatitude{i}, 'Evening', plotFigures, results);
%     
%     results = plotAndAnalyzeChangesByOrbit(morningDensityNoBg{i}, morningMagneticLatitude{i}, averagedDensityNoBg{i},...
%         timestamps1minFixed{i}, morningTimestamps10s{i}, 'Morning', plotFigures, results);
%     results = plotAndAnalyzeChangesByOrbit(eveningDensityNoBg{i}, eveningMagneticLatitude{i}, averagedDensityNoBg{i},...
%         timestamps1minFixed{i}, eveningTimestamps10s{i}, 'Evening', plotFigures, results);
=======
     [results, morningLatitudes, morningCrossingTimes, morningDensityGrid] = plotAndAnalyzeDensityByLatitude(ae{i}, timestamps1min{i}, timestamps1minFixed{i}, ...
         morningDensityNoBg{i}, morningTimestamps10s{i}, morningMagneticLatitude{i}, 'Morning', plotFigures, results);
%     results = plotAndAnalyzeDensityByLatitude(ae{i}, timestamps1min{i}, timestamps1minFixed{i}, ...
%         eveningDensityNoBg{i}, eveningTimestamps10s{i}, eveningMagneticLatitude{i}, 'Evening', plotFigures, results);
%     
    results = plotAndAnalyzeChangesByOrbit(morningDensityNoBg{i}, morningMagneticLatitude{i}, averagedDensityNoBg{i},...
        timestamps1minFixed{i}, morningTimestamps10s{i}, 'Morning', plotFigures, results);
    results = plotAndAnalyzeChangesByOrbit(eveningDensityNoBg{i}, eveningMagneticLatitude{i}, averagedDensityNoBg{i},...
        timestamps1minFixed{i}, eveningTimestamps10s{i}, 'Evening', plotFigures, results);
>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
    
    if plotFigures ~= 0
        simpleColorMap(morningDensityNoBg{i}, morningMsisDensity{i}, morningMagneticLatitude{i}, morningTimestamps10s{i}, 'Morning');
        %simpleColorMap(eveningDensityNoBg{i}, eveningMsisDensity{i}, eveningMagneticLatitude{i}, eveningTimestamps10s{i}, 'Evening');
<<<<<<< HEAD
        plotDensityLatitudeTimeSurf(morningDensityNoBg{i}, morningMsisDensity{i}, morningMagneticLatitude{i}, morningTimestamps10s{i},'Morning');
=======
        plotDensityLatitudeTimeSurf(morningDensityNoBg{i}, morningMsisDensity{i}, morningMagneticLatitude{i}, morningTimestamps10s{i}, ...
            morningLatitudes, morningCrossingTimes, morningDensityGrid, 'Morning');
>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
        %plotDensityLatitudeTimeSurf(eveningDensityNoBg{i}, eveningMsisDensity{i}, eveningMagneticLatitude{i}, eveningTimestamps10s{i},'Evening');
    end
    
    plotFigures = 0;
end

makeSummaryOfResults(results);

end
<<<<<<< HEAD

function initialize()
% initialize()
clear;
format compact
if matlabpool('size') <= 0
    matlabpool open
end

end

function [threshold, plotDates] = processInputArguments(inputArgs, numOfInputArgs)
% [ AEFilename, densityFilename, threshold] = processInputArguments(varargin, nargin)

if numOfInputArgs == 1
    threshold = inputArgs{1};
    plotDates = -1;
elseif numOfInputArgs > 1
    threshold = inputArgs{1};
    plotDates = zeros(numOfInputArgs - 1, 1);
    for i = 2:numOfInputArgs
        plotDates(i - 1) = datenum(inputArgs{i});
    end
else
    fprintf(2, '%s %d %s\n', 'Wrong number of input arguments: ', numOfInputArgs, '. See >>help goceDataAnalyze')
    error('goceDataAnalyze: Argin error')
end

end

function progress(i, cellArrayLength)
%

fprintf('%s%d%s%d\n', 'Now analyzing storm: ', i, '/', cellArrayLength)

end

function [ae, ap, absB, akasofuEpsilon, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = variables(threshold, results)

% [ae, averagedDensity, timestamps] = readAEDensityAndTimestamps( AEFilename, densityFilename )

if exist('GoceVariables.mat', 'file') == 2
    load('GoceVariables.mat')
else
    fprintf(2, '%s', 'Warning: no GoceVariables.mat found in Matlab PATH, now attempting to create new one -> readFiles.m');
end

%compareGoceDensityToMsis(measuredDensity, msisDensityVariableAlt, ae, timestampsAeDatenum, timestampsDensityDatenum, results);

intervalsOfInterest = findInterestingIntervals(ae, timestampsAeDatenum, epsilonQualityFlag, timestampsEpsilonDatenum, densityNoBg, timestampsDensityDatenum, threshold);

[morningTimestamps10s, morningMagneticLatitude, morningDensityNoBg, morningMsisDensity, eveningTimestamps10s, morningTimestampsDatenum, eveningTimetampsDatenum,  ...
    eveningMagneticLatitude, eveningDensityNoBg, eveningMsisDensity] = ...
    splitBySolarTime(timestamps10sFixed, timestampsDensityDatenum, magneticLatitude, densityNoBg, msisDensity270km, solarTime);

[ae, ap, absB, akasofuEpsilon, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, morningTimestampsDatenum, eveningTimetampsDatenum, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = sliceToInterestingIntervals(ae, ap, absB, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, morningTimestampsDatenum, eveningTimetampsDatenum, timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed,...
 density3h, morningMagneticLatitude, eveningMagneticLatitude, timestamps10sFixed, timestamps1min, timestampsAbsB, akasofuEpsilon, timestampsDensityDatenum, intervalsOfInterest);

end

function results = writeTimeIntervalAndMaxAeToResultArray(timestampsDatenum, ae, results)
%

[rowNum, ~] = size(results);
if rowNum == 0
    results(1,:) = {'Interval start', 'Interval end', 'Decimal year', 'Max AE'};
    rowNum = rowNum + 1;
end
rowNum = rowNum + 1;
results{rowNum,1} = datestr(timestampsDatenum(1), 'yyyy-mm-dd');
results{rowNum,2} = datestr(timestampsDatenum(end), 'yyyy-mm-dd');
decimalYear = (timestampsDatenum(1) + 2 - datenum(datestr(timestampsDatenum(1), 'yyyy'), 'yyyy')) / 366;
results{rowNum,3} = decimalYear;
results{rowNum,4} = max(ae);

end

function plotOrNot = userRequestedThisStormToPlot(plotDates, timestampsDatenum)
%

=======

function initialize()
% initialize()
clear;
format compact
if matlabpool('size') <= 0
    matlabpool open
end

end

function [threshold, plotDates] = processInputArguments(inputArgs, numOfInputArgs)
% [ AEFilename, densityFilename, threshold] = processInputArguments(varargin, nargin)

if numOfInputArgs == 1
    threshold = inputArgs{1};
    plotDates = -1;
elseif numOfInputArgs > 1
    threshold = inputArgs{1};
    plotDates = zeros(numOfInputArgs - 1, 1);
    for i = 2:numOfInputArgs
        plotDates(i - 1) = datenum(inputArgs{i});
    end
else
    fprintf(2, '%s %d %s\n', 'Wrong number of input arguments: ', numOfInputArgs, '. See >>help goceDataAnalyze')
    error('goceDataAnalyze: Argin error')
end

end

function progress(i, cellArrayLength)
%

fprintf('%s%d%s%d\n', 'Now analyzing storm: ', i, '/', cellArrayLength)

end

function [ae, ap, absB, akasofuEpsilon, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = variables(threshold, results)

% [ae, averagedDensity, timestamps] = readAEDensityAndTimestamps( AEFilename, densityFilename )

if exist('GoceVariables.mat', 'file') == 2
    load('GoceVariables.mat')
else
    fprintf(2, '%s', 'Warning: no GoceVariables.mat found in Matlab PATH, now attempting to create new one -> readFiles.m');
end

%compareGoceDensityToMsis(measuredDensity, msisDensityVariableAlt, ae, timestampsAeDatenum, timestampsDensityDatenum, results);

intervalsOfInterest = findInterestingIntervals(ae, timestampsAeDatenum, epsilonQualityFlag, timestampsEpsilonDatenum, densityNoBg, timestampsDensityDatenum, threshold);

[morningTimestamps10s, morningMagneticLatitude, morningDensityNoBg, morningMsisDensity, eveningTimestamps10s, morningTimestampsDatenum, eveningTimetampsDatenum,  ...
    eveningMagneticLatitude, eveningDensityNoBg, eveningMsisDensity] = ...
    splitBySolarTime(timestamps10sFixed, timestampsDensityDatenum, magneticLatitude, densityNoBg, msisDensity270km, solarTime);

[ae, ap, absB, akasofuEpsilon, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, morningTimestampsDatenum, eveningTimetampsDatenum, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = sliceToInterestingIntervals(ae, ap, absB, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, morningTimestampsDatenum, eveningTimetampsDatenum, timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed,...
 density3h, morningMagneticLatitude, eveningMagneticLatitude, timestamps10sFixed, timestamps1min, timestampsAbsB, akasofuEpsilon, timestampsDensityDatenum, intervalsOfInterest);

end

function results = writeTimeIntervalAndMaxAeToResultArray(timestampsDatenum, ae, results)
%

[rowNum, ~] = size(results);
if rowNum == 0
    results(1,:) = {'Interval start', 'Interval end', 'Decimal year', 'Max AE'};
    rowNum = rowNum + 1;
end
rowNum = rowNum + 1;
results{rowNum,1} = datestr(timestampsDatenum(1), 'yyyy-mm-dd');
results{rowNum,2} = datestr(timestampsDatenum(end), 'yyyy-mm-dd');
decimalYear = (timestampsDatenum(1) + 2 - datenum(datestr(timestampsDatenum(1), 'yyyy'), 'yyyy')) / 366;
results{rowNum,3} = decimalYear;
results{rowNum,4} = max(ae);

end

function plotOrNot = userRequestedThisStormToPlot(plotDates, timestampsDatenum)
%

>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
plotOrNot = plotDates >= min(timestampsDatenum) & plotDates <= max(timestampsDatenum);

end

function makeSummaryOfResults(results)
%

equinoxesAndSolstices = [-0.02 0.22 0.47 0.73 0.98];
decimalYears = cell2mat(results(2:end, 3));
springStorms = [];
summerStorms = [];
autumnStorms = [];
winterStorms = [];
for i = 1:length(decimalYears)
    [~, season] = min(abs(equinoxesAndSolstices - decimalYears(i)));
    
    switch season
        case 2
            springStorms = [springStorms; i + 1];
        case 3
            summerStorms = [summerStorms; i + 1];
        case 4
            autumnStorms = [autumnStorms; i + 1];
        otherwise
            winterStorms = [winterStorms; i + 1];
    end
end

fprintf('\n%s\n', 'Seasonal distribution:')
fprintf('%s %d\n', 'Spring: ', length(springStorms));
fprintf('%s %d\n', 'Summer: ', length(summerStorms));
fprintf('%s %d\n', 'Autumn: ', length(autumnStorms));
fprintf('%s %d\n\n', 'Winter: ', length(winterStorms));

% nanIndices = cellfun(@isnan, results(2:end, 3:end), 'uniformoutput', 0);
% [row,col] = find(cell2mat(nanIndices));
% means = cellfun(@nanmean, results(2:end, row+2));
% results(nanIndices) = means;
% otherNumericVariables = cellfun(@double, results(2:end, 3:end));
% numericVariablesMean = mean(otherNumericVariables);
% numericVariablesStd = std(otherNumericVariables);
% [rows, ~] = size(results);
% results(rows + 1, 3:end) = num2cell(numericVariablesMean);
% results(rows + 2, 3:end) = num2cell(numericVariablesStd);
% results(rows + 1, 1:2) = {'Mean', 'NaN'};
% results(rows + 2, 1:2) = {'Std', 'NaN'};

smallVariationsColumns = strfind(results(1,:),'SH/NH');
smallVariationsColumns = find(~cellfun(@isempty,smallVariationsColumns));
morningVariations = cell2mat(results(2:end,smallVariationsColumns(1)));
eveningVariations = cell2mat(results(2:end,smallVariationsColumns(2)));

figure;
plot(decimalYears, eveningVariations, 'r.', 'MarkerSize', 10);
xlim([0 1]);
title('South std / North std of small variations')
ylabel('South std / North std');
xlabel('Decimal Year');
legend('Morning', 'Evening')
grid on

cell2csv('goceResults.csv', results);

end

function [morningTimestamps10s, morningMagneticLatitude, morningDensityNoBg, morningMsisDensity, eveningTimestamps10s, morningTimestampsDatenum, eveningTimetampsDatenum, ...
    eveningMagneticLatitude, eveningDensityNoBg, eveningMsisDensity] = ...
    splitBySolarTime(timestamps10s, timestampsDensityDatenum, magneticLatitude, densityNoBg, msisDensity, solarTime)
%

morningIndices = find(solarTime <= 12);
eveningIndices = find(solarTime > 12);
<<<<<<< HEAD

morningTimestamps10s = timestamps10s(morningIndices);
morningTimestampsDatenum = timestampsDensityDatenum(morningIndices);
morningMagneticLatitude = magneticLatitude(morningIndices);
morningDensityNoBg = densityNoBg(morningIndices);
morningMsisDensity = msisDensity(morningIndices);

=======

morningTimestamps10s = timestamps10s(morningIndices);
morningTimestampsDatenum = timestampsDensityDatenum(morningIndices);
morningMagneticLatitude = magneticLatitude(morningIndices);
morningDensityNoBg = densityNoBg(morningIndices);
morningMsisDensity = msisDensity(morningIndices);

>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
eveningTimestamps10s = timestamps10s(eveningIndices);
eveningTimetampsDatenum = timestampsDensityDatenum(eveningIndices);
eveningMagneticLatitude = magneticLatitude(eveningIndices);
eveningDensityNoBg = densityNoBg(eveningIndices);
eveningMsisDensity = msisDensity(eveningIndices);

end

function intervalsOfInterest = findInterestingIntervals(ae, timestampsAeDatenum, epsilonQualityFlag, timestampsEpsilonDatenum, densityNoBg, timestampsDensityDatenum, threshold)
%

intervalsOfInterest = zeros(1, 2);
medianCrossings = findCrossings(ae, 'smooth', 'mean');
calmDays = 2;
if length(medianCrossings) > 1
    for i = 1:length(medianCrossings) - 1
        if ~isempty(find(ae(medianCrossings(i):medianCrossings(i + 1)) >= threshold, 1))
            [intervalBegin, intervalEnd] = addCalmDaysToBothSidesOfPeak(medianCrossings(i), medianCrossings(i + 1), timestampsAeDatenum, calmDays);
            intervalsOfInterest = vertcat(intervalsOfInterest, [intervalBegin intervalEnd]);
        end
    end
end

intervalsOfInterest(1,:) = [];
if isempty(intervalsOfInterest)
    fprintf(2, '%s %d %s %d\n', 'There were no AE values above threshold: ', threshold, '. The maximum of given AE data is ', max(ae))
    error('No AE values above threshold found!')
<<<<<<< HEAD
end

indicesInDay = 24 * 60;
if length(intervalsOfInterest(:,1)) > 1
    lengthOfOverlap = 0;
    newIntervals = intervalsOfInterest;
    indicesToRemove = [];
    for i = 1:length(intervalsOfInterest(:,1)) - 1
        nextStormOverlap = intervalsOfInterest(i, 2) - intervalsOfInterest(i + 1,1);
        nextStormLength = intervalsOfInterest(i + 1, 2) - intervalsOfInterest(i + 1,1);
        if nextStormOverlap / nextStormLength > 0.4
            if lengthOfOverlap > 0
                newIntervals(i - lengthOfOverlap,:) = [newIntervals(i - lengthOfOverlap, 1) intervalsOfInterest(i + 1, 2)];
            else
                newIntervals(i,:) = [intervalsOfInterest(i, 1) intervalsOfInterest(i + 1, 2)];
            end
            indicesToRemove = [indicesToRemove; (i + 1)];
            lengthOfOverlap = lengthOfOverlap + 1;
        else
            lengthOfOverlap = 0;
        end
    end
    indicesToConserve = setdiff(1:length(newIntervals(:,1)), indicesToRemove);
    intervalsOfInterest = newIntervals(indicesToConserve,:);
end

intervalsToRemove = [];
for i = 1:length(intervalsOfInterest(:,1))
    aeIndices = intervalsOfInterest(i,1):intervalsOfInterest(i,2);
    [~, densityIntervalBeginIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(aeIndices(1))));
    [~, densityIntervalEndIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(aeIndices(end))));
    [~, epsilonIntervalBeginIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(aeIndices(1))));
    [~, epsilonIntervalEndIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(aeIndices(end))));
    epsilonIntervalIndices = epsilonIntervalBeginIndex:epsilonIntervalEndIndex;
    densityIntervalIndices = densityIntervalBeginIndex:densityIntervalEndIndex;
    
    stormBeginIndex = aeIndices(1) + calmDays * indicesInDay;
    stormEndIndex = aeIndices(end) - calmDays * indicesInDay;
    [~, densityStormBeginIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(stormBeginIndex)));
    [~, densityStormEndIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(stormEndIndex)));
    [~, epsilonStormBeginIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(stormBeginIndex)));
    [~, epsilonStormEndIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(stormEndIndex)));
    epsilonStormIndices = length(find(epsilonQualityFlag(epsilonStormBeginIndex:epsilonStormEndIndex)));
    densityStormIndices = length(find(ismember(timestampsDensityDatenum(densityStormBeginIndex:densityStormEndIndex), timestampsAeDatenum(stormBeginIndex:stormEndIndex))));
    epsilonBeforeIndices = length(find(epsilonQualityFlag(epsilonIntervalBeginIndex:epsilonStormBeginIndex)));
    densityBeforeIndices = length(find(ismember(timestampsDensityDatenum(densityIntervalBeginIndex:densityStormBeginIndex), timestampsAeDatenum(aeIndices(1):stormBeginIndex))));
    epsilonAfterIndices = length(find(epsilonQualityFlag(epsilonStormEndIndex:epsilonIntervalEndIndex)));
    densityAfterIndices = length(find(ismember(timestampsDensityDatenum(densityStormEndIndex:densityIntervalEndIndex), timestampsAeDatenum(stormEndIndex:aeIndices(end)))));
    
    valuesShouldHaveAfter = length(ae(stormEndIndex:aeIndices(end)));
    valuesShouldHaveBefore = length(ae(aeIndices(1):stormBeginIndex));
    valuesShouldHaveStorm = length(ae(stormBeginIndex:stormEndIndex));
    
    if densityStormIndices / valuesShouldHaveStorm < 0.5 ||...
       (densityAfterIndices / valuesShouldHaveAfter < 0.75 && densityBeforeIndices / valuesShouldHaveBefore < 0.75)
        intervalsToRemove = [intervalsToRemove i];
        fprintf(2, '%s\n', ['Warning: Storm between dates ', datestr(timestampsAeDatenum(aeIndices(1)), 'yyyy-mm-dd'), ...
            ' and ', datestr(timestampsAeDatenum(aeIndices(end)), 'yyyy-mm-dd'), ' has too large density data gaps. It will be omitted.']);
    end
    
    if epsilonStormIndices / valuesShouldHaveStorm < 0.5 ||...
       (epsilonAfterIndices / valuesShouldHaveAfter < 0.3 && epsilonBeforeIndices / valuesShouldHaveBefore < 0.3)
        intervalsToRemove = [intervalsToRemove i];
        fprintf(2, '%s\n', ['Warning: Storm between dates ', datestr(timestampsAeDatenum(aeIndices(1)), 'yyyy-mm-dd'), ...
            ' and ', datestr(timestampsAeDatenum(aeIndices(end)), 'yyyy-mm-dd'), ' has too large Akasofu epsilon data gaps. It will be omitted.']);
    end
      
end

intervalsToConserve = setdiff(1:length(intervalsOfInterest(:,1)), intervalsToRemove);
intervalsOfInterest = intervalsOfInterest(intervalsToConserve,:);

fprintf('%s\n', 'Following storms will be analyzed:')
for i = 1:length(intervalsOfInterest(:,1))
    beginIndex = intervalsOfInterest(i,1);
    endIndex = intervalsOfInterest(i,2);
    fprintf('%s\n', [datestr(timestampsAeDatenum(beginIndex), 'yyyy-mm-dd'), ' to ', datestr(timestampsAeDatenum(endIndex), 'yyyy-mm-dd')])
end

=======
end

indicesInDay = 24 * 60;
if length(intervalsOfInterest(:,1)) > 1
    lengthOfOverlap = 0;
    newIntervals = intervalsOfInterest;
    indicesToRemove = [];
    for i = 1:length(intervalsOfInterest(:,1)) - 1
        nextStormOverlap = intervalsOfInterest(i, 2) - intervalsOfInterest(i + 1,1);
        nextStormLength = intervalsOfInterest(i + 1, 2) - intervalsOfInterest(i + 1,1);
        if nextStormOverlap / nextStormLength > 0.4
            if lengthOfOverlap > 0
                newIntervals(i - lengthOfOverlap,:) = [newIntervals(i - lengthOfOverlap, 1) intervalsOfInterest(i + 1, 2)];
            else
                newIntervals(i,:) = [intervalsOfInterest(i, 1) intervalsOfInterest(i + 1, 2)];
            end
            indicesToRemove = [indicesToRemove; (i + 1)];
            lengthOfOverlap = lengthOfOverlap + 1;
        else
            lengthOfOverlap = 0;
        end
    end
    indicesToConserve = setdiff(1:length(newIntervals(:,1)), indicesToRemove);
    intervalsOfInterest = newIntervals(indicesToConserve,:);
end

intervalsToRemove = [];
for i = 1:length(intervalsOfInterest(:,1))
    aeIndices = intervalsOfInterest(i,1):intervalsOfInterest(i,2);
    [~, densityIntervalBeginIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(aeIndices(1))));
    [~, densityIntervalEndIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(aeIndices(end))));
    [~, epsilonIntervalBeginIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(aeIndices(1))));
    [~, epsilonIntervalEndIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(aeIndices(end))));
    epsilonIntervalIndices = epsilonIntervalBeginIndex:epsilonIntervalEndIndex;
    densityIntervalIndices = densityIntervalBeginIndex:densityIntervalEndIndex;
    
    stormBeginIndex = aeIndices(1) + calmDays * indicesInDay;
    stormEndIndex = aeIndices(end) - calmDays * indicesInDay;
    [~, densityStormBeginIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(stormBeginIndex)));
    [~, densityStormEndIndex] = min(abs(timestampsDensityDatenum - timestampsAeDatenum(stormEndIndex)));
    [~, epsilonStormBeginIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(stormBeginIndex)));
    [~, epsilonStormEndIndex] = min(abs(timestampsEpsilonDatenum - timestampsAeDatenum(stormEndIndex)));
    epsilonStormIndices = length(find(epsilonQualityFlag(epsilonStormBeginIndex:epsilonStormEndIndex)));
    densityStormIndices = length(find(ismember(timestampsDensityDatenum(densityStormBeginIndex:densityStormEndIndex), timestampsAeDatenum(stormBeginIndex:stormEndIndex))));
    epsilonBeforeIndices = length(find(epsilonQualityFlag(epsilonIntervalBeginIndex:epsilonStormBeginIndex)));
    densityBeforeIndices = length(find(ismember(timestampsDensityDatenum(densityIntervalBeginIndex:densityStormBeginIndex), timestampsAeDatenum(aeIndices(1):stormBeginIndex))));
    epsilonAfterIndices = length(find(epsilonQualityFlag(epsilonStormEndIndex:epsilonIntervalEndIndex)));
    densityAfterIndices = length(find(ismember(timestampsDensityDatenum(densityStormEndIndex:densityIntervalEndIndex), timestampsAeDatenum(stormEndIndex:aeIndices(end)))));
    
    valuesShouldHaveAfter = length(ae(stormEndIndex:aeIndices(end)));
    valuesShouldHaveBefore = length(ae(aeIndices(1):stormBeginIndex));
    valuesShouldHaveStorm = length(ae(stormBeginIndex:stormEndIndex));
    
    if densityStormIndices / valuesShouldHaveStorm < 0.5 ||...
       (densityAfterIndices / valuesShouldHaveAfter < 0.75 && densityBeforeIndices / valuesShouldHaveBefore < 0.75)
        intervalsToRemove = [intervalsToRemove i];
        fprintf(2, '%s\n', ['Warning: Storm between dates ', datestr(timestampsAeDatenum(aeIndices(1)), 'yyyy-mm-dd'), ...
            ' and ', datestr(timestampsAeDatenum(aeIndices(end)), 'yyyy-mm-dd'), ' has too large density data gaps. It will be omitted.']);
    end
    
    if epsilonStormIndices / valuesShouldHaveStorm < 0.5 ||...
       (epsilonAfterIndices / valuesShouldHaveAfter < 0.3 && epsilonBeforeIndices / valuesShouldHaveBefore < 0.3)
        intervalsToRemove = [intervalsToRemove i];
        fprintf(2, '%s\n', ['Warning: Storm between dates ', datestr(timestampsAeDatenum(aeIndices(1)), 'yyyy-mm-dd'), ...
            ' and ', datestr(timestampsAeDatenum(aeIndices(end)), 'yyyy-mm-dd'), ' has too large Akasofu epsilon data gaps. It will be omitted.']);
    end
      
end

intervalsToConserve = setdiff(1:length(intervalsOfInterest(:,1)), intervalsToRemove);
intervalsOfInterest = intervalsOfInterest(intervalsToConserve,:);

fprintf('%s\n', 'Following storms will be analyzed:')
for i = 1:length(intervalsOfInterest(:,1))
    beginIndex = intervalsOfInterest(i,1);
    endIndex = intervalsOfInterest(i,2);
    fprintf('%s\n', [datestr(timestampsAeDatenum(beginIndex), 'yyyy-mm-dd'), ' to ', datestr(timestampsAeDatenum(endIndex), 'yyyy-mm-dd')])
end

>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
end

function [intervalBegin, intervalEnd] = addCalmDaysToBothSidesOfPeak(peakBeginIndex, peakEndIndex, timestampsAeDatenum, calmDays)
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

function [ae, ap, absB, akasofuEpsilon, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, morningTimestampsDatenum, eveningTimetampsDatenum, timestamps1minOut, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = sliceToInterestingIntervals(ae, ap, absB, averagedDensity, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, morningTimestampsDatenum, eveningTimetampsDatenum, timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed,...
 density3h, morningMagneticLatitude, eveningMagneticLatitude, timestamps10sFixed, timestamps1min, timestampsAbsB, akasofuEpsilon, timestampsDatenum, intervalsOfInterest)
%

aeTemp = ae; apTemp = ap; averagedDensityTemp = averagedDensity;
averagedDensityNoBgTemp = averagedDensityNoBg; morningDensityNoBgTemp = morningDensityNoBg; 
eveningDensityNoBgTemp = eveningDensityNoBg; morningMsisDensityTemp = morningMsisDensity;
eveningMsisDensityTemp = eveningMsisDensity; morningTimestamps10sTemp = morningTimestamps10s; 
eveningTimestamps10sTemp = eveningTimestamps10s; timestamps1minFixedTemp = timestamps1minFixed; 
morningTimestampsDatenumTemp = morningTimestampsDatenum; eveningTimestampsDatenumTemp = eveningTimetampsDatenum;
timestamps3hTemp = timestamps3h; timestamps3hFixedTemp = timestamps3hFixed; density3hTemp = density3h; morningMagneticLatitudeTemp = morningMagneticLatitude;
eveningMagneticLatitudeTemp = eveningMagneticLatitude; timestampsEpsilonTemp = timestampsEpsilon; akasofuEpsilonTemp = akasofuEpsilon;
timestampsAbsBTemp = timestampsAbsB; absBTemp = absB; timestampsDatenumTemp = timestampsDatenum; 

ae = {}; ap = {}; averagedDensity = {}; averagedDensityNoBg = {}; morningDensityNoBg = {};
eveningDensityNoBg = {}; morningMsisDensity = {}; eveningMsisDensity = {}; morningTimestamps10s = {}; eveningTimetampsDatenum = {};
eveningTimestamps10s = {}; morningTimestampsDatenum = {}; timestamps1minFixed = {}; timestamps3h = {}; density3h = {}; morningMagneticLatitude = {};
eveningMagneticLatitude = {}; timestampsAbsB = {}; akasofuEpsilon = {}; timestampsEpsilon = {}; absB = {}; timestampsDatenum = {};
timestamps3hFixed = {};

cellArrayLength = length(intervalsOfInterest(:,1));

for i = 1:cellArrayLength
    beginIndex = intervalsOfInterest(i,1);
    endIndex = intervalsOfInterest(i,2);
    timestamps1minOut{i} = timestamps1min(beginIndex:endIndex);
    ae{i} = aeTemp(beginIndex:endIndex);
    
    threeHinSec = 3 * 60 * 60;
    [~, beginIndex3h] = min(abs(timestamps3hTemp - threeHinSec * round(timestamps1min(beginIndex) / threeHinSec)));
    [~, endIndex3h] = min(abs(timestamps3hTemp - threeHinSec * round(timestamps1min(endIndex) / threeHinSec)));
    ap{i} = apTemp(beginIndex3h:endIndex3h);
    timestamps3h{i} = timestamps3hTemp(beginIndex3h:endIndex3h);
    
    [~, beginIndex3hFixed] = min(abs(timestamps3hFixedTemp - threeHinSec * round(timestamps1min(beginIndex) / threeHinSec)));
    [~, endIndex3hFixed] = min(abs(timestamps3hFixedTemp - threeHinSec * round(timestamps1min(endIndex) / threeHinSec)));
    density3h{i} = density3hTemp(beginIndex3hFixed:endIndex3hFixed);
    timestamps3hFixed{i} = timestamps3hFixedTemp(beginIndex3hFixed:endIndex3hFixed);
    
    [~, beginIndexAbsB] = min(abs(timestampsAbsBTemp - timestamps1min(beginIndex)));
    [~, endIndexAbsB] = min(abs(timestampsAbsBTemp - timestamps1min(endIndex)));
    absB{i} = absBTemp(beginIndexAbsB:endIndexAbsB);
    timestampsAbsB{i} = timestampsAbsBTemp(beginIndexAbsB:endIndexAbsB);
    
    [~, beginIndexEpsilon] = min(abs(timestampsEpsilonTemp - timestamps1min(beginIndex)));
    [~, endIndexEpsilon] = min(abs(timestampsEpsilonTemp - timestamps1min(endIndex)));
    akasofuEpsilon{i} = akasofuEpsilonTemp(beginIndexEpsilon:endIndexEpsilon);
    timestampsEpsilon{i} = timestampsEpsilonTemp(beginIndexEpsilon:endIndexEpsilon);
    
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
    morningTimestampsDatenum{i} = morningTimestampsDatenumTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningTimestampsDatenum{i} = eveningTimestampsDatenumTemp(beginIndexEvening10s:endIndexEvening10s);
    morningMagneticLatitude{i} = morningMagneticLatitudeTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningMagneticLatitude{i} = eveningMagneticLatitudeTemp(beginIndexEvening10s:endIndexEvening10s);
    
    [~, beginIndex10s] = min(abs(timestamps10sFixed - timestamps1min(beginIndex)));
    [~, endIndex10s] = min(abs(timestamps10sFixed - timestamps1min(endIndex)));
    timestampsDatenum{i} = timestampsDatenumTemp(beginIndex10s:endIndex10s);
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
<<<<<<< HEAD

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
    if maxIndex < length(medianCrossings)
        peakEnd = medianCrossings(maxIndex + 1);
    else
        peakEnd = length(subtractedVal);
    end
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

=======

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
    if maxIndex < length(medianCrossings)
        peakEnd = medianCrossings(maxIndex + 1);
    else
        peakEnd = length(subtractedVal);
    end
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

>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
crossings = find(subtractedVal(1:end-1) .* subtractedVal(2:end) < 0);

end

function plotTimeseries(timestamps1min, timestamps1minFixed, timestampsAbsB, timestamps3h, timestamps3hFixed, ae, ap, absB, ...
    averagedDensityNoBg, density3h)
% plotTimeseries(timestamps, ae, averagedDensity)

global timeseriesFigHandle
timeseriesFigHandle = figure;
subplot(2,2,1)
secondsInDay = 60 * 60 * 24;
beginDay = datenum('2009-11-01', 'yyyy-mm-dd');
timestampsInDays1min = timestamps1min / secondsInDay + beginDay;
timestampsInDays1minFixed = timestamps1minFixed / secondsInDay + beginDay;
timestampsInDaysAbsB = timestampsAbsB / secondsInDay + beginDay;
[hAx,~,~] = plotyy(timestampsInDaysAbsB, absB, timestampsInDays1minFixed, averagedDensityNoBg);
title('Timeseries of IMF |B| and density')
datetick(hAx(1), 'x', 'yyyy-mm-dd', 'keepticks', 'keeplimits')
datetick(hAx(2), 'x', 'yyyy-mm-dd', 'keepticks', 'keeplimits')
rotateticklabel(hAx(1), 50);
rotateticklabel(hAx(2), 50);
ylabel(hAx(1), 'IMF |B| / nT')
ylabel(hAx(2), 'Density')
set(hAx, 'XLim', [min(timestampsInDays1min) max(timestampsInDays1min)]);
grid on

subplot(2,2,3)
timestampsInDays3h = timestamps3h / secondsInDay + beginDay;
timestampsInDays3hFixed = timestamps3hFixed / secondsInDay + beginDay;
[hAx,~,~] = plotyy(timestampsInDays3h, ap, timestampsInDays3hFixed, density3h);
title('Timeseries of ap and density')
datetick(hAx(1), 'x', 'yyyy-mm-dd', 'keepticks', 'keeplimits')
datetick(hAx(2), 'x', 'yyyy-mm-dd', 'keepticks', 'keeplimits')
rotateticklabel(hAx(1), 50);
rotateticklabel(hAx(2), 50);
xlabel('t / days')
ylabel(hAx(1), 'ap')
ylabel(hAx(2), 'previous 3h average density')
set(hAx, 'XLim', [min(timestampsInDays3h) max(timestampsInDays3h)]);
grid on

subplot(2,2,4)
[hAx,~,~] = plotyy(timestampsInDays1min, ae, timestampsInDays1minFixed, averagedDensityNoBg);
title('Timeseries of AE and density')
datetick(hAx(1), 'x', 'yyyy-mm-dd', 'keepticks', 'keeplimits')
datetick(hAx(2), 'x', 'yyyy-mm-dd', 'keepticks', 'keeplimits')
rotateticklabel(hAx(1), 50);
rotateticklabel(hAx(2), 50);
xlabel('t / days')
ylabel(hAx(1), 'AE')
ylabel(hAx(2), 'Density')
set(hAx, 'XLim', [min(timestampsInDays1min) max(timestampsInDays1min)]);
grid on


end

function [densityIndexTimelag] = giveMaxCrossCorrelation(density, geomIndex)
% plotCrossCorrelation(averagedDensityNoBg, ae)

maxLag = 60 * 24;
[correlations, lags] = xcorr(density, geomIndex, maxLag, 'coeff');
correlations = correlations(lags > 0);
lags = lags(lags > 0);
indicesInHour = 60;
densityIndexTimelag = lags(correlations == max(correlations)) / indicesInHour;
<<<<<<< HEAD

end

=======

end

>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
function results = plotAndCalculateCorrelation(timestamps, timestampsFixed, geomIndex, density, indexName, plotFigures, results)
% [r, r2] = plotAndCalculateCorrelation(ae, averagedDensity, timelag)

global timeseriesFigHandle
if exist('timeseriesFigHandle', 'var')
    timeseriesHandle = timeseriesFigHandle;
end

[timelag, timelagInHours, bestIntegral, averageGoodLag, averageIntegral] = compareDensityToGeomIndexIntegral(density, geomIndex, timestamps, timestampsFixed, indexName, plotFigures);

[rowNum, ~] = size(results);
emptyCells = cellfun(@isempty,results);
[~, emptyColPositions] = find(emptyCells);
colNum = min(emptyColPositions);
if rowNum == 2
    colNum = length(results(rowNum,:)) + 1;
    if strcmpi(indexName, 'ae')
        results{1, colNum} = 'Max Aver. AE Int.';
        results{1, colNum + 1} = ['Int window: ',indexName];
    else
        results{1, colNum} = ['Int window: ',indexName];
    end
end

if strcmpi(indexName, 'ae')
    results{rowNum, colNum} = max(averageIntegral);
    results{rowNum, colNum + 1} = timelagInHours;
else
    results{rowNum, colNum} = timelagInHours;
end

timestampsFixed = timestampsFixed(ismember(timestampsFixed, timestamps));
densityBestIntIndices = ismember(timestampsFixed, timestamps(timelag + 1:end));
    
if strcmpi(indexName, 'ae') && plotFigures ~=0
    figure(timeseriesHandle);
    subplot(2,2,2)
    secondsInDay = 24 * 60 * 60;
    beginDay = datenum('2009-11-01', 'yyyy-mm-dd');
    timestampsInDays = timestamps(timelag + 1:end) / secondsInDay + beginDay;
    timestampsInDaysFixed = timestampsFixed(densityBestIntIndices) / secondsInDay + beginDay;
    [hAx,~,~] = plotyy(timestampsInDays, bestIntegral, timestampsInDaysFixed, density(densityBestIntIndices));
    xlabel('t / days')
    ylabel(hAx(1), 'AE Integral')
    ylabel(hAx(2), 'Density')
    title('Raw AE integral vs FFT smoothed density');
    set(hAx, 'XLim', [min(timestampsInDaysFixed) max(timestampsInDaysFixed)]);    
    grid on;
end

bestIntegral = bestIntegral(ismember(timestamps(timelag + 1:end), timestampsFixed));
results = plotCorrelation(bestIntegral, density(densityBestIntIndices), [indexName, ' Best Integral'], 'Density', plotFigures, results);

densityAverageIntIndices = ismember(timestampsFixed, timestamps(averageGoodLag + 1:end));
averageIntegral = averageIntegral(ismember(timestamps(averageGoodLag + 1:end), timestampsFixed));
results = plotCorrelation(averageIntegral, density(densityAverageIntIndices), [indexName, ' Average Integral'], 'Density', plotFigures, results);

geomIndexFixed = geomIndex(ismember(timestamps, timestampsFixed));
densityFixed = density(ismember(timestampsFixed, timestamps));
if strcmpi(indexName, 'ae')
    geomIndexNoBg = removePeriodicBackground(geomIndexFixed, 125, 1, 0);
    geomIndexNoBg = normalize(geomIndexNoBg, geomIndexFixed);
    results = plotCorrelation(geomIndexNoBg, densityFixed, indexName, 'Density at 270 km', plotFigures, results);
elseif strcmpi(indexName, 'Akasofu Epsilon') || ~isempty(strfind(upper(indexName), '|B|'))
    timestamps6hAgo = timestamps(ismember(timestamps, timestamps - 6 * 60 * 60));
    timestampsShorter = timestamps6hAgo + 6 * 60 * 60;
    geomIndex6hAgo = geomIndex(ismember(timestamps, timestamps6hAgo));
    geomIndex6hAgo = geomIndex6hAgo(ismember(timestampsShorter, timestampsFixed));
    densityShorter = density(ismember(timestampsFixed, timestampsShorter));
    results = plotCorrelation(geomIndex6hAgo, densityShorter, indexName, 'Density at 270 km', plotFigures, results);
else
    results = plotCorrelation(geomIndexFixed, densityFixed, indexName, 'Density at 270 km', plotFigures, results);
end

end

function results = plotCorrelation(xvals, yvals, xvalName, yvalName, plotFigures, results)
%

r = corr(xvals, yvals);

[rowNum, ~] = size(results);
if rowNum > 0
    emptyCells = cellfun(@isempty,results);
    [~, emptyColPositions] = find(emptyCells);
    colNum = min(emptyColPositions);
    if rowNum == 2
        colNum = length(results(rowNum,:)) + 1;
        results{1, colNum} = ['r: ',xvalName];
    end
    results{rowNum, colNum} = r;
end

if plotFigures ~= 0
    r2 = r * r;
    fprintf('%s %f\n', 'Pearson correlation: ', r)
    fprintf('%s %d %s\n', 'Thus, ', round(r2 * 100), ['% of variation in ', yvalName, ' can be explained by changes in ', xvalName])
    
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
<<<<<<< HEAD

end

function [bestLag, bestLagInHours, geomIndexBestInt, averageGoodLag, averageLagInt] = compareDensityToGeomIndexIntegral(density, ...
    geomIndex, timestamps, timestampsFixed, indexName, plotFigures)
%
maxDays = 3;
if strcmpi(indexName, 'ap'); maxLag = maxDays * 8;
else maxLag = maxDays * 24 * 60; end

indicesInHour = 60;
if strcmpi(indexName, 'ae')
    averageGoodLag = 22 * indicesInHour;
elseif strcmpi(indexName, 'Akasofu Epsilon')
    averageGoodLag = 38 * indicesInHour;
elseif ~isempty(strfind(upper(indexName), '|B|'))
    averageGoodLag = 34 * indicesInHour;
else
    averageGoodLag = 8;
end

cumulativeGeomIndex = cumsum(geomIndex);
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
if strcmpi(indexName, 'ap'); lagsInHours = lags * 3;
else lagsInHours = lags / 60; end

bestLag = lags(find(correlations == max(correlations), 1, 'last'));
bestLagInHours = lagsInHours(find(correlations == max(correlations), 1, 'last'));
if bestLag > 0
    geomIndexBestInt = cumulativeGeomIndex(bestLag + 1 : end) - cumulativeGeomIndex(1 : end - bestLag);
else
    geomIndexBestInt = geomIndex;
end
averageLagInt = cumulativeGeomIndex(averageGoodLag + 1 : end) - cumulativeGeomIndex(1 : end - averageGoodLag);

if plotFigures ~= 0
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
end
% if strcmpi(indexName, 'ae'); integralWindowSize = integralWindowSize * 60;
% elseif ~isempty(strfind(upper(indexName), '|B|')); integralWindowSize = integralWindowSize * 60 / 4;
% else integralWindowSize = integralWindowSize / 3;end

end

function plotDensityLatitudeTimeSurf(correctedDensity, msisDensity, magneticLatitude, timestamps10s, timeOfDay)
=======

end

function [bestLag, bestLagInHours, geomIndexBestInt, averageGoodLag, averageLagInt] = compareDensityToGeomIndexIntegral(density, ...
    geomIndex, timestamps, timestampsFixed, indexName, plotFigures)
%
maxDays = 3;
if strcmpi(indexName, 'ap'); maxLag = maxDays * 8;
else maxLag = maxDays * 24 * 60; end

indicesInHour = 60;
if strcmpi(indexName, 'ae')
    averageGoodLag = 22 * indicesInHour;
elseif strcmpi(indexName, 'Akasofu Epsilon')
    averageGoodLag = 38 * indicesInHour;
elseif ~isempty(strfind(upper(indexName), '|B|'))
    averageGoodLag = 34 * indicesInHour;
else
    averageGoodLag = 8;
end

cumulativeGeomIndex = cumsum(geomIndex);
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
if strcmpi(indexName, 'ap'); lagsInHours = lags * 3;
else lagsInHours = lags / 60; end

bestLag = lags(find(correlations == max(correlations), 1, 'last'));
bestLagInHours = lagsInHours(find(correlations == max(correlations), 1, 'last'));
if bestLag > 0
    geomIndexBestInt = cumulativeGeomIndex(bestLag + 1 : end) - cumulativeGeomIndex(1 : end - bestLag);
else
    geomIndexBestInt = geomIndex;
end
averageLagInt = cumulativeGeomIndex(averageGoodLag + 1 : end) - cumulativeGeomIndex(1 : end - averageGoodLag);

if plotFigures ~= 0
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
end
% if strcmpi(indexName, 'ae'); integralWindowSize = integralWindowSize * 60;
% elseif ~isempty(strfind(upper(indexName), '|B|')); integralWindowSize = integralWindowSize * 60 / 4;
% else integralWindowSize = integralWindowSize / 3;end

end

function plotDensityLatitudeTimeSurf(correctedDensity, msisDensity, magneticLatitude, timestamps10s, regriddedLatitude, regriddedTime, regriddedDensity, timeOfDay)
>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
% plotDensityLatitudeTimeSurf(averagedDensity, averagedLatitude, timestamps)

figure;
secondsInDay = 60 * 60 * 24;
<<<<<<< HEAD
timestampsInDays = timestamps10s / secondsInDay + datenum('2009-11-01', 'yyyy-mm-dd');
[minLat, maxLat] = findInterpolationLimits(magneticLatitude);
[tInterp, latitudeInterp] = meshgrid(timestampsInDays(1):9500/secondsInDay:timestampsInDays(end), minLat:1:maxLat);

% F = scatteredInterpolant(timestampsInDays, magneticLatitude, correctedDensity);
densityInterp = griddata(timestampsInDays, magneticLatitude, correctedDensity, tInterp, latitudeInterp, 'v4');
% densityInterp = F(tInterp, latitudeInterp);
surf(tInterp, latitudeInterp, densityInterp, 'EdgeColor', 'None')
=======
timestampsInDays = timestamps10s; %/ secondsInDay + datenum('2009-11-01', 'yyyy-mm-dd');
[minLat, maxLat] = findInterpolationLimits(magneticLatitude);

%densityInterp = griddata(timestampsInDays, magneticLatitude, correctedDensity, regriddedTime, regriddedLatitude, 'linear');
surf(regriddedTime, regriddedLatitude, regriddedDensity, 'EdgeColor', 'None')
%surf(regriddedTime, regriddedLatitude, densityInterp, 'EdgeColor', 'None')
>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb

xlim([timestampsInDays(1) timestampsInDays(end)]);
ylim([minLat maxLat]);
caxis([min(correctedDensity) max(correctedDensity)])
colormap jet(30)
colorbar
view(2);
xlabel('t / days')
ylabel('Lat (°)')
zlabel('Density')

end

<<<<<<< HEAD
function results = plotAndAnalyzeDensityByLatitude(ae, timestamps1min, timestamps1minFixed, correctedDensity, ...
=======
function [results, latitudeMatrix, timeMatrix, densityMatrix] = plotAndAnalyzeDensityByLatitude(ae, timestamps1min, timestamps1minFixed, correctedDensity, ...
>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
    timestamps10s, magneticLatitude, timeOfDay, plotFigures, results)
% plotCorrectedDensityLatitudes(ae, timestamps1min, correctedDensity, timestamps10s, latitude, timestampsDatenum, computeLatitudes);
    
[latBeginIndex, latEndIndex] = limitLatitudeToIntegerMultipleOfOrbitalPeriod(magneticLatitude);
limitedLatitude = magneticLatitude(latBeginIndex:latEndIndex);
limitedTimestamps = timestamps10s(latBeginIndex:latEndIndex);
<<<<<<< HEAD
[~, exactOrbitIndices] = splitIntoOrbits(limitedLatitude);
limitedLatitude = limitedLatitude(exactOrbitIndices);
limitedTimestamps = limitedTimestamps(exactOrbitIndices);
[orbits, ~] = splitIntoOrbits(limitedLatitude);
=======
[~, exactOrbitIndices] = splitIntoOrbits(limitedLatitude, limitedTimestamps);
limitedLatitude = limitedLatitude(exactOrbitIndices);
limitedTimestamps = limitedTimestamps(exactOrbitIndices);
[orbits, ~] = splitIntoOrbits(limitedLatitude, limitedTimestamps);
>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb

[minAllowedLatitude, maxAllowedLatitude] = findInterpolationLimits(limitedLatitude);
orbitsToDelete = find(abs(limitedLatitude(orbits(:,2)) - limitedLatitude(orbits(:,1))) < ...
    (maxAllowedLatitude - minAllowedLatitude));
newIndices = 1:length(limitedLatitude);
for i = 1:length(orbitsToDelete)
    newIndices = setdiff(newIndices, (orbits(orbitsToDelete(i), 1) : orbits(orbitsToDelete(i), 2)));
end
limitedTimestamps = limitedTimestamps(newIndices);
limitedLatitude = limitedLatitude(newIndices);

<<<<<<< HEAD
F = scatteredInterpolant(timestamps10s, magneticLatitude, correctedDensity);

oneDegreeStep = minAllowedLatitude:1:maxAllowedLatitude;

parfor i = 1:length(oneDegreeStep)
    %oneDegreeStep(i) = checkOutOfBounds(minAllowedLatitude, maxAllowedLatitude, oneDegreeStep(i));
    crossingTimes(:,i) = latitudeCrossingTimes(limitedLatitude, limitedTimestamps, oneDegreeStep(i));

    latitudeInterp = ones(length(crossingTimes(:,i)), 1) * oneDegreeStep(i);    
    densityByLatitude(:,i) = F(crossingTimes(:,i), latitudeInterp);
=======
oneDegreeStep = minAllowedLatitude:maxAllowedLatitude;
oneQuarterDegreeStep = minAllowedLatitude:0.25:maxAllowedLatitude;

parfor i = 1:length(oneQuarterDegreeStep)
    regriddedTime(:,i) = latitudeCrossingTimes(limitedLatitude, limitedTimestamps, oneQuarterDegreeStep(i));    
end

% regriddedDensity = griddata(timestamps10s, magneticLatitude, correctedDensity, regriddedTime, latitudeMatrix, 'natural');
regriddedDensity = interp1(timestamps10s, correctedDensity, regriddedTime, 'linear');
crossingTimes = regriddedTime(:,1:4:end);
densityByLatitude = regriddedDensity(:,1:4:end);

numOfOrbits = length(regriddedTime(:,1));
numOfValuesInOrbit = length(regriddedTime(1,:));
for i = 1:numOfValuesInOrbit
    timeThisLatitude = regriddedTime(:,i);
    densityThisLatitude = regriddedDensity(:,i);
    
    tInterp = interp1(1:numOfOrbits, timeThisLatitude, 1:1/8:numOfOrbits);
    interpolatedDensity = interp1(timeThisLatitude, densityThisLatitude, tInterp, 'linear');
    
    latitudeMatrix(:,i) = ones(length(tInterp), 1) * oneQuarterDegreeStep(i); 
    densityMatrix(:,i) = interpolatedDensity;
    timeMatrix(:,i) = tInterp;
>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
end

limitedTimestamps = limitedTimestamps(ismember(limitedTimestamps, timestamps1minFixed));

aeShort = ae(ismember(timestamps1min, limitedTimestamps));
if plotFigures ~= 0
    writeAndPlotPeakAnalysis(limitedTimestamps, aeShort, densityByLatitude, crossingTimes, minAllowedLatitude, maxAllowedLatitude, timeOfDay);
end

northIndices = (oneDegreeStep < 90 & oneDegreeStep > 45);
equatorIndices = (oneDegreeStep < 30 & oneDegreeStep > -30);
southIndices = (oneDegreeStep < -45 & oneDegreeStep > -90);

northernDensity = mean(densityByLatitude(:,northIndices), 2);
equatorDensity = mean(densityByLatitude(:,equatorIndices), 2);
southernDensity = mean(densityByLatitude(:,southIndices), 2);

northernError = std(densityByLatitude(:,northIndices), 0, 2);
equatorError = std(densityByLatitude(:,equatorIndices), 0, 2);
southernError = std(densityByLatitude(:,southIndices), 0, 2);

beginDay = datenum('2009-11-01', 'yyyy-mm-dd');

northTimestamps = mean(crossingTimes(:,northIndices), 2);
equatorTimestamps = mean(crossingTimes(:,equatorIndices), 2);
southTimestamps = mean(crossingTimes(:,southIndices), 2);

[~, northPeakBegin, northPeakEnd] = limitToNearPeak(northernDensity, 'noSmooth', 'median');
[~, equatorPeakBegin, equatorPeakEnd] = limitToNearPeak(equatorDensity, 'noSmooth', 'median');
[~, southPeakBegin, southPeakEnd] = limitToNearPeak(southernDensity, 'noSmooth', 'median');

northCalmMean = mean(northernDensity(1:northPeakBegin));
equatorCalmMean = mean(equatorDensity(1:equatorPeakBegin));
southCalmMean = mean(southernDensity(1:southPeakBegin));

northAbsDiff = northernDensity - northCalmMean;
equatorAbsDiff = equatorDensity - equatorCalmMean;
southAbsDiff = southernDensity - southCalmMean;

northRelDiff = northAbsDiff / northCalmMean;
equatorRelDiff = equatorAbsDiff / equatorCalmMean;
southRelDiff = southAbsDiff / southCalmMean;

meanNorthAbsDiff = mean(northAbsDiff(northPeakBegin:northPeakEnd));
meanEquatorAbsDiff = mean(equatorAbsDiff(equatorPeakBegin:equatorPeakEnd));
meanSouthAbsDiff = mean(southAbsDiff(southPeakBegin:southPeakEnd));
meanNorthRelDiff = mean(northRelDiff(northPeakBegin:northPeakEnd));
meanEquatorRelDiff = mean(equatorRelDiff(equatorPeakBegin:equatorPeakEnd));
meanSouthRelDiff = mean(southRelDiff(southPeakBegin:southPeakEnd));

maxNorthAbsDiff = max(northAbsDiff);
maxEquatorAbsDiff = max(equatorAbsDiff);
maxSouthAbsDiff = max(southAbsDiff);
maxNorthRelDiff = max(northRelDiff);
maxEquatorRelDiff = max(equatorRelDiff);
maxSouthRelDiff = max(southRelDiff);

northRelDiffForXcorr = interp1(northTimestamps, northRelDiff, timestamps1min, 'nearest', 'extrap');
equatorRelDiffForXcorr = interp1(equatorTimestamps, equatorRelDiff, timestamps1min, 'nearest', 'extrap');
southRelDiffForXcorr = interp1(southTimestamps, southRelDiff, timestamps1min, 'nearest', 'extrap');

northLag = giveMaxCrossCorrelation(northRelDiffForXcorr, ae);
equatorLag = giveMaxCrossCorrelation(equatorRelDiffForXcorr, ae);
southLag = giveMaxCrossCorrelation(southRelDiffForXcorr, ae);

[rowNum, ~] = size(results);
emptyCells = cellfun(@isempty,results);
[~, emptyColPositions] = find(emptyCells);
colNum = min(emptyColPositions);
if rowNum == 2
    colNum = length(results(rowNum,:)) + 1;
    results{1, colNum}     = ['NH AE-Dens. ', timeOfDay, ' lag'];
    results{1, colNum + 1} = ['EQ AE-Dens. ', timeOfDay, ' lag'];
    results{1, colNum + 2} = ['SH AE-Dens. ', timeOfDay, ' lag'];
    
    results{1, colNum + 3} = ['Mean NH Abs. Diff ', timeOfDay];
    results{1, colNum + 4} = ['Mean EQ Abs. Diff ', timeOfDay];
    results{1, colNum + 5} = ['Mean SH Abs. Diff ', timeOfDay];
    results{1, colNum + 6} = ['Mean NH Rel. Diff ', timeOfDay];
    results{1, colNum + 7} = ['Mean EQ Rel. Diff ', timeOfDay];
    results{1, colNum + 8} = ['Mean SH Rel. Diff ', timeOfDay];
    
    results{1, colNum + 9} = ['Max NH Abs. Diff ', timeOfDay];
    results{1, colNum + 10} = ['Max EQ Abs. Diff ', timeOfDay];
    results{1, colNum + 11} = ['Max SH Abs. Diff ', timeOfDay];
    results{1, colNum + 12} = ['Max NH Rel. Diff ', timeOfDay];
    results{1, colNum + 13} = ['Max EQ Rel. Diff ', timeOfDay];
    results{1, colNum + 14} = ['Max SH Rel. Diff ', timeOfDay];
end
results{rowNum, colNum} =     northLag;
results{rowNum, colNum + 1} = equatorLag;
results{rowNum, colNum + 2} = southLag;

results{rowNum, colNum + 3} = meanNorthAbsDiff;
results{rowNum, colNum + 4} = meanEquatorAbsDiff;
results{rowNum, colNum + 5} = meanSouthAbsDiff;
results{rowNum, colNum + 6} = meanNorthRelDiff;
results{rowNum, colNum + 7} = meanEquatorRelDiff;
results{rowNum, colNum + 8} = meanSouthRelDiff;

results{rowNum, colNum + 9} = maxNorthAbsDiff;
results{rowNum, colNum + 10} = maxEquatorAbsDiff;
results{rowNum, colNum + 11} = maxSouthAbsDiff;
results{rowNum, colNum + 12} = maxNorthRelDiff;
results{rowNum, colNum + 13} = maxEquatorRelDiff;
results{rowNum, colNum + 14} = maxSouthRelDiff;

if plotFigures ~= 0
    secondsInDay = 24 * 60 * 60;

    northTimestamps = northTimestamps / secondsInDay + beginDay;
    equatorTimestamps = equatorTimestamps / secondsInDay + beginDay;
    southTimestamps = southTimestamps / secondsInDay + beginDay;
    
    figure;
    hAx(1) = subplot(3,2,1);
    errorbar(northTimestamps, northAbsDiff, northernError);
    datetick('x', 'dd', 'keepticks', 'keeplimits')
    ylabel('Density [10^{-11} kgm^{-3}]')
    title('Absolute difference, NH')
    grid on

    hAx(3) = subplot(3,2,3);
    errorbar(equatorTimestamps, equatorAbsDiff, equatorError, 'g');
    datetick('x', 'dd', 'keepticks', 'keeplimits')
    ylabel('Density [10^{-11} kgm^{-3}]')
    title('Absolute difference, equator')
    grid on

    hAx(5) = subplot(3,2,5);
    errorbar(southTimestamps, southAbsDiff, southernError, 'r');
    datetick('x', 'dd', 'keepticks', 'keeplimits')
    ylabel('Density [10^{-11} kgm^{-3}]')
    title('Absolute difference, SH')
    grid on

    hAx(2) = subplot(3,2,2);
    errorbar(northTimestamps, northRelDiff * 100, 100 * northernError/northCalmMean);
    datetick('x', 'dd', 'keepticks', 'keeplimits')
    ylabel('%')
    title('Relative difference, NH')
    grid on

    hAx(4) = subplot(3,2,4);
    errorbar(equatorTimestamps, equatorRelDiff * 100, 100 * equatorError/equatorCalmMean, 'g');
    datetick('x', 'dd', 'keepticks', 'keeplimits')
    ylabel('%')
    title('Relative difference, equator')
    grid on

    hAx(6) = subplot(3,2,6);
    errorbar(southTimestamps, southRelDiff * 100, 100 * southernError/southCalmMean, 'r');
    datetick('x', 'dd', 'keepticks', 'keeplimits')
    ylabel('%')
    title('Relative difference, SH')
    grid on

    annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Absolute and relative density changes on hemispheres: ', timeOfDay], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
end

end

<<<<<<< HEAD
function [orbits, exactOrbitIndices] = splitIntoOrbits(latitude)
%

if satelliteIsGoingSouth(latitude)
    orbitEndIndices = [find(latitude(1:end-1) < latitude(2:end)); length(latitude)];
else
    orbitEndIndices = [find(latitude(1:end-1) > latitude(2:end)); length(latitude)];
=======
function [orbits, exactOrbitIndices] = splitIntoOrbits(latitude, timestamps10s)
%
oneHinSec = 60 * 60;
if satelliteIsGoingSouth(latitude)
    orbitEndIndices = [find(latitude(1:end-1) < latitude(2:end) | timestamps10s(1:end -1) - timestamps10s(2:end) < -1 * oneHinSec); length(latitude)];
else
    orbitEndIndices = [find(latitude(1:end-1) > latitude(2:end) | timestamps10s(1:end -1) - timestamps10s(2:end) < -1 * oneHinSec); length(latitude)];
>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
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
<<<<<<< HEAD

function writeAndPlotPeakAnalysis(limitedTimestamps, ae, densityByLatitude, crossingTimes, minAllowedLatitude, maxAllowedLatitude, timeOfDay)
% writeAndPlotPeakAnalysis(timestamps, ae, splinedDensity, computeLatitutes)

persistent aeDensityLagByLatHandle
%global densityPeakMaxTimesByLatHandle

=======

function writeAndPlotPeakAnalysis(limitedTimestamps, ae, densityByLatitude, crossingTimes, minAllowedLatitude, maxAllowedLatitude, timeOfDay)
% writeAndPlotPeakAnalysis(timestamps, ae, splinedDensity, computeLatitutes)

persistent aeDensityLagByLatHandle
%global densityPeakMaxTimesByLatHandle

>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
ae = removePeriodicBackground(ae, 432, 1, 0); % 432 min =^ 0.3 days
% ae = smooth(ae, 431);
[ae, ~, ~] = limitToNearPeak(ae, 'noSmooth', 'mean');
for i = 1:length(densityByLatitude(1,:))
   samplingFreq = 86400 / mean(diff(crossingTimes(:,i)));
   densityByLatitude(:,i) = removePeriodicBackground(densityByLatitude(:,i), 0.3, samplingFreq, 0);
end
<<<<<<< HEAD

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
    end
    
    if catTimelag > 0
        timelag = [timelag mean(lagsAllIntervalLats)];
        errInLag = [errInLag (std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats)))];  
    else
        timelag(length(i:-step:minAllowedLatitude) + 1) = mean(lagsAllIntervalLats);
        errInLag(length(i:-step:minAllowedLatitude) + 1) = (std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats)));
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
    end
    
    timelag(length(i:-step:minAllowedLatitude)) = mean(lagsAllIntervalLats);
    errInLag(length(i:-step:minAllowedLatitude)) = std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats));
end

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
=======

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
    end
    
    if catTimelag > 0
        timelag = [timelag mean(lagsAllIntervalLats)];
        errInLag = [errInLag (std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats)))];  
    else
        timelag(length(i:-step:minAllowedLatitude) + 1) = mean(lagsAllIntervalLats);
        errInLag(length(i:-step:minAllowedLatitude) + 1) = (std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats)));
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
    end
    
    timelag(length(i:-step:minAllowedLatitude)) = mean(lagsAllIntervalLats);
    errInLag(length(i:-step:minAllowedLatitude)) = std(lagsAllIntervalLats) / sqrt(length(lagsAllIntervalLats));
end

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

end

>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb

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

<<<<<<< HEAD

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

=======
>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
function simpleColorMap(correctedDensity, msisDensity, latitude, timestamps10s, timeOfDay)
% simpleDataPlot(correctedDensity, latitude, timestamps10s)

figure;
timestampsInDays = timestamps10s / (60 * 60 * 24) + datenum('2009-11-01', 'yyyy-mm-dd');
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
<<<<<<< HEAD

end

function compareGoceDensityToMsis(goceDensity, msisDensity, ae, timestampsAeDatenum, timestampsDensityDatenum, results)
%

ratio = goceDensity ./ msisDensity;
ratioTrend = removePeriodicBackground(ratio, 125, 6, 0);

figure;
[hAx,hLine1,hLine2] = plotyy(timestampsDensityDatenum, ratio, timestampsAeDatenum, ae);
set(hLine1, 'LineStyle', 'none', 'Marker', '.')
% For some reason Matlab can't handle axis limits properly, so they need to
% be set manually.
ratioYmin = 0.1 * floor(min(ratio) / 0.1);
ratioYmax = 0.1 * ceil(max(ratio) / 0.1);
set(hAx(1), 'YLim', [ratioYmin ratioYmax], 'YTick', ratioYmin:0.1:ratioYmax);
set(hAx, 'XLim', [min(timestampsDensityDatenum) max(timestampsDensityDatenum)], 'XMinorTick', 'on'); 
aeYmax = 500 * ceil(max(ae) / 500);
set(hAx(2),'YLim', [0 aeYmax], 'YTick', 0:250:aeYmax);
datetick(hAx(1), 'x', 'yyyy-mm', 'keepticks', 'keeplimits')
datetick(hAx(2), 'x', 'yyyy-mm', 'keepticks', 'keeplimits')
rotateticklabel(hAx(1), 50);
rotateticklabel(hAx(2), 50);
hold on;
plot(timestampsDensityDatenum, ratioTrend, 'r-');
hold off;
%plotyy(timestampsInDays, correctedDensity, timestampsInDays, msisDensity270km);
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

plotCorrelation(msisDensity, goceDensity, 'NRLMSISE00 density', 'GOCE measured density', 1, results);

end

function results = plotAndAnalyzeChangesByOrbit(densityNoBg, magneticLatitude, averagedDensityNoBg, timestamps1minFixed, ...
    timestamps10s, timeOfDay, plotFigures, results)
% plotAndAnalyzeChangesByOrbit(densityNoBg, magneticLatitude, averagedDensityNoBg)
persistent densityByOrbitFigHandle
persistent residueFigHandle
persistent densityByOrbitAxesHandle
persistent residueAxisHandle

averagedDensityNoBg = averagedDensityNoBg(ismember(timestamps1minFixed, timestamps10s));
timestamps1minFixed = timestamps1minFixed(ismember(timestamps1minFixed, timestamps10s));
[~, peakBeginIndex, peakEndIndex] = limitToNearPeak(averagedDensityNoBg, 'noSmooth', 'mean');
peakBeginIndex = find(timestamps10s == timestamps1minFixed(peakBeginIndex));
peakEndIndex = find(timestamps10s == timestamps1minFixed(peakEndIndex));
orbits = splitIntoOrbits(magneticLatitude);

if plotFigures ~= 0
if ~isempty(strfind(lower(timeOfDay), 'morning')); densityByOrbitFigHandle = figure; subplotNum = 1; else subplotNum = 2; end
    figure(densityByOrbitFigHandle);
    densityByOrbitAxesHandle(subplotNum) = subplot(2,1,subplotNum);
    hold all;
end
linehandles = [];
relativeResidues = nan(size(magneticLatitude));
TADplot = nan(size(magneticLatitude));
calmOrbits = 5;
maxNumOfColorOrbits = 10;
extraColorMapOrbits = 5;
[beginOrbit, endOrbit] = findBeginAndEndOrbits(orbits, peakBeginIndex, peakEndIndex, calmOrbits);
if endOrbit - beginOrbit - 2 * calmOrbits > maxNumOfColorOrbits + extraColorMapOrbits
    endColorOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits - 1;
    endColorMapOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits + extraColorMapOrbits - 1;
elseif endOrbit - beginOrbit - 2 * calmOrbits > maxNumOfColorOrbits
    endColorOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits - 1;
    endColorMapOrbit = endOrbit - calmOrbits;
else
    endColorOrbit = endOrbit;
    endColorMapOrbit = endOrbit;
end

loopOrbits = beginOrbit:endOrbit;
TADPlotIndices = nan(size(magneticLatitude));
for i = 1:length(loopOrbits)
    indices = orbits(loopOrbits(i),1) : orbits(loopOrbits(i),2);
    smoothedDensity150s = smooth(densityNoBg(indices), 15);
    relativeResidues(indices) = (densityNoBg(indices) - smoothedDensity150s) ./ smoothedDensity150s;
    if (loopOrbits(i) <= endColorOrbit || loopOrbits(i) > endOrbit - calmOrbits) && plotFigures ~= 0
        h = plot(magneticLatitude(indices), smoothedDensity150s, 'LineWidth', 2);
        linehandles = [linehandles h];
    end
    
    if loopOrbits(i) <= endColorMapOrbit
        TADPlotIndices(indices) = indices;            
        smoothedDensity10400km = smooth(densityNoBg(indices), 133);
        smoothedDensity2600km = smooth(densityNoBg(indices), 33); 
        TADplot(indices) = (smoothedDensity2600km - smoothedDensity10400km) ./ smoothedDensity10400km;
    end    
end

%hold off;
%(calmOrbits + 1 : length(linehandles) - calmOrbits);
if plotFigures ~= 0
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

end
magneticLatitudeResiduePlot = magneticLatitude(~isnan(relativeResidues));
relativeResidues = relativeResidues(~isnan(relativeResidues));

if plotFigures ~= 0
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
end

stdSouth = std(relativeResidues(magneticLatitudeResiduePlot < 0));
stdNorth = std(relativeResidues(magneticLatitudeResiduePlot > 0));

[rowNum, ~] = size(results);
emptyCells = cellfun(@isempty,results);
[~, emptyColPositions] = find(emptyCells);
colNum = min(emptyColPositions);
if rowNum == 2
    colNum = length(results(rowNum,:)) + 1;
    results{1, colNum}     = ['Std SH ', timeOfDay];
    results{1, colNum + 1} = ['Std NH ', timeOfDay];
    results{1, colNum + 2} = ['Std SH/NH ', timeOfDay];
end
results{rowNum, colNum}     = stdSouth;
results{rowNum, colNum + 1} = stdNorth;
results{rowNum, colNum + 2} = stdSouth / stdNorth;

if plotFigures ~= 0
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
    xlim([min(confIntervalX) max(confIntervalX)]);
    hold all;
    plot(confIntervalX, confIntervalMean, 'k--');
    hold off;
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
    indices = TADPlotIndices(~isnan(TADPlotIndices));
    figure;
    secondsInDay = 24 * 60 *60;
    timestampsInDays = datenum('2009-11-01','yyyy-mm-dd') + timestamps10s / secondsInDay;
    scatter(timestampsInDays(indices), magneticLatitude(indices), 60, TADplot, '.');
    datetick('x', 'yyyy-mm-dd', 'keepticks', 'keeplimits')
    rotateticklabel(gca, 50);
    ylabel('IGRF Magnetic Latitude')
    title(['1300-5200km changes (TADs) [(2600 km smooth - 10400 km smooth) / 10400 km smooth] ', timeOfDay])
    colorbar
end
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
=======
>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb

if magneticLatitude(testIndex1) > magneticLatitude(testIndex2)
   goingSouth = 1; 
else
   goingSouth = 0;
end

<<<<<<< HEAD
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

=======
function compareGoceDensityToMsis(goceDensity, msisDensity, ae, timestampsAeDatenum, timestampsDensityDatenum, results)
%

ratio = goceDensity ./ msisDensity;
ratioTrend = removePeriodicBackground(ratio, 125, 6, 0);

figure;
[hAx,hLine1,hLine2] = plotyy(timestampsDensityDatenum, ratio, timestampsAeDatenum, ae);
set(hLine1, 'LineStyle', 'none', 'Marker', '.')
% For some reason Matlab can't handle axis limits properly, so they need to
% be set manually.
ratioYmin = 0.1 * floor(min(ratio) / 0.1);
ratioYmax = 0.1 * ceil(max(ratio) / 0.1);
set(hAx(1), 'YLim', [ratioYmin ratioYmax], 'YTick', ratioYmin:0.1:ratioYmax);
set(hAx, 'XLim', [min(timestampsDensityDatenum) max(timestampsDensityDatenum)], 'XMinorTick', 'on'); 
aeYmax = 500 * ceil(max(ae) / 500);
set(hAx(2),'YLim', [0 aeYmax], 'YTick', 0:250:aeYmax);
datetick(hAx(1), 'x', 'yyyy-mm', 'keepticks', 'keeplimits')
datetick(hAx(2), 'x', 'yyyy-mm', 'keepticks', 'keeplimits')
rotateticklabel(hAx(1), 50);
rotateticklabel(hAx(2), 50);
hold on;
plot(timestampsDensityDatenum, ratioTrend, 'r-');
hold off;
%plotyy(timestampsInDays, correctedDensity, timestampsInDays, msisDensity270km);
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

plotCorrelation(msisDensity, goceDensity, 'NRLMSISE00 density', 'GOCE measured density', 1, results);

end

function results = plotAndAnalyzeChangesByOrbit(densityNoBg, magneticLatitude, averagedDensityNoBg, timestamps1minFixed, ...
    timestamps10s, timeOfDay, plotFigures, results)
% plotAndAnalyzeChangesByOrbit(densityNoBg, magneticLatitude, averagedDensityNoBg)
persistent densityByOrbitFigHandle
persistent residueFigHandle
persistent densityByOrbitAxesHandle
persistent residueAxisHandle

averagedDensityNoBg = averagedDensityNoBg(ismember(timestamps1minFixed, timestamps10s));
timestamps1minFixed = timestamps1minFixed(ismember(timestamps1minFixed, timestamps10s));
[~, peakBeginIndex, peakEndIndex] = limitToNearPeak(averagedDensityNoBg, 'noSmooth', 'mean');
peakBeginIndex = find(timestamps10s == timestamps1minFixed(peakBeginIndex));
peakEndIndex = find(timestamps10s == timestamps1minFixed(peakEndIndex));
% orbits = splitIntoOrbits(magneticLatitude);

[latBeginIndex, latEndIndex] = limitLatitudeToIntegerMultipleOfOrbitalPeriod(magneticLatitude);
limitedLatitude = magneticLatitude(latBeginIndex:latEndIndex);
limitedTimestamps = timestamps10s(latBeginIndex:latEndIndex);
[~, exactOrbitIndices] = splitIntoOrbits(limitedLatitude, limitedTimestamps);
limitedLatitude = limitedLatitude(exactOrbitIndices);
limitedTimestamps = limitedTimestamps(exactOrbitIndices);
[orbits, ~] = splitIntoOrbits(limitedLatitude, limitedTimestamps);

[minAllowedLatitude, maxAllowedLatitude] = findInterpolationLimits(limitedLatitude);
orbitsToDelete = find(abs(limitedLatitude(orbits(:,2)) - limitedLatitude(orbits(:,1))) < ...
    (maxAllowedLatitude - minAllowedLatitude));
min(abs(limitedLatitude(orbits(:,2)) - limitedLatitude(orbits(:,1))))
newIndices = 1:length(limitedLatitude);
for i = 1:length(orbitsToDelete)
    newIndices = setdiff(newIndices, (orbits(orbitsToDelete(i), 1) : orbits(orbitsToDelete(i), 2)));
>>>>>>> 70ac792af66c5f27a724151c372fd920b36242cb
end
limitedTimestamps = limitedTimestamps(newIndices);
limitedLatitude = limitedLatitude(newIndices);

if plotFigures ~= 0
if ~isempty(strfind(lower(timeOfDay), 'morning')); densityByOrbitFigHandle = figure; subplotNum = 1; else subplotNum = 2; end
    figure(densityByOrbitFigHandle);
    densityByOrbitAxesHandle(subplotNum) = subplot(2,1,subplotNum);
    hold all;
end
linehandles = [];
relativeResidues = nan(size(magneticLatitude));
TADplot = nan(size(magneticLatitude));
calmOrbits = 5;
maxNumOfColorOrbits = 10;
extraColorMapOrbits = 5;
[beginOrbit, endOrbit] = findBeginAndEndOrbits(orbits, peakBeginIndex, peakEndIndex, calmOrbits);
if endOrbit - beginOrbit - 2 * calmOrbits > maxNumOfColorOrbits + extraColorMapOrbits
    endColorOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits - 1;
    endColorMapOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits + extraColorMapOrbits - 1;
elseif endOrbit - beginOrbit - 2 * calmOrbits > maxNumOfColorOrbits
    endColorOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits - 1;
    endColorMapOrbit = endOrbit - calmOrbits;
else
    endColorOrbit = endOrbit;
    endColorMapOrbit = endOrbit;
end

loopOrbits = beginOrbit:endOrbit;
TADPlotIndices = nan(size(magneticLatitude));
for i = 1:length(loopOrbits)
    indices = orbits(loopOrbits(i),1) : orbits(loopOrbits(i),2);
    smoothedDensity150s = smooth(densityNoBg(indices), 15);
    relativeResidues(indices) = (densityNoBg(indices) - smoothedDensity150s) ./ smoothedDensity150s;
    if (loopOrbits(i) <= endColorOrbit || loopOrbits(i) > endOrbit - calmOrbits) && plotFigures ~= 0
        h = plot(magneticLatitude(indices), smoothedDensity150s, 'LineWidth', 2);
        linehandles = [linehandles h];
    end
    
    if loopOrbits(i) <= endColorMapOrbit
        TADPlotIndices(indices) = indices;            
        smoothedDensity10400km = smooth(densityNoBg(indices), 133);
        smoothedDensity2600km = smooth(densityNoBg(indices), 33); 
        TADplot(indices) = (smoothedDensity2600km - smoothedDensity10400km) ./ smoothedDensity10400km;
    end    
end

%hold off;
%(calmOrbits + 1 : length(linehandles) - calmOrbits);
if plotFigures ~= 0
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

function th=rotateticklabel(h,rot,demo)
%ROTATETICKLABEL rotates tick labels
%   TH=ROTATETICKLABEL(H,ROT) is the calling form where H is a handle to
%   the axis that contains the XTickLabels that are to be rotated. ROT is
%   an optional parameter that specifies the angle of rotation. The default
%   angle is 90. TH is a handle to the text objects created. For long
%   strings such as those produced by datetick, you may have to adjust the
%   position of the axes so the labels don't get cut off.
%
%   Of course, GCA can be substituted for H if desired.
%
%   TH=ROTATETICKLABEL([],[],'demo') shows a demo figure.
%
%   Known deficiencies: if tick labels are raised to a power, the power
%   will be lost after rotation.
%
%   See also datetick.

% Copyright (c) 2005, Andrew Bliss
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

%   Written Oct 14, 2005 by Andy Bliss
%   Copyright 2005 by Andy Bliss

%DEMO:
if nargin==3
    x=[now-.7 now-.3 now];
    y=[20 35 15];
    figure
    plot(x,y,'.-')
    datetick('x',0,'keepticks')
    h=gca;
    set(h,'position',[0.13 0.35 0.775 0.55])
    rot=90;
end

%set the default rotation if user doesn't specify
if nargin==1
    rot=90;
end
%make sure the rotation is in the range 0:360 (brute force method)
while rot>360
    rot=rot-360;
end
while rot<0
    rot=rot+360;
end
%get current tick labels
a=get(h,'XTickLabel');
%erase current tick labels from figure
set(h,'XTickLabel',[]);
%get tick label positions
b=get(h,'XTick');
c=get(h,'YTick');
%make new tick labels
if rot<180
    th=text(b,repmat(c(1)-.1*(c(2)-c(1)),length(b),1),a,'HorizontalAlignment','right','rotation',rot);
else
    th=text(b,repmat(c(1)-.1*(c(2)-c(1)),length(b),1),a,'HorizontalAlignment','left','rotation',rot);
end

end

function cell2csv(fileName, cellArray, separator, excelYear, decimal)
% Writes cell array content into a *.csv file.
% 
% CELL2CSV(fileName, cellArray, separator, excelYear, decimal)
%
% fileName     = Name of the file to save. [ i.e. 'text.csv' ]
% cellArray    = Name of the Cell Array where the data is in
% separator    = sign separating the values (default = ';')
% excelYear    = depending on the Excel version, the cells are put into
%                quotes before they are written to the file. The separator
%                is set to semicolon (;)
% decimal      = defines the decimal separator (default = '.')
%
%         by Sylvain Fiedler, KA, 2004
% updated by Sylvain Fiedler, Metz, 06
% fixed the logical-bug, Kaiserslautern, 06/2008, S.Fiedler
% added the choice of decimal separator, 11/2010, S.Fiedler
% 
% Copyright (c) 2004-2010, Sylvain Fiedler
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

% Checking für optional Variables
if ~exist('separator', 'var')
    separator = ',';
end

if ~exist('excelYear', 'var')
    excelYear = 1997;
end

if ~exist('decimal', 'var')
    decimal = '.';
end

% Setting separator for newer excelYears
if excelYear > 2000
    separator = ';';
end

% Write file
datei = fopen(fileName, 'w');

for z=1:size(cellArray, 1)
    for s=1:size(cellArray, 2)
        
        var = eval(['cellArray{z,s}']);
        % If zero, then empty cell
        if size(var, 1) == 0
            var = '';
        end
        % If numeric -> String
        if isnumeric(var)
            var = num2str(var);
            % Conversion of decimal separator (4 Europe & South America)
            % http://commons.wikimedia.org/wiki/File:DecimalSeparator.svg
            if decimal ~= '.'
                var = strrep(var, '.', decimal);
            end
        end
        % If logical -> 'true' or 'false'
        if islogical(var)
            if var == 1
                var = 'TRUE';
            else
                var = 'FALSE';
            end
        end
        % If newer version of Excel -> Quotes 4 Strings
        if excelYear > 2000
            var = ['"' var '"'];
        end
        
        % OUTPUT value
        fprintf(datei, '%s', var);
        
        % OUTPUT separator
        if s ~= size(cellArray, 2)
            fprintf(datei, separator);
        end
    end
    if z ~= size(cellArray, 1) % prevent a empty line at EOF
        % OUTPUT newline
        fprintf(datei, '\n');
    end
end
% Closing file
fclose(datei);
% END
end
magneticLatitudeResiduePlot = magneticLatitude(~isnan(relativeResidues));
relativeResidues = relativeResidues(~isnan(relativeResidues));

if plotFigures ~= 0
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
end

stdSouth = std(relativeResidues(magneticLatitudeResiduePlot < 0));
stdNorth = std(relativeResidues(magneticLatitudeResiduePlot > 0));

[rowNum, ~] = size(results);
emptyCells = cellfun(@isempty,results);
[~, emptyColPositions] = find(emptyCells);
colNum = min(emptyColPositions);
if rowNum == 2
    colNum = length(results(rowNum,:)) + 1;
    results{1, colNum}     = ['Std SH ', timeOfDay];
    results{1, colNum + 1} = ['Std NH ', timeOfDay];
    results{1, colNum + 2} = ['Std SH/NH ', timeOfDay];
end
results{rowNum, colNum}     = stdSouth;
results{rowNum, colNum + 1} = stdNorth;
results{rowNum, colNum + 2} = stdSouth / stdNorth;

if plotFigures ~= 0
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
    xlim([min(confIntervalX) max(confIntervalX)]);
    hold all;
    plot(confIntervalX, confIntervalMean, 'k--');
    hold off;
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
    indices = TADPlotIndices(~isnan(TADPlotIndices));
    TADlatitude = magneticLatitude(indices);
    secondsInDay = 24 * 60 *60;
    timestampsInDays = datenum('2009-11-01','yyyy-mm-dd') + limitedTimestamps / secondsInDay; % EI OIKEIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TADtime = timestampsInDays(indices);

    oneQuarterDegreeStep = minAllowedLatitude:0.25:maxAllowedLatitude;

    figure;
    scatter(TADtime, TADlatitude, 60, TADplot,'.');
    
    regriddedTime(:,1) = latitudeCrossingTimes(TADlatitude, TADtime, oneQuarterDegreeStep(i));
    for i = 2:length(oneQuarterDegreeStep)
        times = latitudeCrossingTimes(TADlatitude, TADtime, oneQuarterDegreeStep(i));
        if length(times) == length(regriddedTime(:,1))
            regriddedTime(:,i) = times; 
        else
            length(times)
        end
    end

    % regriddedDensity = griddata(timestamps10s, magneticLatitude, correctedDensity, regriddedTime, latitudeMatrix, 'natural');
    regriddedDensity = interp1(TADtime, TADplot, regriddedTime, 'spline');

    numOfOrbits = length(regriddedTime(:,1));
    numOfValuesInOrbit = length(regriddedTime(1,:));
    for i = 1:numOfValuesInOrbit
        timeThisLatitude = regriddedTime(:,i);
        densityThisLatitude = regriddedDensity(:,i);

        tInterp = interp1(1:numOfOrbits, timeThisLatitude, 1:1/10:numOfOrbits);
        interpolatedDensity = interp1(timeThisLatitude, densityThisLatitude, tInterp, 'spline');
        thisLatitude = ones(length(tInterp), 1) * oneQuarterDegreeStep(i);
        latitudeMatrix(:,i) = thisLatitude;
        densityMatrix(:,i) = interpolatedDensity;
        timeMatrix(:,i) = tInterp;
    end
    
    figure;
    %surf(timeMatrix, latitudeMatrix, densityMatrix, 'EdgeColor', 'None')
    latitude = reshape(latitudeMatrix, [], 1);
    density = reshape(densityMatrix, [], 1);
    time = reshape(timeMatrix, [], 1);
    density(density > 0.5 | density < -0.5) = nan(1);
    scatter(time, latitude, 60, density,'.');
    view(2);
    ylabel('IGRF Magnetic Latitude')
    title(['1300-5200km changes (TADs) [(2600 km smooth - 10400 km smooth) / 10400 km smooth] ', timeOfDay])
    colorbar
%     colormap jet(500)
    ylim([minAllowedLatitude maxAllowedLatitude]);
    xlim([min(TADtime) max(TADtime)]);

end
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

function th=rotateticklabel(h,rot,demo)
%ROTATETICKLABEL rotates tick labels
%   TH=ROTATETICKLABEL(H,ROT) is the calling form where H is a handle to
%   the axis that contains the XTickLabels that are to be rotated. ROT is
%   an optional parameter that specifies the angle of rotation. The default
%   angle is 90. TH is a handle to the text objects created. For long
%   strings such as those produced by datetick, you may have to adjust the
%   position of the axes so the labels don't get cut off.
%
%   Of course, GCA can be substituted for H if desired.
%
%   TH=ROTATETICKLABEL([],[],'demo') shows a demo figure.
%
%   Known deficiencies: if tick labels are raised to a power, the power
%   will be lost after rotation.
%
%   See also datetick.

% Copyright (c) 2005, Andrew Bliss
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

%   Written Oct 14, 2005 by Andy Bliss
%   Copyright 2005 by Andy Bliss

%DEMO:
if nargin==3
    x=[now-.7 now-.3 now];
    y=[20 35 15];
    figure
    plot(x,y,'.-')
    datetick('x',0,'keepticks')
    h=gca;
    set(h,'position',[0.13 0.35 0.775 0.55])
    rot=90;
end

%set the default rotation if user doesn't specify
if nargin==1
    rot=90;
end
%make sure the rotation is in the range 0:360 (brute force method)
while rot>360
    rot=rot-360;
end
while rot<0
    rot=rot+360;
end
%get current tick labels
a=get(h,'XTickLabel');
%erase current tick labels from figure
set(h,'XTickLabel',[]);
%get tick label positions
b=get(h,'XTick');
c=get(h,'YTick');
%make new tick labels
if rot<180
    th=text(b,repmat(c(1)-.1*(c(2)-c(1)),length(b),1),a,'HorizontalAlignment','right','rotation',rot);
else
    th=text(b,repmat(c(1)-.1*(c(2)-c(1)),length(b),1),a,'HorizontalAlignment','left','rotation',rot);
end

end

function cell2csv(fileName, cellArray, separator, excelYear, decimal)
% Writes cell array content into a *.csv file.
% 
% CELL2CSV(fileName, cellArray, separator, excelYear, decimal)
%
% fileName     = Name of the file to save. [ i.e. 'text.csv' ]
% cellArray    = Name of the Cell Array where the data is in
% separator    = sign separating the values (default = ';')
% excelYear    = depending on the Excel version, the cells are put into
%                quotes before they are written to the file. The separator
%                is set to semicolon (;)
% decimal      = defines the decimal separator (default = '.')
%
%         by Sylvain Fiedler, KA, 2004
% updated by Sylvain Fiedler, Metz, 06
% fixed the logical-bug, Kaiserslautern, 06/2008, S.Fiedler
% added the choice of decimal separator, 11/2010, S.Fiedler
% 
% Copyright (c) 2004-2010, Sylvain Fiedler
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

% Checking für optional Variables
if ~exist('separator', 'var')
    separator = ',';
end

if ~exist('excelYear', 'var')
    excelYear = 1997;
end

if ~exist('decimal', 'var')
    decimal = '.';
end

% Setting separator for newer excelYears
if excelYear > 2000
    separator = ';';
end

% Write file
datei = fopen(fileName, 'w');

for z=1:size(cellArray, 1)
    for s=1:size(cellArray, 2)
        
        var = eval(['cellArray{z,s}']);
        % If zero, then empty cell
        if size(var, 1) == 0
            var = '';
        end
        % If numeric -> String
        if isnumeric(var)
            var = num2str(var);
            % Conversion of decimal separator (4 Europe & South America)
            % http://commons.wikimedia.org/wiki/File:DecimalSeparator.svg
            if decimal ~= '.'
                var = strrep(var, '.', decimal);
            end
        end
        % If logical -> 'true' or 'false'
        if islogical(var)
            if var == 1
                var = 'TRUE';
            else
                var = 'FALSE';
            end
        end
        % If newer version of Excel -> Quotes 4 Strings
        if excelYear > 2000
            var = ['"' var '"'];
        end
        
        % OUTPUT value
        fprintf(datei, '%s', var);
        
        % OUTPUT separator
        if s ~= size(cellArray, 2)
            fprintf(datei, separator);
        end
    end
    if z ~= size(cellArray, 1) % prevent a empty line at EOF
        % OUTPUT newline
        fprintf(datei, '\n');
    end
end
% Closing file
fclose(datei);
% END
end
