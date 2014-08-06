function goceDataAnalyze( varargin )
% goceDataAnalyze( threshold )
%   Detailed explanation goes here

initialize()
results = {};

[threshold, plotDates] = processInputArguments(varargin, nargin);

[ae, ap, absB, akasofuEpsilon, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude,...
 eveningMagneticLatitude, cellArrayLength, firstDatenum] ...
 = variables(threshold, results);

plotFigures = 0;
for i = 1:cellArrayLength
    if userRequestedThisStormToPlot(plotDates, timestampsDatenum{i})
        plotFigures = 1;
    end
    results = writeTimeIntervalAndMaxAeToResultArray(timestampsDatenum{i}, ae{i}, results);
    progress(i, cellArrayLength);
    
    if plotFigures ~= 0
        timeseriesFigHandle = plotTimeseries(firstDatenum, timestamps1min{i}, timestamps1minFixed{i}, timestampsAbsB{i},...
            timestamps3h{i}, timestamps3hFixed{i}, ae{i}, ap{i}, absB{i},averagedDensityNoBg{i}, density3h{i});
    else
        timeseriesFigHandle = nan(1);
    end

    results = plotAndCalculateCorrelation(firstDatenum, timestamps1min{i}, timestamps1minFixed{i}, ae{i}, averagedDensityNoBg{i}, 'AE', plotFigures, results, timeseriesFigHandle); 
    results = plotAndCalculateCorrelation(firstDatenum, timestamps3h{i}, timestamps3hFixed{i}, ap{i}, density3h{i}, 'ap', plotFigures, results, timeseriesFigHandle); 
    results = plotAndCalculateCorrelation(firstDatenum, timestampsAbsB{i}, timestamps1minFixed{i}, absB{i}, averagedDensityNoBg{i}, 'IMF |B|', plotFigures, results, timeseriesFigHandle); 
    results = plotAndCalculateCorrelation(firstDatenum, timestampsEpsilon{i}, timestamps1minFixed{i}, akasofuEpsilon{i}, averagedDensityNoBg{i}, 'Akasofu Epsilon', plotFigures, results, timeseriesFigHandle);
    
    results = plotAndAnalyzeDensityByLatitude(firstDatenum, ae{i}, timestamps1min{i}, timestamps1minFixed{i}, ...
        morningDensityNoBg{i}, morningMsisDensity{i}, morningTimestamps10s{i}, morningMagneticLatitude{i}, 'Morning', plotFigures, results);
    results = plotAndAnalyzeDensityByLatitude(firstDatenum, ae{i}, timestamps1min{i}, timestamps1minFixed{i}, ...
        eveningDensityNoBg{i}, eveningMsisDensity{i}, eveningTimestamps10s{i}, eveningMagneticLatitude{i}, 'Evening', plotFigures, results);
    
    results = plotAndAnalyzeChangesByOrbit(firstDatenum, morningDensityNoBg{i}, morningMagneticLatitude{i}, averagedDensityNoBg{i},...
        timestamps1minFixed{i}, morningTimestamps10s{i}, 'Morning', plotFigures, results);
    results = plotAndAnalyzeChangesByOrbit(firstDatenum, eveningDensityNoBg{i}, eveningMagneticLatitude{i}, averagedDensityNoBg{i},...
        timestamps1minFixed{i}, eveningTimestamps10s{i}, 'Evening', plotFigures, results);
    
    plotFigures = 0;
end

makeSummaryOfResults(results);

end

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

plotOrNot = plotDates >= min(timestampsDatenum) & plotDates <= max(timestampsDatenum);

end

