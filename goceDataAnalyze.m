function goceDataAnalyze( varargin )
% goceDataAnalyze( threshold )
%   Detailed explanation goes here

results = initialize();

[threshold, plotDates] = processInputArguments(varargin, nargin);

[ae, ap, absB, akasofuEpsilon, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude,...
 eveningMagneticLatitude, cellArrayLength, firstDatenum] ...
 = compareToMsisAndGiveVariables(threshold, results);

plotFigures = 0;
for i = 1:cellArrayLength
    if userRequestedThisStormToPlot(plotDates, timestampsDatenum{i})
        plotFigures = 1;
    end
    results = writeTimeIntervalAndMaxAeToResultArray(timestampsDatenum{i}, ae{i}, results);
    progress(i, cellArrayLength);
    
%     if plotFigures ~= 0
%         timeseriesFigHandle = plotTimeseries(firstDatenum, timestamps1min{i}, timestamps1minFixed{i}, timestampsAbsB{i},...
%             timestamps3h{i}, timestamps3hFixed{i}, ae{i}, ap{i}, absB{i},averagedDensityNoBg{i}, density3h{i});
%     else
%         timeseriesFigHandle = nan(1);
%     end
% 
%     [results, aeIntegral, timestampsAeInt] = plotAndCalculateCorrelation(firstDatenum, timestamps1min{i}, timestamps1minFixed{i}, ...
%         ae{i}, averagedDensityNoBg{i}, 'AE', plotFigures, results, timeseriesFigHandle); 
%     results = plotAndCalculateCorrelation(firstDatenum, timestamps3h{i}, timestamps3hFixed{i}, ap{i}, density3h{i}, 'ap',...
%         plotFigures, results, timeseriesFigHandle); 
%     results = plotAndCalculateCorrelation(firstDatenum, timestampsAbsB{i}, timestamps1minFixed{i}, absB{i}, averagedDensityNoBg{i},...
%         'IMF |B|', plotFigures, results, timeseriesFigHandle); 
%     results = plotAndCalculateCorrelation(firstDatenum, timestampsEpsilon{i}, timestamps1minFixed{i}, akasofuEpsilon{i}, ...
%         averagedDensityNoBg{i}, 'Akasofu Epsilon', plotFigures, results, timeseriesFigHandle);
% %     
%     results = plotAndAnalyzeDensityByLatitude(firstDatenum, ae{i}, timestamps1min{i}, aeIntegral, timestampsAeInt, timestamps1minFixed{i}, ...
%         morningDensityNoBg{i}, morningMsisDensity{i}, morningTimestamps10s{i}, morningMagneticLatitude{i}, 'Morning', plotFigures, results);
%     results = plotAndAnalyzeDensityByLatitude(firstDatenum, ae{i}, timestamps1min{i}, aeIntegral, timestampsAeInt, timestamps1minFixed{i}, ...
%         eveningDensityNoBg{i}, eveningMsisDensity{i}, eveningTimestamps10s{i}, eveningMagneticLatitude{i}, 'Evening', plotFigures, results);
    
    results = plotAndAnalyzeChangesByOrbit(firstDatenum, morningDensityNoBg{i}, morningMagneticLatitude{i}, averagedDensityNoBg{i},...
        timestamps1minFixed{i}, morningTimestamps10s{i}, 'Morning', plotFigures, results);
    results = plotAndAnalyzeChangesByOrbit(firstDatenum, eveningDensityNoBg{i}, eveningMagneticLatitude{i}, averagedDensityNoBg{i},...
        timestamps1minFixed{i}, eveningTimestamps10s{i}, 'Evening', plotFigures, results);
    
    plotFigures = 0;
end

makeSummaryOfResults(results);

end

function results = initialize()
% initialize()

format compact
if matlabpool('size') <= 0
    matlabpool open
end
results = {};

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

for i = 1:length(plotDates)
    plotOrNot = plotDates(i) >= min(timestampsDatenum) && plotDates(i) <= max(timestampsDatenum);
    if plotOrNot == 1
        break;
    end
end

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

%results = computeResultMeanAndStd(results);

cell2csv('goceResults.csv', results);

end


function results = computeResultMeanAndStd(results)
%

numericVariables = cellfun(@double, results(2:end, 3:end));
numericVariablesMean = mean(numericVariables);
numericVariablesStd = std(numericVariables);
[rows, ~] = size(results);
results(rows + 1, 3:end) = num2cell(numericVariablesMean);
results(rows + 2, 3:end) = num2cell(numericVariablesStd);
results(rows + 1, 1:2) = {'Mean', 'NaN'};
results(rows + 2, 1:2) = {'Std', 'NaN'};

end