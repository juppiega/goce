function goceDataAnalyze( varargin )
% goceDataAnalyze( threshold )
%   Detailed explanation goes here

results = initialize();

[threshold, plotDates] = processInputArguments(varargin, nargin);

[ae, ap, absB, vBz, akasofuEpsilon, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningAeProxy, eveningAeProxy, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude,...
 eveningMagneticLatitude, cellArrayLength, firstDatenum] ...
 = compareToMsisAndGiveVariables(threshold, results);

plotFigures = 0;
for i = 1:cellArrayLength
    if userRequestedThisStormToPlot(plotDates, timestampsDatenum{i})
        plotFigures = 1;
    end
    results = writeTimeIntervalAndMaxAeToResultArray(timestampsDatenum{i}, ae{i}, results);
    progressBar = progress(i, cellArrayLength);
    
    if plotFigures ~= 0
        timeseriesFigHandle = plotTimeseries(firstDatenum, timestamps1min{i}, timestamps1minFixed{i}, timestampsAbsB{i},...
            timestamps3h{i}, timestamps3hFixed{i}, ae{i}, ap{i}, absB{i},averagedDensityNoBg{i}, density3h{i});
    else
        timeseriesFigHandle = nan(1);
    end

    [results, aeIntegral, timestampsAeInt] = plotAndCalculateCorrelation(firstDatenum, timestamps1min{i}, timestamps1minFixed{i}, ...
        ae{i}, averagedDensityNoBg{i}, 'AE', plotFigures, results, timeseriesFigHandle); 
%     results = plotAndCalculateCorrelation(firstDatenum, timestamps3h{i}, timestamps3hFixed{i}, ap{i}, density3h{i}, 'ap',...
%         plotFigures, results, timeseriesFigHandle); 
%     results = plotAndCalculateCorrelation(firstDatenum, timestampsAbsB{i}, timestamps1minFixed{i}, absB{i}, averagedDensityNoBg{i},...
%         'IMF |B|', plotFigures, results, timeseriesFigHandle); 
%     results = plotAndCalculateCorrelation(firstDatenum, timestampsEpsilon{i}, timestamps1minFixed{i}, akasofuEpsilon{i}, ...
%         averagedDensityNoBg{i}, 'Akasofu Epsilon', plotFigures, results, timeseriesFigHandle);
%     results = plotAndCalculateCorrelation(firstDatenum, timestampsEpsilon{i}, timestamps1minFixed{i}, vBz{i}, ...
%         averagedDensityNoBg{i}, '|V| * Bz', plotFigures, results, timeseriesFigHandle);
%     
    results = plotAndAnalyzeDensityByLatitude(firstDatenum, ae{i}, timestamps1min{i}, aeIntegral, timestampsAeInt, timestamps1minFixed{i}, ...
        morningDensityNoBg{i}, morningMsisDensity{i}, morningAeProxy{i}, morningTimestamps10s{i}, morningMagneticLatitude{i}, 'Morning', plotFigures, results);
    results = plotAndAnalyzeDensityByLatitude(firstDatenum, ae{i}, timestamps1min{i}, aeIntegral, timestampsAeInt, timestamps1minFixed{i}, ...
        eveningDensityNoBg{i}, eveningMsisDensity{i}, eveningAeProxy{i}, eveningTimestamps10s{i}, eveningMagneticLatitude{i}, 'Evening', plotFigures, results);
    
%     results = plotAndAnalyzeChangesByOrbit(firstDatenum, morningDensityNoBg{i}, morningMagneticLatitude{i}, averagedDensityNoBg{i},...
%         timestamps1minFixed{i}, morningTimestamps10s{i}, 'Morning', plotFigures, results);
%     results = plotAndAnalyzeChangesByOrbit(firstDatenum, eveningDensityNoBg{i}, eveningMagneticLatitude{i}, averagedDensityNoBg{i},...
%         timestamps1minFixed{i}, eveningTimestamps10s{i}, 'Evening', plotFigures, results);
    
    plotFigures = 0;
end

progressBar.stop;
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

function p = progress(i, cellArrayLength)
%

persistent progressBar;

if i == 1
    barWidth = 50;
    progressBar = TimedProgressBar( cellArrayLength, barWidth, ...
                    ['Analyzing ', num2str(cellArrayLength) ,' storms, ETA '], ...
                    '. Now at ', ...
                    'Concluded in ' );
end

progressBar.progress;
p = progressBar;

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