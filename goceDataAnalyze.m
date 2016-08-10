function goceDataAnalyze( varargin )
% goceDataAnalyze( threshold )
%   Detailed explanation goes here%

results = initialize();

[threshold, plotDates] = processInputArguments(varargin, nargin);

[ae, ap, absB, vBz, akasofuEpsilon, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, morningJbDensity, morningDtmDensity, eveningDtmDensity, eveningMsisDensity, eveningJbDensity, morningTiegcmDensity, eveningTiegcmDensity, morningAeProxy, eveningAeProxy, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum,  morningMagneticLatitude,...
 eveningMagneticLatitude, latitude, timestamps10sFixed, cellArrayLength, firstDatenum, winds] ...
 = compareToModelsAndGiveVariables(threshold, results);

plotFigures = 0;
for i = 1:cellArrayLength
    if userRequestedThisStormToPlot(plotDates, timestampsDatenum{i})
        plotFigures = 1;
    end
    results = writeTimeIntervalAndMaxAeToResultArray(timestampsDatenum{i}, ae{i}, results);
    progressBar = progress(i, cellArrayLength);
    
    if plotFigures ~= 0
        timeseriesFigHandle = plotTimeseries(firstDatenum, timestamps1min{i}, timestamps1minFixed{i}, timestampsAbsB{i},...
            timestamps3h{i}, timestamps3hFixed{i}, ae{i}, ap{i}, absB{i},averagedDensityNoBg{i}, density3h{i}, morningAeProxy{i}, eveningAeProxy{i}, ...
            morningMsisDensity{i}, morningJbDensity{i}, eveningMsisDensity{i}, eveningJbDensity{i}, morningTimestamps10s{i}, eveningTimestamps10s{i});
    else
        timeseriesFigHandle = nan(1);
    end

    [results, aeIntegral, timestampsAeInt] = plotAndCalculateCorrelation(firstDatenum, timestamps1min{i}, morningTimestamps10s{i}, eveningTimestamps10s{i}, ...
        ae{i}, morningDensityNoBg{i}, eveningDensityNoBg{i}, latitude, timestamps10sFixed, 'AE', plotFigures, results, timeseriesFigHandle); 
%     results = plotAndCalculateCorrelation(firstDatenum, timestamps3h{i}, timestamps3hFixed{i}, timestamps3hFixed{i}, ap{i}, density3h{i}, density3h{i}, latitude, timestamps10sFixed, 'ap',...
%         plotFigures, results, timeseriesFigHandle); 
%     results = plotAndCalculateCorrelation(firstDatenum, timestampsAbsB{i}, morningTimestamps10s{i}, eveningTimestamps10s{i}, absB{i}, morningDensityNoBg{i}, eveningDensityNoBg{i},...
%         latitude, timestamps10sFixed, 'IMF |B|', plotFigures, results, timeseriesFigHandle); 
%     results = plotAndCalculateCorrelation(firstDatenum, timestampsEpsilon{i}, morningTimestamps10s{i}, eveningTimestamps10s{i}, akasofuEpsilon{i}, ...
%          morningDensityNoBg{i}, eveningDensityNoBg{i}, latitude, timestamps10sFixed, 'Akasofu Epsilon', plotFigures, results, timeseriesFigHandle);
%     results = plotAndCalculateCorrelation(firstDatenum, timestampsEpsilon{i}, morningTimestamps10s{i}, eveningTimestamps10s{i}, vBz{i}, ...
%         morningDensityNoBg{i}, eveningDensityNoBg{i}, latitude, timestamps10sFixed, '|V| * Bz', plotFigures, results, timeseriesFigHandle);
%     
% 
%     results = writeCorrelationsToResults(morningMsisDensity{i}, eveningMsisDensity{i}, morningTimestamps10s{i}, eveningTimestamps10s{i},...
%         morningDensityNoBg{i}, eveningDensityNoBg{i}, latitude, timestamps10sFixed, timestamps1min{i}, ae{i}, results, 'MSIS');
%     results = writeCorrelationsToResults(morningJbDensity{i}, eveningJbDensity{i}, morningTimestamps10s{i}, eveningTimestamps10s{i},...
%         morningDensityNoBg{i}, eveningDensityNoBg{i}, latitude, timestamps10sFixed, timestamps1min{i}, ae{i}, results, 'JB');
%     results = writeCorrelationsToResults(morningAeProxy{i}, eveningAeProxy{i}, morningTimestamps10s{i}, eveningTimestamps10s{i},...
%         morningDensityNoBg{i}, eveningDensityNoBg{i}, latitude, timestamps10sFixed, timestamps1min{i}, ae{i}, results, 'AE Int.');
%     results = writeCorrelationsToResults(morningDtmDensity{i}, eveningDtmDensity{i}, morningTimestamps10s{i}, eveningTimestamps10s{i},...
%         morningDensityNoBg{i}, eveningDensityNoBg{i}, latitude, timestamps10sFixed, timestamps1min{i}, ae{i}, results, 'DTM');
%     results = writeCorrelationsToResults(morningTiegcmDensity{i}, eveningTiegcmDensity{i}, morningTimestamps10s{i}, eveningTimestamps10s{i},...
%         morningDensityNoBg{i}, eveningDensityNoBg{i}, latitude, timestamps10sFixed, timestamps1min{i}, ae{i}, results, 'TIEGCM');
    
    results = plotAndAnalyzeDensityByLatitude(firstDatenum, ae{i}, timestamps1min{i}, aeIntegral, timestampsAeInt, timestamps1minFixed{i}, ...
        morningDensityNoBg{i}, morningMsisDensity{i}, morningJbDensity{i}, morningDtmDensity{i}, morningTiegcmDensity{i}, morningAeProxy{i}, morningTimestamps10s{i}, morningMagneticLatitude{i}, i, winds, 'Morning', plotFigures, results);
    results = plotAndAnalyzeDensityByLatitude(firstDatenum, ae{i}, timestamps1min{i}, aeIntegral, timestampsAeInt, timestamps1minFixed{i}, ...
        eveningDensityNoBg{i}, eveningMsisDensity{i}, eveningJbDensity{i}, eveningDtmDensity{i}, eveningTiegcmDensity{i}, eveningAeProxy{i}, eveningTimestamps10s{i}, eveningMagneticLatitude{i}, i, winds, 'Evening', plotFigures, results);
    
%     results = plotAndAnalyzeChangesByOrbit(firstDatenum, morningDensityNoBg{i}, morningMagneticLatitude{i}, averagedDensityNoBg{i},...
%         timestamps1minFixed{i}, morningTimestamps10s{i}, 'Morning', 'GOCE', plotFigures, results);
%     results = plotAndAnalyzeChangesByOrbit(firstDatenum, eveningDensityNoBg{i}, eveningMagneticLatitude{i}, averagedDensityNoBg{i},...
%         timestamps1minFixed{i}, eveningTimestamps10s{i}, 'Evening', 'GOCE', plotFigures, results);
%     
%     results = plotAndAnalyzeChangesByOrbit(firstDatenum, morningTiegcmDensity{i}, morningMagneticLatitude{i}, averagedDensityNoBg{i},...
%         timestamps1minFixed{i}, morningTimestamps10s{i}, 'Morning', 'TIEGCM', plotFigures, results);
%     results = plotAndAnalyzeChangesByOrbit(firstDatenum, eveningTiegcmDensity{i}, eveningMagneticLatitude{i}, averagedDensityNoBg{i},...
%         timestamps1minFixed{i}, eveningTimestamps10s{i}, 'Evening', 'TIEGCM', plotFigures, results);
    
    plotFigures = 0;
end

progressBar.stop;
makeSummaryOfResults(results);

end

function results = initialize()
% initialize()

% poolobj = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(poolobj)
%     parpool(16);
% end

if(matlabpool('size')==0)
    matlabpool(8);
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

function results = writeCorrelationsToResults(morningDensity, eveningDensity, morningTimestamps10s, eveningTimestamps10s, goceMorning, goceEvening, latitude, timestamps10sFixed, timestamps1min, ae, results, modelName)

persistent goceTitleCellWritten;

[timestamps10s, order] = unique([morningTimestamps10s; eveningTimestamps10s]);
latitude = latitude(ismember(timestamps10sFixed, timestamps10s));
modelDensity = [morningDensity; eveningDensity];
goceDensity = [goceMorning; goceEvening];
modelDensity = modelDensity(order);
goceDensity = goceDensity(order);

firstDayIndices = timestamps10s < timestamps10s(1) + 86400;
modelMultiplier = mean(goceDensity(firstDayIndices) ./ modelDensity(firstDayIndices));
modelDensity = modelMultiplier * modelDensity;

timestampsNorth = timestamps10s(latitude > 80);
timesOrbAver = timestampsNorth(find(diff(timestampsNorth) > 45 * 60) + 1);
ind = ismember(timestamps10s, timesOrbAver);
modelSmooth = smooth(modelDensity, 541);
modelOrbAver = modelSmooth(ind);
goceSmooth = smooth(goceDensity, 541);
goceOrbAver = goceSmooth(ind);

[modelEfoldRising, modelEfoldFalling, modelRiseBegin] = giveEfoldingTimes(timesOrbAver, modelOrbAver, timestamps1min);
modelForXcorr = interp1(timesOrbAver, modelOrbAver, timestamps1min, 'linear', mean([modelOrbAver(1), modelOrbAver(end)]));
modelForXcorr = deleteDescendingParts(limitToDay(modelForXcorr, timestamps1min, modelRiseBegin));

[goceEfoldRising, goceEfoldFalling, goceRiseBegin] = giveEfoldingTimes(timesOrbAver, goceOrbAver, timestamps1min);
goceForXcorr = interp1(timesOrbAver, goceOrbAver, timestamps1min, 'linear', mean([goceOrbAver(1), goceOrbAver(end)]));
goceForXcorr = deleteDescendingParts(limitToDay(goceForXcorr, timestamps1min, goceRiseBegin));

aeSmooth = deleteDescendingParts(smooth(ae, 91));
aeForModel = limitToDay(aeSmooth, timestamps1min, modelRiseBegin);
aeForGoce = limitToDay(aeSmooth, timestamps1min, goceRiseBegin);
modelLag = giveMaxCrossCorrelation(modelForXcorr, aeForModel, 60);
goceLag = giveMaxCrossCorrelation(goceForXcorr, aeForGoce, 60);

equatorIndices = -30 < latitude & latitude < 30;
modelEquator = modelDensity(equatorIndices);
goceEquator = goceDensity(equatorIndices);

polarIndices = abs(latitude) > 45;
modelPolar = modelDensity(polarIndices);
gocePolar = goceDensity(polarIndices);

modelCorr = corr(modelDensity, goceDensity);
modelEqCorr = corr(modelEquator, goceEquator);
modelPolarCorr = corr(modelPolar, gocePolar);
modelOrbAverCorr = corr(modelOrbAver, goceOrbAver);

goceModelRatio = goceOrbAver ./ modelOrbAver;
modelMeanRatio = mean(goceModelRatio);
modelStdRatio = std(goceModelRatio);

[rowNum, ~] = size(results);
emptyCells = cellfun(@isempty,results);
[~, emptyColPositions] = find(emptyCells);
colNum = min(emptyColPositions);
if isempty(colNum)
    [~, colNum] = size(results);
    colNum = colNum + 1;
end

if isempty(goceTitleCellWritten)
    results{1, colNum} = 'Goce lag';
    results{1, colNum + 1} = 'Goce efold rising';
    results{1, colNum + 2} = 'Goce efold falling';
    goceTitleCellWritten = 1;
end 

    function isSame = stringCompare(cellString)
        isSame = strcmpi(cellString, 'Goce lag');
    end

goceColNum = find(cellfun(@stringCompare, results(1,:)));
if colNum == goceColNum
    results{rowNum, colNum} = goceLag;
    results{rowNum, colNum + 1} = goceEfoldRising;
    results{rowNum, colNum + 2} = goceEfoldFalling;
else
    colNum = colNum - 3;
end

results{1, colNum + 3} = [modelName, ' corr.'];
results{1, colNum + 4} = [modelName, ' eq.'];
results{1, colNum + 5} = [modelName, ' polar.'];
results{1, colNum + 6} = [modelName, ' O/M.'];
results{1, colNum + 7} = [modelName, ' std O/M.'];
results{1, colNum + 8} = [modelName, ' Orb. Aver. corr.'];
results{1, colNum + 9} = [modelName, ' lag'];
results{1, colNum + 10} = [modelName, ' efold rising'];
results{1, colNum + 11} = [modelName, ' efold falling'];

results{rowNum, colNum + 3} = modelCorr;
results{rowNum, colNum + 4} = modelEqCorr;
results{rowNum, colNum + 5} = modelPolarCorr;
results{rowNum, colNum + 6} = modelMeanRatio;
results{rowNum, colNum + 7} = modelStdRatio;
results{rowNum, colNum + 8} = modelOrbAverCorr;
results{rowNum, colNum + 9} = modelLag;
results{rowNum, colNum + 10} = modelEfoldRising;
results{rowNum, colNum + 11} = modelEfoldFalling;

end


% function results = writeCorrelationsToResults(morningMsisDensity, eveningMsisDensity, morningJbDensity, eveningJbDensity, morningAeProxy, eveningAeProxy, morningTimestamps10s, eveningTimestamps10s,...
%         morningDensityNoBg, eveningDensityNoBg, latitude, timestamps10sFixed, timestamps1min, ae, results)
% %
% 
% [timestamps10s, order] = unique([morningTimestamps10s; eveningTimestamps10s]);
% latitude = latitude(ismember(timestamps10sFixed, timestamps10s));
% aePredictedDensity = [morningAeProxy; eveningAeProxy];
% msisDensity = [morningMsisDensity; eveningMsisDensity];
% jbDensity = [morningJbDensity; eveningJbDensity];
% goceDensity = [morningDensityNoBg; eveningDensityNoBg];
% 
% aePredictedDensity = aePredictedDensity(order);
% msisDensity = msisDensity(order);
% jbDensity = jbDensity(order);
% goceDensity = goceDensity(order);
% 
% firstDayIndices = timestamps10s < timestamps10s(1) + 86400;
% msisMultiplier = mean(goceDensity(firstDayIndices) ./ msisDensity(firstDayIndices));
% aeMultiplier = mean(goceDensity(firstDayIndices) ./ aePredictedDensity(firstDayIndices));
% jbMultiplier = mean(goceDensity(firstDayIndices) ./ jbDensity(firstDayIndices));
% 
% msisDensity = msisMultiplier * msisDensity;
% aePredictedDensity = aeMultiplier * aePredictedDensity;
% jbDensity = jbMultiplier * jbDensity;
% 
% timestampsNorth = timestamps10s(latitude > 80);
% timesOrbAver = timestampsNorth(find(diff(timestampsNorth) > 45 * 60) + 1);
% ind = ismember(timestamps10s, timesOrbAver);
% msisSmooth = smooth(msisDensity, 541);
% aeIntegSmooth = smooth(aePredictedDensity, 541);
% jbSmooth = smooth(jbDensity, 541);
% goceSmooth = smooth(goceDensity, 541);
% msisOrbAver = msisSmooth(ind);
% aeOrbAver = aeIntegSmooth(ind);
% jbOrbAver = jbSmooth(ind);
% goceOrbAver = goceSmooth(ind);
% 
% [msisEfoldRising, msisEfoldFalling, msisRiseBegin] = giveEfoldingTimes(timesOrbAver, msisOrbAver, timestamps1min);
% [jbEfoldRising, jbEfoldFalling, jbRiseBegin] = giveEfoldingTimes(timesOrbAver, jbOrbAver, timestamps1min);
% [aeEfoldRising, aeEfoldFalling, aeRiseBegin] = giveEfoldingTimes(timesOrbAver, aeOrbAver, timestamps1min);
% [goceEfoldRising, goceEfoldFalling, goceRiseBegin] = giveEfoldingTimes(timesOrbAver, goceOrbAver, timestamps1min);
% 
% msisForXcorr = interp1(timesOrbAver, msisOrbAver, timestamps1min, 'linear', mean([msisOrbAver(1), msisOrbAver(end)]));
% jbForXcorr = interp1(timesOrbAver, jbOrbAver, timestamps1min, 'linear', mean([jbOrbAver(1), jbOrbAver(end)]));
% aeIntegForXcorr = interp1(timesOrbAver, aeOrbAver, timestamps1min, 'linear', mean([aeOrbAver(1), aeOrbAver(end)]));
% goceForXcorr = interp1(timesOrbAver, goceOrbAver, timestamps1min, 'linear', mean([goceOrbAver(1), goceOrbAver(end)]));
% 
% msisForXcorr = deleteDescendingParts(limitToDay(msisForXcorr, timestamps1min, msisRiseBegin));
% jbForXcorr = deleteDescendingParts(limitToDay(jbForXcorr, timestamps1min, jbRiseBegin));
% aeIntegForXcorr = deleteDescendingParts(limitToDay(aeIntegForXcorr, timestamps1min, aeRiseBegin));
% goceForXcorr = deleteDescendingParts(limitToDay(goceForXcorr, timestamps1min, goceRiseBegin));
% 
% aeSmooth = deleteDescendingParts(smooth(ae, 91));
% 
% aeForMsis = limitToDay(aeSmooth, timestamps1min, msisRiseBegin);
% aeForJb = limitToDay(aeSmooth, timestamps1min, jbRiseBegin);
% aeForAeInteg = limitToDay(aeSmooth, timestamps1min, aeRiseBegin);
% aeForGoce = limitToDay(aeSmooth, timestamps1min, goceRiseBegin);
% 
% msisLag = giveMaxCrossCorrelation(msisForXcorr, aeForMsis, 60);
% jbLag = giveMaxCrossCorrelation(jbForXcorr, aeForJb, 60);
% aeIntegLag = giveMaxCrossCorrelation(aeIntegForXcorr, aeForAeInteg, 60);
% goceLag = giveMaxCrossCorrelation(goceForXcorr, aeForGoce, 60);
% 
% equatorIndices = -30 < latitude & latitude < 30;
% msisEquator = msisDensity(equatorIndices);
% aeEquator = aePredictedDensity(equatorIndices);
% jbEquator = jbDensity(equatorIndices);
% goceEquator = goceDensity(equatorIndices);
% 
% polarIndices = abs(latitude) > 45;
% msisPolar = msisDensity(polarIndices);
% aePolar = aePredictedDensity(polarIndices);
% jbPolar = jbDensity(polarIndices);
% gocePolar = goceDensity(polarIndices);
% 
% aeModelCorr = corr(aePredictedDensity, goceDensity);
% msisCorr = corr(msisDensity, goceDensity);
% jbCorr = corr(jbDensity, goceDensity);
% msisEqCorr = corr(msisEquator, goceEquator);
% aeEqCorr = corr(aeEquator, goceEquator);
% jbEqCorr = corr(jbEquator, goceEquator);
% msisPolarCorr = corr(msisPolar, gocePolar);
% aePolarCorr = corr(aePolar, gocePolar);
% jbPolarCorr = corr(jbPolar, gocePolar);
% msisOrbAverCorr = corr(msisOrbAver, goceOrbAver);
% aeOrbAverCorr = corr(aeOrbAver, goceOrbAver);
% jbOrbAverCorr = corr(jbOrbAver, goceOrbAver);
% 
% % goceAeRatio = goceDensity ./ aePredictedDensity;
% % goceMsisRatio = goceDensity ./ msisDensity;
% % goceJbRatio = goceDensity ./ jbDensity;
% 
% goceAeRatio = goceOrbAver ./ aeOrbAver;
% goceMsisRatio = goceOrbAver ./ msisOrbAver;
% goceJbRatio = goceOrbAver ./ jbOrbAver;
% 
% aeModelMeanRatio = mean(goceAeRatio);
% msisMeanRatio = mean(goceMsisRatio);
% jbMeanRatio = mean(goceJbRatio);
% 
% aeModelStdRatio = std(goceAeRatio);
% msisStdRatio = std(goceMsisRatio);
% jbStdRatio = std(goceJbRatio);
% 
% [rowNum, ~] = size(results);
% emptyCells = cellfun(@isempty,results);
% [~, emptyColPositions] = find(emptyCells);
% colNum = min(emptyColPositions);
% if isempty(colNum)
%     [~, colNum] = size(results);
%     colNum = colNum + 1;
% end
% 
% results{1, colNum} = 'Msis corr.';
% results{1, colNum + 1} = 'AE Int. model corr.';
% results{1, colNum + 2} = 'JB2008 corr.';
% 
% results{1, colNum + 3} = 'Msis eq.';
% results{1, colNum + 4} = 'AE eq.';
% results{1, colNum + 5} = 'JB2008 eq.';
% 
% results{1, colNum + 6} = 'Msis polar.';
% results{1, colNum + 7} = 'AE polar.';
% results{1, colNum + 8} = 'JB2008 polar.';
% 
% results{1, colNum + 9} = 'Msis O/M.';
% results{1, colNum + 10} = 'AE Int. model O/M.';
% results{1, colNum + 11} = 'JB2008 O/M.';
% 
% results{1, colNum + 12} = 'Msis std O/M.';
% results{1, colNum + 13} = 'AE Int. model std O/M.';
% results{1, colNum + 14} = 'JB2008 std O/M.';
% 
% results{1, colNum + 15} = 'Msis Orb. Aver. corr.';
% results{1, colNum + 16} = 'AE Orb. Aver. corr.';
% results{1, colNum + 17} = 'JB2008 Orb. Aver. corr.';
% 
% results{1, colNum + 18} = 'Msis lag';
% results{1, colNum + 19} = 'AE Int. lag';
% results{1, colNum + 20} = 'JB2008 lag';
% 
% results{1, colNum + 21} = 'Msis efold rising';
% results{1, colNum + 22} = 'AE efold rising';
% results{1, colNum + 23} = 'JB2008 efold rising';
% 
% results{1, colNum + 24} = 'Msis efold falling';
% results{1, colNum + 25} = 'AE efold falling';
% results{1, colNum + 26} = 'JB2008 efold falling';
% 
% results{1, colNum + 27} = 'Goce lag';
% results{1, colNum + 28} = 'Goce efold rising';
% results{1, colNum + 29} = 'Goce efold falling';
% 
% results{rowNum, colNum} = msisCorr;
% results{rowNum, colNum + 1} = aeModelCorr;
% results{rowNum, colNum + 2} = jbCorr;
% 
% results{rowNum, colNum + 3} = msisEqCorr;
% results{rowNum, colNum + 4} = aeEqCorr;
% results{rowNum, colNum + 5} = jbEqCorr;
% 
% results{rowNum, colNum + 6} = msisPolarCorr;
% results{rowNum, colNum + 7} = aePolarCorr;
% results{rowNum, colNum + 8} = jbPolarCorr;
% 
% results{rowNum, colNum + 9} = msisMeanRatio;
% results{rowNum, colNum + 10} = aeModelMeanRatio;
% results{rowNum, colNum + 11} = jbMeanRatio;
% 
% results{rowNum, colNum + 12} = msisStdRatio;
% results{rowNum, colNum + 13} = aeModelStdRatio;
% results{rowNum, colNum + 14} = jbStdRatio;
% 
% results{rowNum, colNum + 15} = msisOrbAverCorr;
% results{rowNum, colNum + 16} = aeOrbAverCorr;
% results{rowNum, colNum + 17} = jbOrbAverCorr;
% 
% results{rowNum, colNum + 18} = msisLag;
% results{rowNum, colNum + 19} = aeIntegLag;
% results{rowNum, colNum + 20} = jbLag;
% 
% results{rowNum, colNum + 21} = msisEfoldRising;
% results{rowNum, colNum + 22} = aeEfoldRising;
% results{rowNum, colNum + 23} = jbEfoldRising;
% 
% results{rowNum, colNum + 24} = msisEfoldFalling;
% results{rowNum, colNum + 25} = aeEfoldFalling;
% results{rowNum, colNum + 26} = jbEfoldFalling;
% 
% results{rowNum, colNum + 27} = goceLag;
% results{rowNum, colNum + 28} = goceEfoldRising;
% results{rowNum, colNum + 29} = goceEfoldFalling;
% 
% end

function value = limitToDay(value, timestamps, ind)

[~,dayBeforeInd] = min(abs(timestamps - timestamps(ind) + 86400));
[~,dayAfterInd] = min(abs(timestamps - timestamps(ind) - 2*86400));
value = value(dayBeforeInd:dayAfterInd);

end

function value = deleteDescendingParts(value)

dVal = diff(value);
dVal = [dVal; dVal(end)];
value(dVal < 0) = 0;

end

function [densityIndexTimelag] = giveMaxCrossCorrelation(density, geomIndex, indicesInHour)
% plotCrossCorrelation(averagedDensityNoBg, ae)

maxLag = round(indicesInHour * 10.5);
[correlations, lags] = xcorr(density, geomIndex, maxLag, 'coeff');
correlations = correlations(lags >= 0);
lags = lags(lags >= 0);
densityIndexTimelag = lags(correlations == max(correlations)) / indicesInHour;

end
