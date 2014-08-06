function results = plotAndAnalyzeByHemisphere(firstDatenum, densityByLatitude, ae, crossingTimes, timestamps1min, timeOfDay, oneDegreeStep, plotFigures, results)
%

northIndices = (oneDegreeStep < 90 & oneDegreeStep > 45);
equatorIndices = (oneDegreeStep < 30 & oneDegreeStep > -30);
southIndices = (oneDegreeStep < -45 & oneDegreeStep > -90);

northernDensity = mean(densityByLatitude(:,northIndices), 2);
equatorDensity = mean(densityByLatitude(:,equatorIndices), 2);
southernDensity = mean(densityByLatitude(:,southIndices), 2);

northernError = std(densityByLatitude(:,northIndices), 0, 2);
equatorError = std(densityByLatitude(:,equatorIndices), 0, 2);
southernError = std(densityByLatitude(:,southIndices), 0, 2);

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

    northTimestamps = northTimestamps / secondsInDay + firstDatenum;
    equatorTimestamps = equatorTimestamps / secondsInDay + firstDatenum;
    southTimestamps = southTimestamps / secondsInDay + firstDatenum;
    
    figHandle = figure;
    hAx(1) = subplot(3,2,1);
    errorbar(northTimestamps, northAbsDiff, northernError);
    datetick('x', 'dd', 'keeplimits')
    ylabel('Density [10^{-11} kgm^{-3}]')
    title('Absolute difference, NH')
    grid on

    hAx(3) = subplot(3,2,3);
    errorbar(equatorTimestamps, equatorAbsDiff, equatorError, 'g');
    datetick('x', 'dd', 'keeplimits')
    ylabel('Density [10^{-11} kgm^{-3}]')
    title('Absolute difference, equator')
    grid on

    hAx(5) = subplot(3,2,5);
    errorbar(southTimestamps, southAbsDiff, southernError, 'r');
    datetick('x', 'dd', 'keeplimits')
    ylabel('Density [10^{-11} kgm^{-3}]')
    title('Absolute difference, SH')
    grid on

    hAx(2) = subplot(3,2,2);
    errorbar(northTimestamps, northRelDiff * 100, 100 * northernError/northCalmMean);
    datetick('x', 'dd', 'keeplimits')
    ylabel('%')
    title('Relative difference, NH')
    grid on

    hAx(4) = subplot(3,2,4);
    errorbar(equatorTimestamps, equatorRelDiff * 100, 100 * equatorError/equatorCalmMean, 'g');
    datetick('x', 'dd', 'keeplimits')
    ylabel('%')
    title('Relative difference, equator')
    grid on

    hAx(6) = subplot(3,2,6);
    errorbar(southTimestamps, southRelDiff * 100, 100 * southernError/southCalmMean, 'r');
    datetick('x', 'dd', 'keeplimits')
    ylabel('%')
    title('Relative difference, SH')
    grid on
    
    axis(hAx, 'tight')

    annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Absolute and relative density changes on hemispheres: ', timeOfDay], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

    tightfig(figHandle);
end

end

function [densityIndexTimelag] = giveMaxCrossCorrelation(density, geomIndex)
% plotCrossCorrelation(averagedDensityNoBg, ae)

maxLag = 60 * 24;
[correlations, lags] = xcorr(density, geomIndex, maxLag, 'coeff');
correlations = correlations(lags > 0);
lags = lags(lags > 0);
indicesInHour = 60;
densityIndexTimelag = lags(correlations == max(correlations)) / indicesInHour;

end