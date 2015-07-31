function [efoldRising, efoldFalling, risingBegin] = giveEfoldingTimes(timestampsDensity, densityNoBg, timestamps1min)

densityNoBg = interp1(timestampsDensity, densityNoBg, timestamps1min, 'linear', 0);
densityNoBg = removeLeadingAndTrailingZeros(densityNoBg);

smoothRange = 630;
densityNoBg = smooth(densityNoBg, smoothRange);

[risingBegin, risingEnd, fallingBegin, fallingEnd] = findRisingAndFallingLimbs(densityNoBg);

efoldRising = computeEfold(timestamps1min, risingBegin, risingEnd, densityNoBg);
efoldFalling = computeEfold(timestamps1min, fallingBegin, fallingEnd, densityNoBg);

end

function efoldTime = computeEfold(timestamps1min, beginInd, endInd, densityNoBg)

e = 2.71828183;
ind = beginInd:endInd; 
modifTime = (timestamps1min(ind) - timestamps1min(beginInd)) / 3600;
% startB = sign(densityNoBg(endInd) - densityNoBg(beginInd)) * 1/20;
% startC = min([densityNoBg(beginInd), densityNoBg(endInd)])*0.9;
% startA = densityNoBg(beginInd) - startC;
% startVals = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', [startA, startB, startC]);
% expFit = fittype('a*exp(b*t)+c', 'independent', 't', 'coefficients', {'a','b','c'}, 'options', startVals);
% f = fit(modifTime, densityNoBg(ind), expFit);
% 
% efoldTime = abs(1/f.b);

signMult = sign(densityNoBg(endInd) - densityNoBg(beginInd));
valDiff = abs(densityNoBg(endInd) - densityNoBg(beginInd));
if signMult > 0
    targetVal = densityNoBg(beginInd) + (1-1/e) * valDiff;
else
    targetVal = densityNoBg(endInd) + 1/e * valDiff;
end
[~,timeInd] = min(abs(densityNoBg(ind) - targetVal));
timeInd = timeInd(1);
efoldTime = modifTime(timeInd);

%if efoldTime > 100
%     plot(modifTime, densityNoBg(ind));
%     hold all;
%     line([efoldTime, efoldTime], get(gca,'ylim'), 'color', 'k')
%     hold off
%end

end

function value = removeLeadingAndTrailingZeros(value)
%

firstNonZero = find(value > 0, 1, 'first');
firstNonZeroValue = value(firstNonZero);
lastNonZero = find(value > 0, 1, 'last');
lastNonZeroValue = value(lastNonZero);

value(1:firstNonZero) = firstNonZeroValue;
value(lastNonZero:end) = lastNonZeroValue;

end

function [risingBegin, risingEnd, fallingBegin, fallingEnd] = findRisingAndFallingLimbs(value)
%

dValue = diff(value);
dValue = smooth(dValue, 720);
threshold = 0;
[~, peakApproxBegin, peakApproxEnd] = limitToNearPeak(value, 'noSmooth', 'median');

dValueRising = dValue(1:peakApproxEnd - 1);
risingIndices = find(dValueRising > threshold);
risingIndices = ismember(1:length(dValueRising), risingIndices);
edgeArray = diff([0 risingIndices 0]);
indexLimitsAboveZero = [find(edgeArray > 0)' (find(edgeArray < 0)-1)'];

[~, risingBeginApprox, risingEndApprox] = findLargestSum(dValueRising, indexLimitsAboveZero);
risingEndApprox = risingEndApprox + 1;
[~,risingBegin] = min(value(risingBeginApprox:risingEndApprox));
[~,risingEnd] = max(value(risingBeginApprox:risingEndApprox));
risingBegin = risingBegin + risingBeginApprox - 1;
risingEnd = risingEnd + risingBeginApprox - 1;

dValueFalling = dValue;
fallingIndices = find(dValueFalling < threshold);
fallingIndices = ismember(1:length(dValueFalling), fallingIndices);
fallingIndices(1:risingEnd - 1) = 0;
edgeArray = diff([0 fallingIndices 0]);
indexLimitsBelowZero = [find(edgeArray > 0)' (find(edgeArray < 0)-1)'];

try
    [~, fallingBeginApprox, fallingEndApprox] = findLargestSum(-1 * dValueFalling, indexLimitsBelowZero);
catch
    [~, fallingBeginApprox, fallingEndApprox] = findLargestSum(-1 * dValueFalling, indexLimitsBelowZero);
end
fallingEndApprox = fallingEndApprox + 1;
[~,fallingBegin] = max(value(fallingBeginApprox:fallingEndApprox));
[~,fallingEnd] = min(value(fallingBeginApprox:fallingEndApprox));
fallingBegin = fallingBegin + fallingBeginApprox - 1;
fallingEnd = fallingEnd + fallingBeginApprox - 1;

% figure;
% plot(1:length(value), value);
% ylims = get(gca, 'ylim');
% line([risingBegin risingBegin], [ylims(1) ylims(2)], 'LineStyle', '--')
% line([risingEnd risingEnd], [ylims(1) ylims(2)], 'LineStyle', '--')
% line([fallingBegin fallingBegin], [ylims(1) ylims(2)], 'LineStyle', '--')
% line([fallingEnd fallingEnd], [ylims(1) ylims(2)], 'LineStyle', '--')

end

function [largestSum, largestSumBeginIndex, largestSumEndIndex] = findLargestSum(value, indices)
%

intervalBegin = indices(:,1);
intervalEnd = indices(:,2);

largestSum = 0;

for i = 1:length(intervalBegin)
    intervalIndices = intervalBegin(i):intervalEnd(i);
    intervalSum = sum(value(intervalIndices));
    
    if intervalSum > largestSum
        largestSumBeginIndex = intervalBegin(i);
        largestSumEndIndex = intervalEnd(i);
        largestSum = intervalSum;
    end
end

end