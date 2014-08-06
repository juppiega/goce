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