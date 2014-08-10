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

crossings = find(subtractedVal(1:end-1) .* subtractedVal(2:end) <= 0);

end