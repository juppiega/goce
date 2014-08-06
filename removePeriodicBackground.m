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