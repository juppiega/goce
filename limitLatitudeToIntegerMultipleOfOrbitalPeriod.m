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