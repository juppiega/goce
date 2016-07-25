function [magneticLocalTime] = computeMagneticTime(magneticLongitude, doy, timestampsDatenum)
%

hours = (timestampsDatenum - floor(timestampsDatenum)) * 24;

subSolarLat = 23.5 * sin(0.0172 * doy - 1.405);
subSolarLon = 15 * (12 - hours);

altitude = ones(size(doy)) * 270e3;
[~, magneticSubSolarLon] = convertToMagneticCoordinates(subSolarLat, subSolarLon, altitude);

magneticLocalTime = 12 + (magneticLongitude - magneticSubSolarLon) / 15;
tooBigTimes = magneticLocalTime >= 24;
tooSmallTimes = magneticLocalTime < 0;
magneticLocalTime(tooBigTimes) = magneticLocalTime(tooBigTimes) - 24;
magneticLocalTime(tooSmallTimes) = magneticLocalTime(tooSmallTimes) + 24;

end