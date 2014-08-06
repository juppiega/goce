function [ae, ap, absB, akasofuEpsilon, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, ...
 morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength, firstDatenum] ...
 = variables(threshold, results)

% [ae, averagedDensity, timestamps] = readAEDensityAndTimestamps( AEFilename, densityFilename )

if exist('GoceVariables.mat', 'file') == 2
    load('GoceVariables.mat')
else
    fprintf(2, '%s', 'Warning: no GoceVariables.mat found in Matlab PATH, now attempting to create new one -> readFiles.m');
end

%compareGoceDensityToMsis(measuredDensity, msisDensityVariableAlt, ae, timestampsAeDatenum, timestampsDensityDatenum, results);

intervalsOfInterest = findInterestingIntervals(ae, timestampsAeDatenum, epsilonQualityFlag, timestampsEpsilonDatenum, densityNoBg, timestampsDensityDatenum, threshold);

[morningTimestamps10s, morningMagneticLatitude, morningDensityNoBg, morningMsisDensity, eveningTimestamps10s, ...
    eveningMagneticLatitude, eveningDensityNoBg, eveningMsisDensity] = ...
    splitBySolarTime(timestamps10sFixed, magneticLatitude, densityNoBg, msisDensity270km, solarTime);

[ae, ap, absB, akasofuEpsilon, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1min, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = sliceToInterestingIntervals(ae, ap, absB, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed,...
 density3h, morningMagneticLatitude, eveningMagneticLatitude, timestamps10sFixed, timestamps1min, timestampsAbsB, akasofuEpsilon, timestampsDensityDatenum, intervalsOfInterest);

end


function [morningTimestamps10s, morningMagneticLatitude, morningDensityNoBg, morningMsisDensity, eveningTimestamps10s, ...
    eveningMagneticLatitude, eveningDensityNoBg, eveningMsisDensity] = ...
    splitBySolarTime(timestamps10s, magneticLatitude, densityNoBg, msisDensity, solarTime)
%

morningIndices = find(solarTime <= 12);
eveningIndices = find(solarTime > 12);

morningTimestamps10s = timestamps10s(morningIndices);
morningMagneticLatitude = magneticLatitude(morningIndices);
morningDensityNoBg = densityNoBg(morningIndices);
morningMsisDensity = msisDensity(morningIndices);

eveningTimestamps10s = timestamps10s(eveningIndices);
eveningMagneticLatitude = magneticLatitude(eveningIndices);
eveningDensityNoBg = densityNoBg(eveningIndices);
eveningMsisDensity = msisDensity(eveningIndices);

end


function [ae, ap, absB, akasofuEpsilon, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1minOut, timestampsAbsB, ...
 timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed, density3h, timestampsDatenum, morningMagneticLatitude, eveningMagneticLatitude, cellArrayLength] ...
 = sliceToInterestingIntervals(ae, ap, absB, averagedDensityNoBg, morningDensityNoBg, eveningDensityNoBg, ...
 morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s, timestamps1minFixed, timestampsEpsilon, timestamps3h, timestamps3hFixed,...
 density3h, morningMagneticLatitude, eveningMagneticLatitude, timestamps10sFixed, timestamps1min, timestampsAbsB, akasofuEpsilon, timestampsDatenum, intervalsOfInterest)
% 

aeTemp = ae; apTemp = ap;
averagedDensityNoBgTemp = averagedDensityNoBg; morningDensityNoBgTemp = morningDensityNoBg; 
eveningDensityNoBgTemp = eveningDensityNoBg; morningMsisDensityTemp = morningMsisDensity;
eveningMsisDensityTemp = eveningMsisDensity; morningTimestamps10sTemp = morningTimestamps10s; 
eveningTimestamps10sTemp = eveningTimestamps10s; timestamps1minFixedTemp = timestamps1minFixed; 
timestamps3hTemp = timestamps3h; timestamps3hFixedTemp = timestamps3hFixed; density3hTemp = density3h; morningMagneticLatitudeTemp = morningMagneticLatitude;
eveningMagneticLatitudeTemp = eveningMagneticLatitude; timestampsEpsilonTemp = timestampsEpsilon; akasofuEpsilonTemp = akasofuEpsilon;
timestampsAbsBTemp = timestampsAbsB; absBTemp = absB; timestampsDatenumTemp = timestampsDatenum; 

ae = {}; ap = {};  averagedDensityNoBg = {}; morningDensityNoBg = {};
eveningDensityNoBg = {}; morningMsisDensity = {}; eveningMsisDensity = {}; morningTimestamps10s = {};
eveningTimestamps10s = {}; timestamps1minFixed = {}; timestamps3h = {}; density3h = {}; morningMagneticLatitude = {};
eveningMagneticLatitude = {}; timestampsAbsB = {}; akasofuEpsilon = {}; timestampsEpsilon = {}; absB = {}; timestampsDatenum = {};
timestamps3hFixed = {};

cellArrayLength = length(intervalsOfInterest(:,1));

for i = 1:cellArrayLength
    beginIndex = intervalsOfInterest(i,1);
    endIndex = intervalsOfInterest(i,2);
    timestamps1minOut{i} = timestamps1min(beginIndex:endIndex);
    ae{i} = aeTemp(beginIndex:endIndex);
    
    threeHinSec = 3 * 60 * 60;
    [~, beginIndex3h] = min(abs(timestamps3hTemp - threeHinSec * round(timestamps1min(beginIndex) / threeHinSec)));
    [~, endIndex3h] = min(abs(timestamps3hTemp - threeHinSec * round(timestamps1min(endIndex) / threeHinSec)));
    ap{i} = apTemp(beginIndex3h:endIndex3h);
    timestamps3h{i} = timestamps3hTemp(beginIndex3h:endIndex3h);
    
    [~, beginIndex3hFixed] = min(abs(timestamps3hFixedTemp - threeHinSec * round(timestamps1min(beginIndex) / threeHinSec)));
    [~, endIndex3hFixed] = min(abs(timestamps3hFixedTemp - threeHinSec * round(timestamps1min(endIndex) / threeHinSec)));
    density3h{i} = density3hTemp(beginIndex3hFixed:endIndex3hFixed);
    timestamps3hFixed{i} = timestamps3hFixedTemp(beginIndex3hFixed:endIndex3hFixed);
    
    [~, beginIndexAbsB] = min(abs(timestampsAbsBTemp - timestamps1min(beginIndex)));
    [~, endIndexAbsB] = min(abs(timestampsAbsBTemp - timestamps1min(endIndex)));
    absB{i} = absBTemp(beginIndexAbsB:endIndexAbsB);
    timestampsAbsB{i} = timestampsAbsBTemp(beginIndexAbsB:endIndexAbsB);
    
    [~, beginIndexEpsilon] = min(abs(timestampsEpsilonTemp - timestamps1min(beginIndex)));
    [~, endIndexEpsilon] = min(abs(timestampsEpsilonTemp - timestamps1min(endIndex)));
    akasofuEpsilon{i} = akasofuEpsilonTemp(beginIndexEpsilon:endIndexEpsilon);
    timestampsEpsilon{i} = timestampsEpsilonTemp(beginIndexEpsilon:endIndexEpsilon);
    
    [~, beginIndex1min] = min(abs(timestamps1minFixedTemp - timestamps1min(beginIndex)));
    [~, endIndex1min] = min(abs(timestamps1minFixedTemp - timestamps1min(endIndex)));
    averagedDensityNoBg{i} = averagedDensityNoBgTemp(beginIndex1min:endIndex1min);
    timestamps1minFixed{i} = timestamps1minFixedTemp(beginIndex1min:endIndex1min);
    
    [~, beginIndexMorning10s] = min(abs(morningTimestamps10sTemp - timestamps1min(beginIndex)));
    [~, endIndexMorning10s] = min(abs(morningTimestamps10sTemp - timestamps1min(endIndex)));   
    [~, beginIndexEvening10s] = min(abs(eveningTimestamps10sTemp - timestamps1min(beginIndex)));
    [~, endIndexEvening10s] = min(abs(eveningTimestamps10sTemp - timestamps1min(endIndex)));
    morningMsisDensity{i} = morningMsisDensityTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningMsisDensity{i} = eveningMsisDensityTemp(beginIndexEvening10s:endIndexEvening10s);
    morningDensityNoBg{i} = morningDensityNoBgTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningDensityNoBg{i} = eveningDensityNoBgTemp(beginIndexEvening10s:endIndexEvening10s);
    morningTimestamps10s{i} = morningTimestamps10sTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningTimestamps10s{i} = eveningTimestamps10sTemp(beginIndexEvening10s:endIndexEvening10s);
    morningMagneticLatitude{i} = morningMagneticLatitudeTemp(beginIndexMorning10s:endIndexMorning10s);
    eveningMagneticLatitude{i} = eveningMagneticLatitudeTemp(beginIndexEvening10s:endIndexEvening10s);
    
    [~, beginIndex10s] = min(abs(timestamps10sFixed - timestamps1min(beginIndex)));
    [~, endIndex10s] = min(abs(timestamps10sFixed - timestamps1min(endIndex)));
    timestampsDatenum{i} = timestampsDatenumTemp(beginIndex10s:endIndex10s);
end

end