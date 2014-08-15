function results = plotAndAnalyzeChangesByOrbit(firstDatenum, densityNoBg, magneticLatitude, averagedDensityNoBg, timestamps1minFixed, ...
    timestamps10s, timeOfDay, plotFigures, results)
% 

[limitedTimestamps, limitedLatitude, limitedDensity, orbits, minAllowedLatitude, maxAllowedLatitude] = ...
    giveInterpolationQualifyingOrbits(timestamps10s, magneticLatitude, densityNoBg);

results = loopAroundOrbits(orbits, timestamps1minFixed, limitedTimestamps, limitedLatitude, averagedDensityNoBg, limitedDensity, ...
    plotFigures, timeOfDay, firstDatenum, minAllowedLatitude, maxAllowedLatitude, results);

end


function [beginOrbit, beginStormOrbit, endOrbit, endColorOrbit, endColorMapOrbit, calmOrbits] = findBeginAndEndOrbits(orbits, timestamps1minFixed, limitedTimestamps, averagedDensityNoBg)
%

calmOrbits = 5;

averagedDensityNoBg = averagedDensityNoBg(ismember(timestamps1minFixed, limitedTimestamps));
timestamps1minFixed = timestamps1minFixed(ismember(timestamps1minFixed, limitedTimestamps));
[~, peakBeginIndex, peakEndIndex] = limitToNearPeak(averagedDensityNoBg, 'noSmooth', 'mean');
peakBeginIndex = find(limitedTimestamps == timestamps1minFixed(peakBeginIndex));
peakEndIndex = find(limitedTimestamps == timestamps1minFixed(peakEndIndex));

beginStormOrbit = find(orbits(:,1) <= peakBeginIndex, 1, 'last');
if beginStormOrbit > calmOrbits
    beginOrbit = beginStormOrbit - calmOrbits;
else
    beginOrbit = 1;
end

endOrbit = find(orbits(:,2) >= peakEndIndex, 1, 'first');
if endOrbit < length(orbits(:,2)) - calmOrbits
    endOrbit = endOrbit + calmOrbits;
else
    endOrbit = length(orbits(:,2));
end

maxNumOfColorOrbits = 10;
extraColorMapOrbits = 5;

if endOrbit - beginOrbit - 2 * calmOrbits > maxNumOfColorOrbits + extraColorMapOrbits
    endColorOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits - 1;
    endColorMapOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits + extraColorMapOrbits - 1;
elseif endOrbit - beginOrbit - 2 * calmOrbits > maxNumOfColorOrbits
    endColorOrbit = beginOrbit + calmOrbits + maxNumOfColorOrbits - 1;
    endColorMapOrbit = endOrbit - calmOrbits;
else
    endColorOrbit = endOrbit;
    endColorMapOrbit = endOrbit;
end

end

function [limitedTimestamps, limitedLatitude, limitedDensity, orbits, minAllowedLatitude, maxAllowedLatitude] = ...
    giveInterpolationQualifyingOrbits(timestamps10s, magneticLatitude, densityNoBg)
%

[latBeginIndex, latEndIndex] = limitLatitudeToIntegerMultipleOfOrbitalPeriod(magneticLatitude);
limitedLatitude = magneticLatitude(latBeginIndex:latEndIndex);
limitedTimestamps = timestamps10s(latBeginIndex:latEndIndex);
limitedDensity = densityNoBg(latBeginIndex:latEndIndex);
[~, exactOrbitIndices] = splitIntoOrbits(limitedLatitude, limitedTimestamps);
limitedLatitude = limitedLatitude(exactOrbitIndices);
limitedTimestamps = limitedTimestamps(exactOrbitIndices);
limitedDensity = limitedDensity(exactOrbitIndices);
[orbits, ~] = splitIntoOrbits(limitedLatitude, limitedTimestamps);

[minAllowedLatitude, maxAllowedLatitude] = findInterpolationLimits(limitedLatitude);

if satelliteIsGoingSouth(limitedLatitude)
    orbitsEndingTooNorth = find(limitedLatitude(orbits(:,2)) > minAllowedLatitude);
    orbitsBeginningTooSouth = find(limitedLatitude(orbits(:,1)) < maxAllowedLatitude);
    orbitsToDelete = unique(vertcat(orbitsBeginningTooSouth, orbitsEndingTooNorth));
else
    orbitsEndingTooSouth = find(limitedLatitude(orbits(:,2)) < maxAllowedLatitude);
    orbitsBeginningTooNorth = find(limitedLatitude(orbits(:,1)) > minAllowedLatitude);
    orbitsToDelete = unique(vertcat(orbitsBeginningTooNorth, orbitsEndingTooSouth));
end

newIndices = 1:length(limitedLatitude);
for i = 1:length(orbitsToDelete)
    newIndices = setdiff(newIndices, (orbits(orbitsToDelete(i), 1) : orbits(orbitsToDelete(i), 2)));
end
limitedTimestamps = limitedTimestamps(newIndices);
limitedLatitude = limitedLatitude(newIndices);
limitedDensity = limitedDensity(newIndices);
[orbits, ~] = splitIntoOrbits(limitedLatitude, limitedTimestamps);

end

function results = loopAroundOrbits(orbits, timestamps1minFixed, limitedTimestamps, limitedLatitude, averagedDensityNoBg, ...
    limitedDensity, plotFigures, timeOfDay, firstDatenum, minAllowedLatitude, maxAllowedLatitude, results)
%

persistent densityByOrbitFigHandle
persistent densityByOrbitAxesHandle

if plotFigures ~= 0
    if ~isempty(strfind(lower(timeOfDay), 'morning')); densityByOrbitFigHandle = figure; subplotNum = 1; else subplotNum = 2; end
    figure(densityByOrbitFigHandle);
    densityByOrbitAxesHandle(subplotNum) = subplot(2,1,subplotNum);
    hold all;
end

linehandles = [];
relativeResidues = nan(size(limitedLatitude));
TADdensity = nan(size(limitedLatitude));

[beginOrbit, beginStormOrbit, endOrbit, endColorOrbit, endColorMapOrbit, calmOrbits] = findBeginAndEndOrbits(orbits, timestamps1minFixed, limitedTimestamps, averagedDensityNoBg);

loopOrbits = beginOrbit:endOrbit;
TADplotIndices = nan(size(limitedLatitude));
for i = 1:length(loopOrbits)
    indices = orbits(loopOrbits(i),1) : orbits(loopOrbits(i),2);
    smoothedDensity150s = smooth(limitedDensity(indices), 15);
    relativeResidues(indices) = (limitedDensity(indices) - smoothedDensity150s) ./ smoothedDensity150s;
    if (loopOrbits(i) <= endColorOrbit || loopOrbits(i) > endOrbit - calmOrbits) && plotFigures ~= 0
        h = plot(limitedLatitude(indices), smoothedDensity150s, 'LineWidth', 2);
        linehandles = [linehandles h];
    end
    
    if loopOrbits(i) <= endColorMapOrbit
        TADplotIndices(indices) = indices;            
        smoothedDensity11800km = smooth(limitedDensity(indices), 151);
        smoothedDensity2600km = smooth(limitedDensity(indices), 33); 
        TADdensity(indices) = (smoothedDensity2600km - smoothedDensity11800km) ./ smoothedDensity11800km;
    end    
end

TADforSpectralDensity = TADdensity;
TADforSpectralDensity = TADforSpectralDensity - nanmean(TADdensity);
TADforSpectralDensity(isnan(TADforSpectralDensity)) = 0;
[intervalPower, maxPower, scaleAtMax, averageScale] = estimatePowerSpectralDensity(TADforSpectralDensity, ...
    orbits(beginStormOrbit:endColorMapOrbit,:), timeOfDay, plotFigures);

[results, magneticLatitudeForResidue, relativeResidues] = writeSmallWaveAndTADdataToResults(limitedLatitude, relativeResidues, intervalPower, maxPower, scaleAtMax, averageScale, timeOfDay, results);

if plotFigures ~= 0

    plotResiduePlot(magneticLatitudeForResidue, relativeResidues, timeOfDay)
    plotTADcolormap(limitedTimestamps, limitedLatitude, TADdensity, TADplotIndices, minAllowedLatitude, maxAllowedLatitude, ...
        orbits, beginOrbit, calmOrbits, endColorMapOrbit, firstDatenum, timeOfDay)
    plotDensityByOrbit(linehandles, calmOrbits, timeOfDay, densityByOrbitAxesHandle, densityByOrbitFigHandle)
end

end

function plotDensityByOrbit(linehandles, calmOrbits, timeOfDay, densityByOrbitAxesHandle, densityByOrbitFigHandle)
%
figure(densityByOrbitFigHandle);

stormOrbits = (calmOrbits + 1 : length(linehandles) - calmOrbits);
[~,cmap] = cmapline('colormap',jet,'filled', 'lines', linehandles(stormOrbits)');
colormap(cmap)
colorbar

set(linehandles(1:calmOrbits), 'LineStyle', '--')
set(linehandles(1:calmOrbits), 'Color', 'k')
endLineHandles = length(linehandles) - calmOrbits + 1 : length(linehandles);
set(linehandles(endLineHandles), 'LineStyle', '-')
set(linehandles(endLineHandles), 'Color', 'k')
hold off;
xlabel('Geomagnetic latitude')
ylabel('Density at 270 km')
xlim([-90 90])
title(timeOfDay)
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '150-s smoothed density along orbit', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
if ~isempty(strfind(lower(timeOfDay), 'evening'))
   linkaxes(densityByOrbitAxesHandle) 
end

end

function plotResiduePlot(magneticLatitudeForResidue, relativeResidues, timeOfDay)
%

persistent residueFigHandle
persistent residueAxisHandle

confIntervalStep = 10;
confIntervalX = [];
confIntervalMean = [];
confIntervalLower = [];
confIntervalUpper = [];
for i = -90:confIntervalStep:90 - confIntervalStep
    indicesInInterval = find(magneticLatitudeForResidue < i + 5 & magneticLatitudeForResidue >= i);
    if ~isempty(indicesInInterval)
        confIntervalX = [confIntervalX mean([i (i + confIntervalStep)])];
        confIntervalMean = [confIntervalMean mean(relativeResidues(indicesInInterval))];
        [~, ~, ci, ~] = ttest(relativeResidues(indicesInInterval));        
        confIntervalLower = [confIntervalLower ci(1)];
        confIntervalUpper = [confIntervalUpper ci(2)];
    end
end

if ~isempty(strfind(lower(timeOfDay), 'morning')); 
    residueFigHandle = figure('Color','w');
    residueAxisHandle = [];
    confIntervalColor = 'b'; 
    axisNum = 1;
    %mlock;
else
    confIntervalColor = 'r';
    axisNum = 2;
    %munlock;
end
figure(residueFigHandle);
residueAxisHandle = [residueAxisHandle ciplot(confIntervalLower, confIntervalUpper, confIntervalX, confIntervalColor)];
xlim([min(confIntervalX) max(confIntervalX)]);
hold all;
plot(confIntervalX, confIntervalMean, 'k--');
hold off;
xlabel('Geomagnetic latitude')
ylabel('Relative residues 95 % confidence')
title('100-1000 km variations')
if axisNum == 1
    hold on
else
    hold off
    legend(residueAxisHandle, '<-Morning', 'Evening->');
end

end


function plotTADcolormap(limitedTimestamps, limitedLatitude, TADdensity, TADplotIndices, minAllowedLatitude, maxAllowedLatitude, orbits, ...
    beginOrbit, calmOrbits, endColorMapOrbit, firstDatenum, timeOfDay)
%

[timeMatrix, latitudeMatrix, densityMatrix] = advancedScatteredInterpolate(limitedTimestamps, limitedLatitude, TADdensity, TADplotIndices,...
    minAllowedLatitude, maxAllowedLatitude);

secondsInDay = 24 * 60 * 60;
timeMatrix = timeMatrix / secondsInDay + firstDatenum;
referenceDay = datestr(min(timeMatrix(:)), 'mmmm dd, yyyy');
timeMatrix = timeMatrix - datenum(referenceDay, 'mmmm dd, yyyy');

limitedTimestamps = limitedTimestamps / secondsInDay + firstDatenum;
limitedTimestamps = limitedTimestamps - datenum(referenceDay, 'mmmm dd, yyyy');

firstOrbitTrackInd = orbits(beginOrbit + calmOrbits, 1);
lastOrbitTrackInd = orbits(endColorMapOrbit,2);
orbitTrackLatitude = limitedLatitude(firstOrbitTrackInd:lastOrbitTrackInd);
orbitTrackTime = limitedTimestamps(firstOrbitTrackInd:lastOrbitTrackInd);
orbitTrackPlotHeight = ones(size(orbitTrackLatitude)) * max(densityMatrix(:));
figure;
surf(timeMatrix, latitudeMatrix, densityMatrix,'EdgeColor', 'none');
view(2);
hold on;
h = plot3(orbitTrackTime, orbitTrackLatitude, orbitTrackPlotHeight, '.k');
set(h, 'MarkerSize', 3)
ylabel('Magnetic Latitude')
title(['1300-5200km changes (TADs) [(2600 km smooth - 10400 km smooth) / 10400 km smooth] ', timeOfDay])
colorbar
%colormap jet(500)
ylim([minAllowedLatitude maxAllowedLatitude]);
xlim([min(timeMatrix(:)) max(timeMatrix(:))]);
xlabel(['Days since the UTC beginning of ', referenceDay])

end

function [timeMatrix, latitudeMatrix, densityMatrix] = advancedScatteredInterpolate(limitedTimestamps, limitedLatitude, TADplot, TADplotIndices,...
    minAllowedLatitude, maxAllowedLatitude)
%

TADplot = TADplot(~isnan(TADplot));
indices = TADplotIndices(~isnan(TADplotIndices));
TADlatitude = limitedLatitude(indices);
TADtime = limitedTimestamps(indices);

oneQuarterDegreeStep = minAllowedLatitude:0.25:maxAllowedLatitude;

for i = 1:length(oneQuarterDegreeStep)
    regriddedTime(:,i) = latitudeCrossingTimes(TADlatitude, TADtime, oneQuarterDegreeStep(i)); 
end

regriddedDensity = interp1(TADtime, TADplot, regriddedTime, 'spline');

numOfOrbits = length(regriddedTime(:,1));
numOfValuesInOrbit = length(regriddedTime(1,:));
for i = 1:numOfValuesInOrbit
    timeThisLatitude = regriddedTime(:,i);
    densityThisLatitude = regriddedDensity(:,i);

    tInterp = interp1(1:numOfOrbits, timeThisLatitude, 1:1/60:numOfOrbits);
    interpolatedDensity = interp1(timeThisLatitude, densityThisLatitude, tInterp, 'spline');
    thisLatitude = ones(length(tInterp), 1) * oneQuarterDegreeStep(i);
    latitudeMatrix(:,i) = thisLatitude;
    densityMatrix(:,i) = interpolatedDensity;
    timeMatrix(:,i) = tInterp;
end

indicesToRemove = findMatrixIndicesInDatagap(timeMatrix, TADtime);
timeMatrix(indicesToRemove) = nan(1);
latitudeMatrix(indicesToRemove) = nan(1);
densityMatrix(indicesToRemove) = 0;

end

function [results, magneticLatitudeForResidue, relativeResidues] = writeSmallWaveAndTADdataToResults(limitedLatitude, ...
    relativeResidues, intervalPower, maxPower, scaleAtMax, averageScale, timeOfDay, results)
%

magneticLatitudeForResidue = limitedLatitude(~isnan(relativeResidues));
relativeResidues = relativeResidues(~isnan(relativeResidues));

stdSouth = std(relativeResidues(magneticLatitudeForResidue < 0));
stdNorth = std(relativeResidues(magneticLatitudeForResidue > 0));

[rowNum, ~] = size(results);
emptyCells = cellfun(@isempty,results);
[~, emptyColPositions] = find(emptyCells);
colNum = min(emptyColPositions);
if rowNum == 2
    colNum = length(results(rowNum,:)) + 1;
    results{1, colNum}     = ['Std SH ', timeOfDay];
    results{1, colNum + 1} = ['Std NH ', timeOfDay];
    
    results{1, colNum + 2} = ['TAD interval power ', timeOfDay];
    results{1, colNum + 3} = ['TAD max power  ', timeOfDay];
    results{1, colNum + 4} = ['TAD scale at max ', timeOfDay];
    results{1, colNum + 5} = ['TAD mean scale ', timeOfDay];
end
results{rowNum, colNum}     = stdSouth;
results{rowNum, colNum + 1} = stdNorth;

results{rowNum, colNum + 2} = intervalPower;
results{rowNum, colNum + 3} = maxPower;
results{rowNum, colNum + 4} = scaleAtMax;
results{rowNum, colNum + 5} = averageScale;

end

function [intervalPower, maxPower, scaleAtMax, averageScale] = estimatePowerSpectralDensity(filteredDensity, orbits, timeOfDay, plotFigures)
%

persistent spectralDensityFigHandle

samplingFreq = 6;
if mod(length(filteredDensity), 2) == 0
    numOfSamples = length(filteredDensity); % Use fixed length for all FFTs
else
    numOfSamples = length(filteredDensity) - 1;
end
orbitsToRemove = mod(length(orbits(:,1)), 3);
if orbitsToRemove > 0
    orbits(end - orbitsToRemove + 1 : end,:) = [];
end

loopIteration = 1;
for i = 1:3:length(orbits(:,1))
    beginIndex = orbits(i,1);
    endIndex = orbits(i + 2, 2);
    densityPsd(:,loopIteration) = (1 / (samplingFreq * numOfSamples)) * ...
        abs(fft(filteredDensity(beginIndex:endIndex), numOfSamples)) .^ 2;
    loopIteration = loopIteration + 1;
end

densityPsd = mean(densityPsd, 2);
densityPsd = densityPsd(1 : numOfSamples / 2 + 1);
densityPsd(2:end - 1) = 2 * densityPsd(2:end - 1);

freq = 0 : samplingFreq / numOfSamples : samplingFreq / 2;
periodInMinutes = 1 ./ freq;

earthCircumferenceAtGoceOrbit = 2 * pi * (6371 + 270);
goceOrbitPeriod = 89;
speed = earthCircumferenceAtGoceOrbit / goceOrbitPeriod;
periodInKm = speed * periodInMinutes;

lowerBound = 0.1 * round(min(periodInKm) / 0.1);
upperBound = 14000;

interval = periodInKm >= lowerBound & periodInKm <= upperBound;
periodInInterval = periodInKm(interval);
psdInInterval = densityPsd(interval);
cumPower = cumtrapz(densityPsd);
intervalBegin = find(interval, 1, 'first');
intervalEnd = find(interval, 1, 'last');

% Power in a frequency band is twice the contribution of power from zero frequency to interval
% begin frequency (largest scale!) plus once the power between largest and
% smallest scales. Keep in mind we have multiplied the original full spectrum by two.
intervalPower = cumPower(intervalBegin) + 0.5 * (cumPower(intervalEnd) - cumPower(intervalBegin));
[maxPower, maxIndex] = max(psdInInterval);
scaleAtMax = periodInInterval(maxIndex);
averageScale = sum(periodInInterval' .* psdInInterval) / sum(psdInInterval);

if plotFigures ~= 0

    if ~isempty(strfind(lower(timeOfDay), 'morning')); 
        spectralDensityFigHandle = figure;
        lineColor = 'b'; 
        axisNum = 1;
    else
        lineColor = 'r';
        axisNum = 2;
    end
    
    figure(spectralDensityFigHandle);
    plot(periodInKm, densityPsd, lineColor);
    xlim([lowerBound upperBound])
    set(gca, 'XTick', [lowerBound 2000:2000:14000])
    title(['TAD Power Spectral Density ', timeOfDay])
    ylabel('Power / Frequency [ 1 / min^{-1} ]')
    xlabel('Horizontal Scale [km]');

    if axisNum == 1
        hold on
    else
        hold off
        legend('Morning', 'Evening');
        
        ylimits = get(gca, 'ylim');
        line([2600 2600], [ylimits(1) ylimits(2)], 'LineStyle', '--', 'Color', 'k');
        line([11800 11800], [ylimits(1) ylimits(2)], 'LineStyle', '--', 'Color', 'k');
    end
end

end

function h = ciplot(lower,upper,x,colour)
     
% ciplot(lower,upper)       
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.

% Raymond Reynolds 24/11/06

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end

if nargin<4
    colour='b';
end

if nargin<3
    x=1:length(lower);
end

% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end

h = fill([x fliplr(x)],[upper fliplr(lower)],colour);
%camlight; lighting gouraud;
alpha(0.5)

end


function varargout = cmapline(varargin)
% CMAPLINE - Apply a colormap to lines in plot
%
%   CMAPLINE finds all lines in an axis and specifies their
%   colors according to a colormap. Also accepts custom 
%   colormaps in the form of a n x 3 matrix. 
%    
% OPTIONS and SYNTAX
%
%   cmapline - with no inputs, cmapline finds all lines in the 
%   current axis and applies the colormap 'jet'.
% 
%   cmapline('ax',gca,'colormap','hot') - will find all lines
%   in the specified axis (in this case, the current axis)
%   and applies the colormap 'hot'.
%
%   cmapline('lines',handles) - applies colormap values to line
%   objects with specified handles.   
%
%   cmapline('filled') - will fill markers (if included in the 
%   line) with corresponding colormap colors.
%
%   lineh=cmapline - The optional output variable returns the
%   handles to the line objects.
%
%   [lineh, cmap]=cmapline - Two optional outputs returns both the 
%   the handles to the line objects and the applied colormap. 
%
% EXAMPLE 1 - color lines in two subplots according to different colormaps
%  
%   %generate some data
%   x=(0:0.3:2*pi);
%   m=10;
%   exdata=bsxfun(@plus,repmat(10.*sin(x),[m 1]),[1:m]');
%   
%   figure
%   subplot(121);
%   plot(x,exdata,'o-','linewidth',2)
%   cmapline('colormap','jet');
%   set(gca,'color','k')
%   title('jet colormap')
%
%   subplot(122);
%   plot(x,exdata,'o-','linewidth',2)
%   custommap=flipud(hot);
%   cmapline('colormap',custommap,'filled')
%   set(gca,'color','k')
%   title('reverse hot colormap, filled markers')  
%
% EXAMPLE 2 (uses data from example 1) - add a colorbar to your plot
%
%   figure
%   plot(x,exdata,'linewidth',2)
%   [lh,cmap]=cmapline('colormap','jet');
%   colormap(cmap)
%   colorbar
%
% SEE ALSO  colormap 

% Copyright (c) 2010, Andrew Stevens
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

% Andrew Stevens @ USGS, 8/15/2008
% astevens@usgs.gov

%default values
ax=gca;
cmap=@jet;
fillflag=0;
lh=[];

%parse inputs and do some error-checking
if nargin>0
    [m,n]=size(varargin);
    opts={'ax','lines','colormap','filled'};

    for i=1:n;
        indi=strcmpi(varargin{i},opts);
        ind=find(indi==1);
        if isempty(ind)~=1
            switch ind
                case 1
                    %make sure input is an axes handle, sort of
                    ax=varargin{i+1};
                    if ~ishandle(ax)
                        error(['Specified axes',...
                            ' must be a valid axis handle'])
                    end
                case 2
                    lh=varargin{i+1};
                    if ~all(ishandle(lh))
                        error('Invalid line handle')
                    else
                        lh=num2cell(lh);
                    end
                    
                case 3
                    cmap=varargin{i+1};
                    if isa(cmap,'function_handle')
                        cmap= func2str(cmap);
                    end
                    %check size of numeric colormap input
                    if isa(cmap,'numeric')
                        [m,n]=size(cmap);
                        if n~=3
                            error('Custom colormap must have 3 columns.')
                        end
                    end
                case 4
                    fillflag=1;

            end
        else
        end
    end
end

%find lines in axes
if isempty(lh)
    lh=num2cell(findobj(ax,'type','line'));
end

numlines=numel(lh);
if isempty(lh)
    fprintf('No lines present in specified axes.\n')
end

if isa(cmap,'numeric')
    %if needed, interpolate colormap to number of lines
    if numlines~=m
        int=m/numlines;
        ivec=1:m;
        ovec=1:int:1+(numlines-1)*int;

        cmap=num2cell(cmap,1);
        cmap=cellfun(@(x)(interp1(ivec,x,ovec,...
            'linear','extrap')'),cmap,'uni',0);
        colrs=num2cell(cell2mat(cmap),2);
    else
        colrs=num2cell(cmap,2);
    end
else
    %if standard colormap is supplied
    colrs=num2cell(feval(cmap,numlines),2);
end

%apply colors to lines
cellfun(@(x,y)(set(x,'color',y)),lh,colrs);

if strcmpi(get(lh{1},'marker'),'none')~=1 && ...
        fillflag==1;
    cellfun(@(x,y)(set(x,'markerfacecolor',y)),...
        lh,colrs);
end

%output 
if nargout>0
    varargout{1}=cell2mat(lh);
end
if nargout>1
    varargout{2}=cell2mat(colrs);
end

end

