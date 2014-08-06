function compareGoceDensityToMsis(goceDensity, msisDensity, ae, timestampsAeDatenum, timestampsDensityDatenum, results)
%

ratio = goceDensity ./ msisDensity;
ratioTrend = removePeriodicBackground(ratio, 125, 6, 0);

figure;
[hAx,hLine1,hLine2] = plotyy(timestampsDensityDatenum, ratio, timestampsAeDatenum, ae);
set(hLine1, 'LineStyle', 'none', 'Marker', '.')
% For some reason Matlab can't handle axis limits properly, so they need to
% be set manually.
ratioYmin = 0.1 * floor(min(ratio) / 0.1);
ratioYmax = 0.1 * ceil(max(ratio) / 0.1);
set(hAx(1), 'YLim', [ratioYmin ratioYmax], 'YTick', ratioYmin:0.1:ratioYmax);
set(hAx, 'XLim', [min(timestampsDensityDatenum) max(timestampsDensityDatenum)], 'XMinorTick', 'on'); 
aeYmax = 500 * ceil(max(ae) / 500);
set(hAx(2),'YLim', [0 aeYmax], 'YTick', 0:250:aeYmax);
datetick(hAx(1), 'x', 'yyyy-mm', 'keepticks', 'keeplimits')
datetick(hAx(2), 'x', 'yyyy-mm', 'keepticks', 'keeplimits')
rotateticklabel(hAx(1), 50);
rotateticklabel(hAx(2), 50);
hold on;
plot(timestampsDensityDatenum, ratioTrend, 'r-');
hold off;
%plotyy(timestampsInDays, correctedDensity, timestampsInDays, msisDensity270km);
ylabel(hAx(1), 'GOCE / NRLMSISE00 density')
ylabel(hAx(2), 'AE')
title('Ratio of GOCE NON-normalized densities to NRLMSISE00 prediction')

meanRatio = mean(ratio);
errRatio = std(ratio);
meanReciprocal = mean(1 ./ ratio);
errReciprocal = std(1 ./ ratio);
textString1 = ['Mean GOCE / NRLMSISE00 density: ', num2str(meanRatio, '%.2f'), ' ± ', num2str(errRatio, '%.2f')];
textString2 = ['Mean NRLMSISE00 / GOCE density: ', num2str(meanReciprocal, '%.2f'), ' ± ', num2str(errReciprocal, '%.2f')];
     
ylimits = get(gca, 'ylim');
xlimits = get(gca, 'xlim');
textYLocation = ylimits(2) - 0.05 * (ylimits(2) - ylimits(1));
textXLocation = xlimits(2) - 0.05 * (xlimits(2) - xlimits(1));
text(textXLocation, textYLocation, textString1, 'FontSize', 9, ...
    'VerticalAlignment','top', 'HorizontalAlignment','right');

textYLocation = ylimits(2) - 0.1 * (ylimits(2) - ylimits(1));
textXLocation = xlimits(2) - 0.05 * (xlimits(2) - xlimits(1));
text(textXLocation, textYLocation, textString2, 'FontSize', 9, ...
    'VerticalAlignment','top', 'HorizontalAlignment','right');

plotCorrelation(msisDensity, goceDensity, 'NRLMSISE00 density', 'GOCE measured density', 1, results);

figure;
ratio = sort(ratio);
minX = floor(min(ratio) * 100) / 100;
maxX = ceil(max(ratio) * 100) / 100;

[numOfElementsInBin, x] = hist(ratio, minX:0.01:maxX);

pdNormal = fitdist(ratio, 'Normal');
pdLognormal = fitdist(ratio, 'Lognormal');
pdTlocationScale = fitdist(ratio, 'tlocationscale');

pdfNormal = pdf(pdNormal, ratio);
pdfLognormal = pdf(pdLognormal, ratio);
pdfTlocationScale = pdf(pdTlocationScale, ratio);

dx = diff(x(1:2));
bar(x,numOfElementsInBin / sum(numOfElementsInBin * dx))
title('Goce/Msis histogram with pdf fits')
xlabel('Goce/NRLMSISE00 density ratio')
ylabel('Percentage of data');
xlim([0.25 1.25])

hold all
normalHandle = plot(ratio, pdfNormal, 'r', 'Linewidth', 2);
lognHandle = plot(ratio, pdfLognormal, 'y', 'Linewidth', 2);
tlocHandle = plot(ratio, pdfTlocationScale, 'g', 'Linewidth', 2);
hold off

h = [normalHandle, lognHandle, tlocHandle];
legend(h, 'Normal', 'Lognormal', 't Location-Scale'); 

figure;
xValuesLognormal = random(pdLognormal, length(ratio), 1);
xValuesTlocScale = random(pdTlocationScale, length(ratio), 1);
xValuesTlocScale(xValuesTlocScale <= 0) = [];

tlocHandle = normplot(xValuesTlocScale);
hold all;
lognHandle = normplot(xValuesLognormal);
ratioHandle = normplot(ratio);

set(tlocHandle(1), 'Color', 'g')
set(lognHandle(1), 'Color', 'y')
set(tlocHandle(2), 'LineStyle', 'none')
set(lognHandle(2), 'LineStyle', 'none')
set(tlocHandle(3), 'LineStyle', 'none')
set(lognHandle(3), 'LineStyle', 'none')

h = [ratioHandle(3), lognHandle(1), tlocHandle(1), ratioHandle(1)];
legend(h, 'Normal', 'Lognormal', 't Location-Scale', 'Ratio', 'Location', 'NorthWest'); 

end
