function timeseriesFigHandle = plotTimeseries(firstDatenum, timestamps1min, timestamps1minFixed, timestampsAbsB, timestamps3h, timestamps3hFixed, ae, ap, absB, ...
    averagedDensityNoBg, density3h, morningAeProxy, eveningAeProxy, morningMsisDensity, eveningMsisDensity, morningTimestamps10s, eveningTimestamps10s)
% plotTimeseries(timestamps, ae, averagedDensity)

[timestamps10s, order] = unique([morningTimestamps10s; eveningTimestamps10s]);
msisDensity = [morningMsisDensity; eveningMsisDensity];
predictedDensity = [morningAeProxy; eveningAeProxy];
msisDensity = msisDensity(order);
predictedDensity = predictedDensity(order);
mean(msisDensity)
mean(predictedDensity)

msisDensity = smooth(msisDensity, 7);
msisDensity = msisDensity(ismember(timestamps10s, timestamps1minFixed));
predictedDensity = smooth(predictedDensity, 7);
predictedDensity = predictedDensity(ismember(timestamps10s, timestamps1minFixed));

averagedMsis = removePeriodicBackground(msisDensity, 125, 1, 0);
averagedPrediction = removePeriodicBackground(predictedDensity, 125, 1, 0);

averagedMsis = normalize(averagedMsis, msisDensity);
averagedPrediction = normalize(averagedPrediction, predictedDensity);


timeseriesFigHandle = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1)
secondsInDay = 60 * 60 * 24;
timestampsInDays1min = timestamps1min / secondsInDay + firstDatenum;
timestampsInDays1minFixed = timestamps1minFixed / secondsInDay + firstDatenum;
timestampsInDaysAbsB = timestampsAbsB / secondsInDay + firstDatenum;
%[hAx,~,~] = plotyy(timestampsInDaysAbsB, absB, timestampsInDays1minFixed, averagedDensityNoBg);
[hAx,~,~] = plotyy(timestampsInDays1minFixed, averagedDensityNoBg, timestampsInDaysAbsB, absB);
title('Timeseries of IMF |B| and density')
datetick(hAx(1), 'x', 'dd')
datetick(hAx(2), 'x', 'dd')
ylabel(hAx(1), 'IMF |B| / nT')
ylabel(hAx(2), 'Density')
set(hAx, 'XLim', [min(timestampsInDays1min) max(timestampsInDays1min)]);
hold all;
plot(timestampsInDays1minFixed, averagedMsis);
plot(timestampsInDays1minFixed, averagedPrediction);
hold off;
axis(hAx(1), 'auto y')
grid on

subplot(2,2,3)
timestampsInDays3h = timestamps3h / secondsInDay + firstDatenum;
timestampsInDays3hFixed = timestamps3hFixed / secondsInDay + firstDatenum;
[hAx,~,~] = plotyy(timestampsInDays3h, ap, timestampsInDays3hFixed, density3h);
title('Timeseries of ap and density')
datetick(hAx(1), 'x', 'dd')
datetick(hAx(2), 'x', 'dd')
xlabel('t / days')
ylabel(hAx(1), 'ap')
ylabel(hAx(2), 'previous 3h average density')
set(hAx, 'XLim', [min(timestampsInDays3h) max(timestampsInDays3h)]);
grid on

subplot(2,2,4)
[hAx,~,~] = plotyy(timestampsInDays1min, ae, timestampsInDays1minFixed, averagedDensityNoBg);
title('Timeseries of AE and density')
datetick(hAx(1), 'x', 'dd')
datetick(hAx(2), 'x', 'dd')
xlabel('t / days')
ylabel(hAx(1), 'AE')
ylabel(hAx(2), 'Density')
set(hAx, 'XLim', [min(timestampsInDays1min) max(timestampsInDays1min)]);
grid on

end