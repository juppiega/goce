function timeseriesFigHandle = plotTimeseries(firstDatenum, timestamps1min, timestamps1minFixed, timestampsAbsB, timestamps3h, timestamps3hFixed, ae, ap, absB, ...
    averagedDensityNoBg, density3h)
% plotTimeseries(timestamps, ae, averagedDensity)

timeseriesFigHandle = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1)
secondsInDay = 60 * 60 * 24;
timestampsInDays1min = timestamps1min / secondsInDay + firstDatenum;
timestampsInDays1minFixed = timestamps1minFixed / secondsInDay + firstDatenum;
timestampsInDaysAbsB = timestampsAbsB / secondsInDay + firstDatenum;
[hAx,~,~] = plotyy(timestampsInDaysAbsB, absB, timestampsInDays1minFixed, averagedDensityNoBg);
title('Timeseries of IMF |B| and density')
datetick(hAx(1), 'x', 'dd')
datetick(hAx(2), 'x', 'dd')
ylabel(hAx(1), 'IMF |B| / nT')
ylabel(hAx(2), 'Density')
set(hAx, 'XLim', [min(timestampsInDays1min) max(timestampsInDays1min)]);
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