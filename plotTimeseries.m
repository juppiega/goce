function timeseriesFigHandle = plotTimeseries(firstDatenum, timestamps1min, timestamps1minFixed, timestampsAbsB, timestamps3h, timestamps3hFixed, ae, ap, absB, ...
    averagedDensityNoBg, density3h, morningAeProxy, eveningAeProxy, morningMsisDensity, morningJbDensity, eveningMsisDensity, eveningJbDensity, morningTimestamps10s, eveningTimestamps10s)
% plotTimeseries(timestamps, ae, averagedDensity)

[timestamps10s, order] = unique([morningTimestamps10s; eveningTimestamps10s]);
msisDensity = [morningMsisDensity; eveningMsisDensity];
jb2008Density = [morningJbDensity; eveningJbDensity];
predictedDensity = [morningAeProxy; eveningAeProxy];
msisDensity = msisDensity(order);
jb2008Density = jb2008Density(order);
predictedDensity = predictedDensity(order);

msisDensity = smooth(msisDensity, 7);
msisDensity = msisDensity(ismember(timestamps10s, timestamps1minFixed));
jb2008Density = smooth(jb2008Density, 7);
jb2008Density = jb2008Density(ismember(timestamps10s, timestamps1minFixed));
predictedDensity = smooth(predictedDensity, 7);
predictedDensity = predictedDensity(ismember(timestamps10s, timestamps1minFixed));

averagedMsis = removePeriodicBackground(msisDensity, 125, 1, 0);
averagedJb = removePeriodicBackground(jb2008Density, 125, 1, 0);
averagedPrediction = removePeriodicBackground(predictedDensity, 125, 1, 0);

averagedMsis = normalize(averagedMsis, msisDensity);
averagedJb = normalize(averagedJb, jb2008Density);
averagedPrediction = normalize(averagedPrediction, predictedDensity);

figure;
secondsInDay = 60 * 60 * 24;
timestampsInDays1minFixed = timestamps1minFixed / secondsInDay + firstDatenum;
plot(timestampsInDays1minFixed, averagedDensityNoBg * 1.22);
hold all;
plot(timestampsInDays1minFixed, averagedMsis);
plot(timestampsInDays1minFixed, averagedJb);
plot(timestampsInDays1minFixed, averagedPrediction);
hold off;
legend('GOCE', 'MSIS', 'JB2008', 'AE Model')
datetick('x', 'dd')
title('Model low pass filtered density comparison')
xlabel(['Date on ', datestr(timestampsInDays1minFixed(1), 'mmm yyyy')])
ylabel('Density [10^{-12} kg/m^3]')

timeseriesFigHandle = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1)
timestampsInDays1min = timestamps1min / secondsInDay + firstDatenum;
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