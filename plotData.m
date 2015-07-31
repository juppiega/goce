function plotData(beginDay, endDay)

load('goceVariables.mat', 'timestampsAeDatenum', 'ae', 'timestampsEpsilonDatenum', 'akasofuEpsilon')
beginDay = datenum(beginDay);
endDay = datenum(endDay);

aeInd = beginDay <= timestampsAeDatenum & timestampsAeDatenum <= endDay;
epsInd = beginDay <= timestampsEpsilonDatenum & timestampsEpsilonDatenum <= endDay;
tAe = timestampsAeDatenum(aeInd);
tEps = timestampsEpsilonDatenum(epsInd);

epsilon = (9*6371E3).^2 * akasofuEpsilon(epsInd) * 1E-2;
ae = ae(aeInd);

figure;
subplot(2,1,1);
plot(tAe, ae, 'linewidth', 1.5, 'color', 'k')
set(gca, 'fontsize', 13);
title('AE');
xlabel('Date')
xlim([min(tAe), max(tAe)])
datetick

subplot(2,1,2)
plot(tEps, epsilon, 'linewidth', 1.5, 'color', 'k')
set(gca, 'fontsize', 13);
title('Akasofu epsilon')
xlabel('Date')
xlim([min(tEps), max(tEps)])
datetick

end