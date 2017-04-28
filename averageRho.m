function [S, removeInd] = averageRho(S, isRho, averTime)
% [averTime] = sec

if nargin <= 2 || averTime <= 0
    averTime = 120; % sec
end

dt = diff(S.timestamps); dt = [dt; dt(end)] * 86400;
averInterval = abs(averTime ./ dt); 
removeInd = true(length(S.data), 1);
disconts = find(averInterval < 1); disconts = [0; disconts; length(averInterval)];

targetCount = round(length(disconts) / 100);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Averaging rho, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

for i = 2:length(disconts)
    ind = disconts(i-1)+1 : disconts(i);
    if isempty(ind)
        continue
    end
    interval = max(mode(round(averInterval(ind))), 1);
    if interval > 1 && mod(interval, 2) == 0; 
        interval = interval - 1; 
    end
    
    S.latitude(ind) = smooth(S.latitude(ind), interval);
    S.longitude(ind) = smooth(S.longitude(ind), interval);
    S.solarTime(ind) = smooth(S.solarTime(ind), interval);
    S.altitude(ind) = smooth(S.altitude(ind), interval);
    S.F(ind) = smooth(S.F(ind), interval);
    S.FA(ind) = smooth(S.FA(ind), interval);
    S.apNow(ind) = smooth(S.apNow(ind), interval);
    S.ap3h(ind) = smooth(S.ap3h(ind), interval);
    S.ap6h(ind) = smooth(S.ap6h(ind), interval);
    S.ap9h(ind) = smooth(S.ap9h(ind), interval);
    S.ap12To33h(ind) = smooth(S.ap12To33h(ind), interval);
    S.ap36To57h(ind) = smooth(S.ap36To57h(ind), interval);
    S.Ap(ind) = smooth(S.Ap(ind), interval);
    S.data(ind) = smooth(S.data(ind), interval);
    S.timestamps(ind) = smooth(S.timestamps(ind), interval);
    
    conserve = min((interval+1)/2, length(ind)) : interval : length(ind);
    conserve = conserve + ind(1) - 1;
    removeInd(conserve) = false;
    
    if mod(i, 100) == 0
        p.progress;
    end
end
p.stop;

if isRho
    S = removeDataPoints(S, removeInd, false, true, false, false);
else
    S = removeDataPoints(S, removeInd);
end

end