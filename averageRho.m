function rhoStruct = averageRho(rhoStruct)

averTime = 120; % sec

dt = diff(rhoStruct.timestamps); dt = [dt; dt(end)] * 86400;
averInterval = abs(averTime ./ dt); 
removeInd = true(length(rhoStruct.data), 1);
disconts = find(averInterval < 1); disconts = [0; disconts];

targetCount = round(length(disconts) / 100);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Averaging rho, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

for i = 2:length(disconts)
    ind = disconts(i-1)+1 : disconts(i);
    interval = max(max(round(averInterval(ind))), 1);
    if interval > 1 && mod(interval, 2) == 0; 
        interval = interval - 1; 
    end
    
    rhoStruct.latitude(ind) = smooth(rhoStruct.latitude(ind), interval);
    rhoStruct.longitude(ind) = smooth(rhoStruct.longitude(ind), interval);
    rhoStruct.solarTime(ind) = smooth(rhoStruct.solarTime(ind), interval);
    rhoStruct.altitude(ind) = smooth(rhoStruct.altitude(ind), interval);
    rhoStruct.F(ind) = smooth(rhoStruct.F(ind), interval);
    rhoStruct.FA(ind) = smooth(rhoStruct.FA(ind), interval);
    rhoStruct.apNow(ind) = smooth(rhoStruct.apNow(ind), interval);
    rhoStruct.ap3h(ind) = smooth(rhoStruct.ap3h(ind), interval);
    rhoStruct.ap6h(ind) = smooth(rhoStruct.ap6h(ind), interval);
    rhoStruct.ap9h(ind) = smooth(rhoStruct.ap9h(ind), interval);
    rhoStruct.ap12To33h(ind) = smooth(rhoStruct.ap12To33h(ind), interval);
    rhoStruct.ap36To57h(ind) = smooth(rhoStruct.ap36To57h(ind), interval);
    rhoStruct.Ap(ind) = smooth(rhoStruct.Ap(ind), interval);
    rhoStruct.data(ind) = smooth(rhoStruct.data(ind), interval);
    rhoStruct.timestamps(ind) = smooth(rhoStruct.timestamps(ind), interval);
    
    conserve = min((interval+1)/2, length(ind)) : interval : length(ind);
    conserve = conserve + ind(1) - 1;
    removeInd(conserve) = false;
    
    if mod(i, 100) == 0
        p.progress;
    end
end
p.stop;

rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, false, false);

end