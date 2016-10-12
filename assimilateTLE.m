function [] = assimilateTLE(beginDateStr, endDateStr, assimilationWindow, intWindow, plotID)
% [assimilationWindow] = days, [intWindow] = days

beginDate = datenum(beginDateStr);
endDate = datenum(endDateStr);

load Bfactors.dat
objectIDs = Bfactors(:,1);

tleMap = downloadTLEs(objectIDs, beginDate, endDate);

M = floor((beginDate-endDate)/assimilationWindow) + 3;
plotTimes = zeros(M,1);
plotOM = zeros(M,1);

date = beginDate;
oldTLEs = selectTLEs(tleMap, 'oldest');
k = 1;
while date <= endDate
    assimilatableTLEs = findAssimilatableTLEs(tleMap, oldTLEs, date, date + assimilationWindow, intWindow);
    S = computeBiRhoAndIntTerms(zeros(11,1), @dummyThermosphere, oldTLEs, assimilatableTLEs, 0.5, 100);
    ind = S.objectIDs == plotID;
    OM = S.rhoObs(ind)./S.rhoModel(ind);
    plotTimes(k) = date + assimilationWindow/2;
    plotOM(k) = OM;
    date = date + assimilationWindow;
    k = k + 1;
end
ind = plotTimes > 0;
plotTimes = plotTimes(ind);
plotOM = plotOM(ind);

figure;
plot(plotTimes, plotOM,'linewidth', 2.0)
datetick('x')


end

function assimilatableTLEs = findAssimilatableTLEs(tleMap, oldTLEs, beginDate, endDate, intWindow)

[year,month,day,hh,mm,ss] = datevec(beginDate);
beginJul = jday(year,month,day,hh,mm,ss);
[year,month,day,hh,mm,ss] = datevec(endDate);
endJul = jday(year,month,day,hh,mm,ss);

TLEsInWindow = containers.Map('KeyType', 'double', 'ValueType', 'any');
obj = keys(tleMap);
for i = 1:length(obj)
    tles = tleMap(obj{i});
    ind = [];
    for j = 1:length(tles.sgp4info)
        sgp4tle = tles.sgp4info(j);
        if beginJul <= sgp4tle.jdsatepoch && sgp4tle.jdsatepoch <= endJul
            ind = [ind;j];
        end
    end
    if ~isempty(ind)
        TLEsInWindow(obj{i}) = struct('Btrue',tles.Btrue,'sig_Btrue',tles.sig_Btrue,...
                                      'sgp4info', tles.sgp4info(ind));
    end
end

assimilatableTLEs = containers.Map('KeyType', 'double', 'ValueType', 'any');
obj = keys(TLEsInWindow);
for i = 1:length(obj)
    tles = TLEsInWindow(obj{i});
    oldTime = oldTLEs(obj{i}).sgp4info.jdsatepoch;
    ind = [];
    for j = 1:length(tles.sgp4info)
        sgp4tle = tles.sgp4info(j);
        if sgp4tle.jdsatepoch >= oldTime + intWindow
            ind = j;
            break;
        end
    end
    if ~isempty(ind)
        assimilatableTLEs(obj{i}) = struct('Btrue',tles.Btrue,'sig_Btrue',tles.sig_Btrue,...
                                      'sgp4info', tles.sgp4info(ind));
    end
end

end