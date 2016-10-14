function [] = assimilateTLE(beginDateStr, endDateStr, assimilationWindow, intWindow, plotID)
% [assimilationWindow] = days, [intWindow] = days

if iscolumn(plotID)
    plotID = plotID';
end

beginDate = datenum(beginDateStr);
endDate = datenum(endDateStr);

load Bfactors.dat
objectIDs = Bfactors(:,1);

if ~all(ismember(plotID, objectIDs))
    error(['Could not find requested object(s): ', plotID(~ismember(plotID, objectIDs)),' in Bfactors.dat'])
end

tleMap = downloadTLEs(objectIDs, beginDate, endDate);

M = floor((endDate-beginDate)/assimilationWindow) + 3;
plotTimes = zeros(M,1);
plotOM = zeros(M,length(plotID));

targetCount = round(M);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running TLE asimilation, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

date = beginDate;
oldTLEs = selectTLEs(tleMap, 'oldest');
k = 1;
while date <= endDate
    assimilatableTLEs = findAssimilatableTLEs(tleMap, oldTLEs, date, date + assimilationWindow, intWindow);
    if ~isempty(keys(assimilatableTLEs))
        S = computeBiRhoAndIntTerms(zeros(11,1), @dummyThermosphere, oldTLEs, assimilatableTLEs, 0.5, 100);
        ind = S.objectIDs == plotID;
        OM = S.rhoObs(ind)./S.rhoModel(ind);
        plotTimes(k) = date + assimilationWindow/2;
        plotOM(k,:) = OM;
    end
    date = date + assimilationWindow;
    k = k + 1;
    p.progress;
end
p.stop;
ind = plotTimes > 0;
plotTimes = plotTimes(ind);
plotOM = plotOM(ind,:);

figure;
plot(repmat(plotTimes,1,length(plotID)), plotOM,'linewidth', 2.0)
legend(strsplit(num2str(plotID)));
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