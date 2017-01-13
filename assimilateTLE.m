function [] = assimilateTLE(beginDateStr, endDateStr, modelString, ...
    assimilationWindow, intWindow, plotID, TexStd)
% [assimilationWindow] = days, [intWindow] = days

if iscolumn(plotID)
    plotID = plotID';
end

plotID = unique(plotID);

beginDate = datenum(beginDateStr);
endDate = datenum(endDateStr);

load Bfactors.dat
objectIDs = Bfactors(:,1);

if ~all(ismember(plotID, objectIDs))
    error(['Could not find requested object(s): ', num2str(plotID(~ismember(plotID, objectIDs))),' in Bfactors.dat'])
end

tleMap = downloadTLEs(objectIDs, beginDate, endDate);

M = floor((endDate-beginDate)/assimilationWindow) + 3;
plotTimes = zeros(M,1);
plotOM = zeros(M,length(plotID));
plotMap = containers.Map('keytype','double','valuetype','double');
for i = 1:length(plotID);
    plotMap(plotID(i)) = i;
end
targetCount = round(M);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running TLE asimilation, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
               
ensembleSize = 100;
ensemble = createInitialEnsemble(modelString, ensembleSize);
obsOperator = @tle_model_operator;
if strcmpi(modelString,'full')
    modelOperator = @il_model_operator;
else
    modelOperator = @dummyThermosphere;
end

date = beginDate;
oldTLEs = selectTLEs(tleMap, 'oldest');
k = 1;
while date <= endDate
    assimilatableTLEs = findAssimilatableTLEs(tleMap, oldTLEs, date, date + assimilationWindow, intWindow);
    if ~isempty(keys(assimilatableTLEs))
        S = computeBiRhoAndIntTerms(ensemble, modelOperator, oldTLEs, assimilatableTLEs, 0.5, 100);
        S.sigma = S.sig_rho;
        S.data = S.rhoObs;

        ensMean = mean(ensemble, 2);
        ensStd = std(ensemble, 0, 2);     
        if k > 1
            ensemble(3,:) = (TexStd)/ensStd(3) * (ensemble(3,:)-ensMean(3)) + ensMean(3);
        end
        [ensemble] = assimilateDataAndUpdateEnsemble(ensemble, obsOperator, S, false);
        
        ind = ismember(S.objectIDs,plotID);
        OM = S.rhoObs(ind)./mean(S.rhoModel(ind,:),2);
        plotTimes(k) = date + assimilationWindow/2;
        for j = 1:length(plotID)
            objInd = find(S.objectIDs == plotID(j));
            if ~isempty(objInd)
                plotOM(k,j) = OM(objInd);
            end
        end
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
hold all;
ylim([min(plotOM(:)), max(plotOM(:))])
for i = 1:length(plotID)
    ind = plotOM(:,i) > 0;
    pt = plotTimes(ind);
    pOM = plotOM(ind,i);
    plot(pt, pOM,'linewidth', 2.0)
end
title('\rho_{obs} / \rho_{model}','fontsize',15)
legend(strsplit(num2str(plotID)));
datetick('x')
set(gca,'fontsize',15)

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