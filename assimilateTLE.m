function [] = assimilateTLE(beginDateStr, endDateStr, modelString, ...
    assimilationWindow, intWindow, independentID, TexStd)
% [assimilationWindow] = days, [intWindow] = days

if ~verLessThan('matlab','8.6')
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        parpool();
    end
else
    if matlabpool('size') == 0
        matlabpool('open',feature('numCores'))
    end
end

if iscolumn(independentID)
    independentID = independentID';
end

independentID = unique(independentID);

beginDate = datenum(beginDateStr);
endDate = datenum(endDateStr);

load Bfactors.dat
objectIDs = Bfactors(:,1);

if ~all(ismember(independentID, objectIDs))
    error(['Could not find requested object(s): ', num2str(independentID(~ismember(independentID, objectIDs))),' in Bfactors.dat'])
end

tleMap = downloadTLEs(objectIDs, beginDate, endDate);

M = floor((endDate-beginDate)/assimilationWindow) + 3;
plotTimes = zeros(M,1);
plotOM = zeros(M,length(independentID),2);
plotMap = containers.Map('keytype','double','valuetype','double');
for i = 1:length(independentID);
    plotMap(independentID(i)) = i;
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
    load('ilData.mat','rhoStruct')
    [Ftimes,order] = unique(rhoStruct.timestamps);
    F = rhoStruct.F(order);
    FA = rhoStruct.FA(order);
    aeInt = rhoStruct.aeInt(order,:);
    load('assimiStruct.mat')
else
    modelOperator = @dummyThermosphere;
end

obsRank = [];
date = beginDate;
oldTLEs = selectTLEs(tleMap, 'oldest');
k = 1;
while date <= endDate
    assimilatableTLEs = findAssimilatableTLEs(tleMap, oldTLEs, date, date + assimilationWindow, intWindow);
    if ~isempty(keys(assimilatableTLEs))
        if strcmpi(modelString,'full')
            S = computeBiRhoAndIntTerms(ensemble, modelOperator, oldTLEs, assimilatableTLEs, 0.5, 100, Ftimes,F,FA,aeInt,assimiStruct,false,[],true);
        else
            S = computeBiRhoAndIntTerms(ensemble, modelOperator, oldTLEs, assimilatableTLEs, 0.5, 100);
        end
        
        ind = ismember(S.objectIDs,independentID);
        if sum(ind) > 0
            OM_DA = exp(S.rhoObs(ind))./exp(mean(S.rhoModel_DA(ind,:),2));
            OM_IL = exp(S.rhoObs(ind))./exp(S.rhoModel_IL(ind));
            plotTimes(k) = date + assimilationWindow/2;
            this_computed_satell = S.objectIDs(ind);
            for j = 1:length(independentID)
                objInd = find(this_computed_satell == independentID(j));
                if ~isempty(objInd)
                    plotOM(k,j,1) = OM_DA(objInd);
                    plotOM(k,j,2) = OM_IL(objInd);
                end
            end
        end
        
%         conserveInd = ~ismember(S.objectIDs, independentID);
%         S.Bi = S.Bi(conserveInd,:);
%         S.Bratio = S.Bratio(conserveInd,:);
%         S.objectIDs = S.objectIDs(conserveInd,:);
%         S.rhoObs = S.rhoObs(conserveInd,:);
%         S.rhoModel_DA = S.rhoModel_DA(conserveInd,:);
%         S.rhoModel_IL = S.rhoModel_IL(conserveInd,:);
%         S.sig_rho = S.sig_rho(conserveInd,:);
        S.sigma = S.sig_rho;
        S.data = S.rhoObs;

        ensMean = mean(ensemble, 2);
        ensStd = std(ensemble, 0, 2);     
        if k > 1
            ensemble(1,:) = (TexStd)/ensStd(1) * (ensemble(1,:)-ensMean(1)) + ensMean(1);
        end
        [ensemble,~,~,~,~,obsRank_this] = ...
            assimilateDataAndUpdateEnsemble(ensemble, obsOperator, S, false, true);
        if k > 3
            obsRank = [obsRank; obsRank_this];
        end
        
    end
    
    assimilatedObj = keys(assimilatableTLEs);
    for i = 1:length(assimilatedObj)
        oldTLEs(assimilatedObj{i}) = assimilatableTLEs(assimilatedObj{i});
    end
    
    date = date + assimilationWindow;
    k = k + 1;
    p.progress;
end
p.stop;
ind = plotTimes > 0;
plotTimes = plotTimes(ind);
plotOM = plotOM(ind,:,:);

figure;
hold all;
ylim([min(plotOM(:)), max(plotOM(:))])
hAx = zeros(size(independentID));
for i = 1:length(independentID)
    ind = plotOM(:,i,1) > 0;
    pt = plotTimes(ind);
    pOM_DA = plotOM(ind,i,1);
    pOM_IL = plotOM(ind,i,2);
    if isempty(pOM_DA); continue; end
    [hAx(i)] = plot(pt, pOM_DA,'linewidth', 2.0);
    h_IL = plot(pt, pOM_IL,'--','linewidth', 2.0);
    set(h_IL,'color',get(hAx(i),'color'));
end
title(['\rho_{hav.} / \rho_{malli}',num2str(TexStd)],'fontsize',15)
legend(hAx(hAx~=0),strsplit(num2str(independentID(hAx~=0))));
datetick('x')
set(gca,'fontsize',15)
grid on

figure;
hist(obsRank, 20);
title(['Sijoituslukujen jakauma',num2str(TexStd)],'fontsize',15);

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