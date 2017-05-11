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
conserveInd = Bfactors(:,4) < 550;
objectIDs = objectIDs(conserveInd);

if ~all(ismember(independentID, objectIDs))
    error(['Could not find requested object(s): ', num2str(independentID(~ismember(independentID, objectIDs))),' in Bfactors.dat'])
end

load('ilData.mat','originalRhoStruct');

tleMap = downloadTLEs(objectIDs, beginDate, endDate);

M = floor((endDate-beginDate)/assimilationWindow) + 3;
plotTimes = zeros(M,1);
plotOM = zeros(M,length(independentID),3);
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

rhoStruct = originalRhoStruct;

Tex_series = [];
dT_series = [];
T0_series = [];
obsRank = [];
date = beginDate;
oldTLEs = selectTLEs(tleMap, 'oldest');
k = 1;
assTimes = [];
while date <= endDate
    assimilatableTLEs = findAssimilatableTLEs(tleMap, oldTLEs, date, date + assimilationWindow, intWindow);
    if ~isempty(keys(assimilatableTLEs))
        fprintf('Num obj.: %d\n', length(keys(assimilatableTLEs)))
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
        assTimes = [assTimes; date + assimilationWindow/2];
%         removeInd = rhoStruct.timestamps < date | rhoStruct.timestamps >= date+assimilationWindow;
%         CH_GR = removeDataPoints(rhoStruct, removeInd,false,true,true,true);
%         CH_GR = averageRho(CH_GR, true, 240);
%         consInd = false(size(CH_GR.data)); consInd(1:6:end) = true;
%         CH_GR = removeDataPoints(CH_GR, ~consInd,false,true,true,true);
%         CH_GR.sigma = CH_GR.sigma ./ CH_GR.data;
%         CH_GR.data = log(CH_GR.data);
%         CH_GR = computeVariablesForFit(CH_GR);
%         fprintf('Num CH_GR.: %d\n', length(CH_GR.data))
%               
%         CH_GR.dTCoeff = assimiStruct.dTCoeff;
%         CH_GR.T0Coeff = assimiStruct.T0Coeff;
%         CH_GR.TexCoeff = assimiStruct.TexCoeff;
%         CH_GR.OCoeff = assimiStruct.OCoeff;
%         CH_GR.N2Coeff = assimiStruct.N2Coeff;
%         CH_GR.HeCoeff = assimiStruct.HeCoeff;
%         CH_GR.ArCoeff = assimiStruct.ArCoeff;
%         CH_GR.O2Coeff = assimiStruct.O2Coeff;
%         CH_GR.O_numBiases = assimiStruct.O_numBiases;
%         CH_GR.N2_numBiases = assimiStruct.N2_numBiases;
%         CH_GR.He_numBiases = assimiStruct.He_numBiases;
%         CH_GR.Ar_numBiases = assimiStruct.Ar_numBiases;
%         CH_GR.O2_numBiases = assimiStruct.O2_numBiases;
%         [ensemble] = ...
%             assimilateDataAndUpdateEnsemble(ensemble, @il_model_operator, CH_GR, false, false);
        
        S_orig = S;
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
        [ensemble,~,~,~,~,obsRank_this] = ...
            assimilateDataAndUpdateEnsemble(ensemble, obsOperator, S, false, true);
        if k > 1
            ensemble(1,:) = (TexStd)/ensStd(1) * (ensemble(1,:)-ensMean(1)) + ensMean(1);
        end
        if k > 3
            obsRank = [obsRank; obsRank_this];
        end
        
        % Compute analysis for the independent objects
        S = S_orig;
        if sum(ind) > 0
            analysisTLEs = containers.Map('KeyType', 'double', 'ValueType', 'any');
            for i = 1:length(independentID)
                if sum(S.objectIDs == independentID(i)) > 0
                    analysisTLEs(independentID(i)) = assimilatableTLEs(independentID(i));
                end
            end

            if strcmpi(modelString,'full')
                S = computeBiRhoAndIntTerms(ensemble, modelOperator, oldTLEs, analysisTLEs, 0.5, 100, Ftimes,F,FA,aeInt,assimiStruct,false,[],true);
            else
                S = computeBiRhoAndIntTerms(ensemble, modelOperator, oldTLEs, analysisTLEs, 0.5, 100);
            end
            
            ind = ismember(S.objectIDs,independentID);
            OM_DA = exp(S.rhoObs(ind))./exp(mean(S.rhoModel_DA(ind,:),2));
            this_computed_satell = S.objectIDs(ind);
            for j = 1:length(independentID)
                objInd = find(this_computed_satell == independentID(j));
                if ~isempty(objInd)
                    plotOM(k,j,3) = OM_DA(objInd);
                end
            end
        end
        
        Tex_series = vertcat(Tex_series, [quantile(ensemble(1,:),0.1), mean(ensemble(1,:)), quantile(ensemble(1,:),0.9)]);
        T0_series = vertcat(T0_series, [quantile(ensemble(2,:),0.1), mean(ensemble(2,:)), quantile(ensemble(2,:),0.9)]);
        dT_series = vertcat(dT_series, [quantile(ensemble(3,:),0.1), mean(ensemble(3,:)), quantile(ensemble(3,:),0.9)]);
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
subplot(3,1,1)
%plot(assTimes,Tex_series(:,1),'r--','linewidth',2.0)
%hold on
plot(assTimes,Tex_series(:,2),'b','linewidth',2.0)
%plot(assTimes,Tex_series(:,3),'r--','linewidth',2.0)
set(gca,'fontsize',15)
[mean(Tex_series(:,2)),std(Tex_series(:,2))]
title('\Delta T_{ex}','fontsize',15)
xlim([beginDate, endDate]);
set(gca,'xticklabel',[])

subplot(3,1,2)
%plot(assTimes,T0_series(:,1),'r--','linewidth',2.0)
%hold on
plot(assTimes,T0_series(:,2),'b','linewidth',2.0)
%plot(assTimes,T0_series(:,3),'r--','linewidth',2.0)
set(gca,'fontsize',15)
[mean(T0_series(:,2)),std(T0_series(:,2))]
title('\Delta T_{0}','fontsize',15)
xlim([beginDate, endDate]);
set(gca,'xticklabel',[])

subplot(3,1,3)
%plot(assTimes,dT_series(:,1),'r--','linewidth',2.0)
%hold on
plot(assTimes,dT_series(:,2),'b','linewidth',2.0)
%plot(assTimes,dT_series(:,3),'r--','linewidth',2.0)
set(gca,'fontsize',15)
[mean(dT_series(:,2)),std(dT_series(:,2))]
title('\Delta T','fontsize',15)
datetick('x')
xlim([beginDate, endDate]);



figure;
%hold all;
ylim([min(plotOM(:)), max(plotOM(:))])
hAx = zeros(size(independentID));
N_ind = length(independentID);
for i = 1:length(independentID)
    ind = plotOM(:,i,1) > 0;
    pt = plotTimes(ind);
    pOM_anal = plotOM(ind,i,3); % analysis = 3, bg. = 1
    pOM_IL = plotOM(ind,i,2);
    pOM_bg = plotOM(ind,i,1); % analysis = 3, bg. = 1
    subplot(N_ind,1,i);
    if isempty(pOM_anal); continue; end
    [hAx(i)] = plot(pt, pOM_anal,'linewidth', 2.0);
    hold all;
    h_bg = plot(pt, pOM_bg,'--','linewidth', 2.0);
    set(h_bg,'color',get(hAx(i),'color'));
    h_IL = plot(pt, pOM_IL,':','linewidth', 2.0);
    set(h_IL,'color',get(hAx(i),'color'));
    hold off;
    %set(gca,'xticklabel',[])
    datetick('x');
    xlim([beginDate, endDate]);
    
    set(gca,'fontsize',15)
    %grid on
    title([num2str(independentID(i))],'fontsize',15);
    %ylabel('\rho_{hav.} / \rho_{malli}','fontsize',15)
    
    k = 3:length(pOM_IL);
    [analBetter, p_analWorse] = ttest((pOM_IL(k)-1).^2, (pOM_anal(k)-1).^2,'tail','right');
    [bgBetter, p_bgWorse] = ttest((pOM_IL(k)-1).^2, (pOM_bg(k)-1).^2,'tail','right');
    fprintf('%d: anal. better: %d (%f), bg. better: %d (%f)\n', independentID(i), analBetter, p_analWorse, bgBetter, p_bgWorse);
    
    rms_IL = rms(pOM_IL(k)-1);
    rms_anal = rms(pOM_anal(k)-1);
    rms_bg = rms(pOM_bg(k)-1);
    fprintf('%d: IL_rms: %f, bg_rms: %f, anal_rms: %f\n',independentID(i),rms_IL,rms_bg,rms_anal)
end
%title(['\rho_{hav.} / \rho_{malli}'],'fontsize',15)
%legend(hAx(hAx~=0),strsplit(num2str(independentID(hAx~=0))));


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