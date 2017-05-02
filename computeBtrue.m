function computeBtrue(beginDateStr, endDateStr, excludeBegin, excludeEnd, windowLen, plotIDs)
% [windowLen] = days

if ~verLessThan('matlab','8.6')
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        parpool();
    end
end

modelString = 'full';
beginDate = datenum(beginDateStr);
endDate = datenum(endDateStr);
excludeBegin = datenum(excludeBegin);
excludeEnd = datenum(excludeEnd);

IDs = load('assimilationSatellites.dat');
if isrow(IDs)
    IDs = IDs';
end
tleMap = downloadTLEs(IDs, beginDate, endDate, true); % Tahan IDss?

assimilationWindow = windowLen;
intWindow = assimilationWindow;
M = floor((endDate-beginDate)/assimilationWindow) + 3;
computeTimes = zeros(M,1);
Bi = zeros(M,length(IDs));
assAlt = zeros(M,length(IDs));
targetCount = round(M);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running TLE Computations, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
               
obsOperator = @tle_model_operator;
if strcmpi(modelString,'full')
    modelOperator = @il_model_operator;
else
    modelOperator = @dummyThermosphere;
end

load('ilData.mat','rhoStruct')
[Ftimes,order] = unique(rhoStruct.timestamps);
F = rhoStruct.F(order);
FA = rhoStruct.FA(order);
aeInt = rhoStruct.aeInt(order,:);
load('assimiStruct.mat')

ensemble = zeros(size(createInitialEnsemble('full', 1)));

date = beginDate;
oldTLEs = selectTLEs(tleMap, 'oldest');
k = 1;
while date <= endDate
    assimilatableTLEs = findAssimilatableTLEs(tleMap, oldTLEs, date, date + assimilationWindow, intWindow);
    if ~isempty(keys(assimilatableTLEs))
        % KORJAA JA TARKISTA:
        S = computeBiRhoAndIntTerms(ensemble, modelOperator, oldTLEs, assimilatableTLEs, 0.5, 1E30,...
            Ftimes, F, FA, aeInt, assimiStruct, false, [], false);
        
        computeTimes(k) = date + assimilationWindow/2;
        for j = 1:length(IDs)
            objInd = find(S.objectIDs == IDs(j));
            if ~isempty(objInd)
                Bi(k,j) = S.Bi(objInd);
                assAlt(k,j) = S.assocAlt(objInd);
            end
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
ind = computeTimes > 0;
computeTimes = computeTimes(ind);
Bi = Bi(ind,:);
assAlt = assAlt(ind,:);

figure; % TODO: Piirra exludejen valiin pystyviivat
hAx = zeros(size(plotIDs));
for i = 1:length(plotIDs)
    k = find(IDs == plotIDs(i));
    ind = Bi(:,k) > 0;
    pt = computeTimes(ind);
    p_Bi = Bi(ind,k);
    if isempty(p_Bi); continue; end
    [hAx(i)] = plot(pt, p_Bi,'linewidth', 2.0);
end
title('B_i','fontsize',15)
legend(hAx(hAx~=0),strsplit(num2str(plotIDs(hAx~=0))));
datetick('x')
ylabel('B_i [m^2 / kg]','fontsize',15)
set(gca,'fontsize',15)
grid on
line([excludeBegin, excludeBegin], get(gca,'ylim'));
line([excludeEnd, excludeEnd], get(gca,'ylim'));

figure;
hAx = zeros(size(plotIDs));
for i = 1:length(plotIDs)
    k = find(IDs == plotIDs(i));
    ind = assAlt(:,k) > 0;
    pt = computeTimes(ind);
    p_assAlt = assAlt(ind,k);
    if isempty(p_assAlt); continue; end
    [hAx(i)] = plot(pt, p_assAlt,'linewidth', 2.0);
end
title('Tiheyskorkeus','fontsize',15)
legend(hAx(hAx~=0),strsplit(num2str(plotIDs(hAx~=0))));
datetick('x')
ylabel('Korkeus [km]','fontsize',15)
set(gca,'fontsize',15)
grid on
line([excludeBegin, excludeBegin], get(gca,'ylim'));
line([excludeEnd, excludeEnd], get(gca,'ylim'));

Btrue = zeros(size(IDs));
Btrue_sig = zeros(size(IDs));
assAltAver = zeros(size(IDs));
for i = 1:length(IDs)
    ind = Bi(:,i) > 0;
    ind(excludeBegin <= computeTimes & computeTimes <= excludeEnd) = false;
    Bi_this = Bi(ind,i);
    assAlt_this = assAlt(ind,i);
    
    Btrue(i) = mean(Bi_this);
    Btrue_sig(i) = std(Bi_this)/sqrt(sum(ind));
    assAltAver(i) = mean(assAlt_this);
    if sum(plotIDs == IDs(i)) > 0
        figure;
        hist(Bi_this)
        set(gca,'fontsize',15)
    end
    
end

rmInd = ~isfinite(Btrue);
IDs(rmInd) = [];
Btrue(rmInd) = [];
Btrue_sig(rmInd) = [];
assAltAver(rmInd) = [];
saveForm = [IDs,Btrue,Btrue_sig,assAltAver];
filename = ['Bfactors_',datestr(beginDate),'_',datestr(endDate),'.dat'];
csvwrite(filename,saveForm);
    
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