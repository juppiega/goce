function [] = computeOptimalLagAndTau(satellite, globalMean)

load aeData.mat
ae = aeInterp;
t_ae = tInterp;

load('ilData.mat','rhoStruct');
if ~isfield(rhoStruct,'dst')
    rhoStruct = computeDst(rhoStruct);
    save('ilData.mat','rhoStruct','-append')
end

load('optCoeff.mat');

N = length(rhoStruct.data);
if strcmpi(satellite,'goce')
    removeInd = ~ismember(1:N, rhoStruct.goce); scaleFac = 0.9635;
elseif strcmpi(satellite,'champ')
    removeInd = ~ismember(1:N, rhoStruct.champ); scaleFac = 0.9578;
elseif strcmpi(satellite,'grace')
    removeInd = ~ismember(1:N, rhoStruct.grace); scaleFac = 0.9362;
elseif strcmpi(satellite,'all')
    removeInd = ismember(1:N, rhoStruct.swarm); scaleFac = 0.9034;
end

rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, false, false);
%[rhoStruct] = removeAndFixData(rhoStruct, 0, [], [], [], [], [], [], [],[], false); %orbAver

numBiasesStruct = struct('O', 5, 'N2', 6,...
     'He', 5, 'Ar', 2, 'O2', 0);

TexStruct.coeffInd = TexInd;
coeffStruct = struct('TexCoeff' , optCoeff(TexInd),... 
'OCoeff', optCoeff(OInd),...
'N2Coeff' , optCoeff(N2Ind),...
'HeCoeff' , optCoeff(HeInd),...
'O2Coeff' , optCoeff(O2Ind),...
'ArCoeff' , optCoeff(ArInd),...
'dTCoeff', dTCoeffs,...
'T0Coeff', T0Coeffs);

aeInt = rhoStruct.aeInt;
rhoStruct.aeInt = zeros(size(rhoStruct.aeInt))+20;
[modelRho] = computeComparisonData(rhoStruct, coeffStruct, numBiasesStruct);
modelRho = modelRho * scaleFac;

if globalMean
    [modelRho, times] = computeOrbitAverage(modelRho, rhoStruct.latitude, rhoStruct.timestamps);%orbAver

    obsRho = computeOrbitAverage(rhoStruct.data, rhoStruct.latitude, rhoStruct.timestamps); % orbAver
    aeIntAver = zeros(length(times),size(aeInt,2));
    for i = 1:size(aeInt,2)
        aeIntAver(:,i) = computeOrbitAverage(aeInt(:,i), rhoStruct.latitude, rhoStruct.timestamps);
    end
    dstAver = computeOrbitAverage(rhoStruct.dst, rhoStruct.latitude, rhoStruct.timestamps);

    % orbAver
    quietInd = dstAver > 75; %all(aeIntAver < 200, 2);
    modelRho(quietInd) = [];
    obsRho(quietInd) = [];

    rhoDiff = obsRho ./ modelRho - 1;
    rmInd = rhoDiff < 0;
    rhoDiff(rmInd) = [];

    maxLag = 7; 
    eTime = (1:1:50);
    crossCorrs = zeros(length(eTime), maxLag+1);

    for i = 1:length(eTime)
        aeInt = computeAEintegralExp(ae, t_ae, eTime(i))';
        aeInt = interp1(t_ae', aeInt, rhoStruct.timestamps);
        aeInt = computeOrbitAverage(aeInt, rhoStruct.latitude, rhoStruct.timestamps);
        aeInt(quietInd) = [];
        aeInt(rmInd) = [];

        crossCorr_this = xcorr(rhoDiff, aeInt, maxLag,'coeff');
        crossCorrs(i,:) = crossCorr_this(maxLag+1:end);
        if eTime(i) == 24
            a=1;
        end
    end

    [X,Y] = meshgrid(0:1:maxLag, eTime);
    Z = bsxfun(@rdivide, crossCorrs, crossCorrs(:,1));
    figure;
    surf(X,Y,Z,'edgecolor','none')
    view(2);
    colorbar;
    axis tight;
    set(gca,'fontsize',15);
    ylabel('e-fold time [h]','fontsize',15);
    xlabel('Lag [# orbits]','fontsize',15);

else
    zones = 0:15:90;
    maxLag = 6; 
    eTime = (1:3);
    crossCorrs = zeros(length(eTime), maxLag+1, length(zones)-1);
    
    quietInd = all(aeInt < 200, 2);
    rhoStruct = removeDataPoints(rhoStruct, quietInd, false, true, false, false);
    modelRho(quietInd) = [];
    for i = 1:length(zones)-1
        ind = zones(i) <= abs(rhoStruct.latitude) & abs(rhoStruct.latitude) < zones(i+1);
        rhoDiff = rhoStruct.data(ind) ./ modelRho(ind) - 1;
        rmInd = rhoDiff < 0;
        rhoDiff(rmInd) = [];
        
        for j = 1:length(eTime)
            aeInt = computeAEintegralExp(ae, t_ae, eTime(j))';
            aeInt = interp1(t_ae', aeInt, rhoStruct.timestamps(ind));
            aeInt(rmInd) = [];

            crossCorr_this = xcorr(rhoDiff, aeInt, maxLag,'coeff');
            crossCorrs(j,:,i) = crossCorr_this(maxLag+1:end);
        end
        
        [X,Y] = meshgrid(0:1:maxLag, eTime);
        figure;
        surf(X,Y,crossCorrs(:,:,i),'edgecolor','none')
        view(2);
        colorbar;
        axis tight;
        set(gca,'fontsize',15);
        title(num2str(mean(zones([i,i+1]))))
        ylabel('e-fold time [h]','fontsize',15);
        xlabel('Lag [# orbits]','fontsize',15);
    end
end

end