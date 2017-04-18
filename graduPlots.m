function [] = graduPlots()

load GRACE_2002-08-01_2003-07-31.mat
load('ilData.mat','rhoStruct')
load rhoStruct_il.mat

ensExample = ensMeans(:,end);
timeExample = assTimes(length(assTimes));
assWindow = round((assTimes(2)-assTimes(1)) * 24);
dt = 6; % hours
if mod(24,dt) ~= 0
    error('dt must be a harmonic of 24 hours!');
end
if mod(assWindow,dt) ~= 0
    error('dt must be a harmonic of the assimilation window length!');
end

duration = 2*7; %days
endTime = timeExample + duration;
assRMS = [];
ilRMS = [];
plotTimes = [];
t = timeExample;

rhoStruct.dTCoeff = assimiStruct.dTCoeff;
rhoStruct.T0Coeff = assimiStruct.T0Coeff;
rhoStruct.TexCoeff = assimiStruct.TexCoeff;
rhoStruct.OCoeff = assimiStruct.OCoeff;
rhoStruct.N2Coeff = assimiStruct.N2Coeff;
rhoStruct.HeCoeff = assimiStruct.HeCoeff;
rhoStruct.ArCoeff = assimiStruct.ArCoeff;
rhoStruct.O2Coeff = assimiStruct.O2Coeff;
rhoStruct.O_numBiases = assimiStruct.O_numBiases;
rhoStruct.N2_numBiases = assimiStruct.N2_numBiases;
rhoStruct.He_numBiases = assimiStruct.He_numBiases;
rhoStruct.Ar_numBiases = assimiStruct.Ar_numBiases;
rhoStruct.O2_numBiases = assimiStruct.O2_numBiases;

while t < endTime
    ind = t <= rhoStruct.timestamps & rhoStruct.timestamps < t+dt/24;
    removeInd = ~ind;
    S = removeDataPoints(rhoStruct, removeInd);
    S = computeVariablesForFit(S);
    prediction = exp(il_model_operator(ensExample, S, 1));
    rAss = rms(S.data./prediction-1);
    rIL = rms(S.data./ilRho(ind)-1);
    assRMS = [assRMS, rAss];
    ilRMS = [ilRMS, rIL];
    plotTimes = [plotTimes, t];
    
    t = t + dt/24;
end
figure;
plot(plotTimes, assRMS, plotTimes, ilRMS, 'linewidth',2.0);
legend('Data-assimilaatio','IL')
datetick('x')
ylabel('RMSE','fontsize',15)
set(gca,'fontsize',15)

Nduration = round(duration*24/dt);
Ndt = round(assWindow/dt) * (length(assTimes)-1) + 1;
N = Ndt + Nduration-1;
Sarray = repmat(S,N,1);
ilRMS = zeros(N,1);
t0 = assTimes(1);
ilRMS = [];
removeInd = false(size(rhoStruct.data));
removeInd(rhoStruct.champ) = true;
rhoStruct = removeDataPoints(rhoStruct, removeInd);
for i = 1:N
    t = t0 + (i-1)*dt/24;
    ind = t <= rhoStruct.timestamps & rhoStruct.timestamps < t+dt/24;
    removeInd = ~ind;
    S = removeDataPoints(rhoStruct, removeInd);
    S = computeVariablesForFit(S);
    Sarray(i) = S;
    
    rIL = rms(S.data./ilRho(ind)-1);
    ilRMS = [ilRMS, rIL];
end

leadTimes = [];

for i = 1:length(assTimes)
    t = assTimes(i);
    ens = ensMeans(:,i);
    
    assRMS = [];
    times = [];
    beginInd = round(assWindow/dt)*(i-1) + 1;
    endInd = beginInd + Nduration - 1;
    
    for k = beginInd : endInd
        prediction = exp(il_model_operator(ensExample, Sarray(k), 1));
        rAss = rms(Sarray(k).data./prediction-1);
        assRMS = [assRMS, rAss];
        times = [times, t];

        t = t + dt/24;
    end
    
    ilRMS_this = ilRMS(beginInd : endInd);
    
    leadTimeInd = find(assRMS > ilRMS_this, 1);
    if isempty(leadTimeInd)
        leadTime = duration;
    elseif leadTimeInd == 1
        leadTime = 0;
    else
        leadTime = times(leadTimeInd-1) - assTimes(i);
    end
    leadTimes = [leadTimes, leadTime*24];
end

figure;
hist(leadTimes)
xlabel('Ennusteaika [h]','fontsize',15)
set(gca,'fontsize',15)

end