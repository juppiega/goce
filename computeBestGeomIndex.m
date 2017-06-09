function [] = computeBestGeomIndex(satellite)

load('ilData.mat','rhoStruct');
load('optCoeff.mat');

N = length(rhoStruct.data);
if strcmpi(satellite,'goce')
    removeInd = ~ismember(1:N, rhoStruct.goce); scaleFac = 0.9635;
elseif strcmpi(satellite,'champ')
    removeInd = ~ismember(1:N, rhoStruct.champ); scaleFac = 0.9578;
elseif strcmpi(satellite,'grace')
    removeInd = ~ismember(1:N, rhoStruct.grace); scaleFac = 0.9362;
elseif strcmpi(satellite,'all')
    removeInd = false(N,1); scaleFac = 0.9034;
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
[modelRho, times] = computeOrbitAverage(modelRho, rhoStruct.latitude, rhoStruct.timestamps);%orbAver

% orbAver
aeIntAver = zeros(length(times),size(aeInt,2));
for i = 1:size(aeInt,2)
    aeIntAver(:,i) = computeOrbitAverage(aeInt(:,i), rhoStruct.latitude, rhoStruct.timestamps);
end

obsRho = computeOrbitAverage(rhoStruct.data, rhoStruct.latitude, rhoStruct.timestamps); % orbAver

% orbAver
quietInd = all(aeIntAver < 500, 2);
modelRho(quietInd) = [];
aeIntAver(quietInd,:) = [];
obsRho(quietInd) = [];
times(quietInd) = [];

rhoDiff = obsRho ./ modelRho - 1;

aeIntAver(rhoDiff < 0,:) = [];
times(rhoDiff < 0) = [];
rhoDiff(rhoDiff < 0) = [];

mdl = stepwiselm(aeIntAver, rhoDiff, 'upper','quadratic','intercept',true,'PEnter',0.33,'PRemove',0.5)
predDiff = predict(mdl,aeIntAver);
figure;
plot(times,predDiff, times,rhoDiff);
legend('pred','obs')
datetick('x')

end
