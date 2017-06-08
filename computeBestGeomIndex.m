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
end

rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, false, false);
[rhoStruct] = removeAndFixData(rhoStruct, 0, [], [], [], [], [], [], [], [], false);

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
[quietRho] = computeComparisonData(rhoStruct, coeffStruct, numBiasesStruct);
quietRho = quietRho * scaleFac;

rhoStruct.aeInt = aeInt;

rhoDiff = rhoStruct.data ./ quietRho - 1;

aeInt(rhoDiff < 0,:) = [];
rhoDiff(rhoDiff < 0) = [];

end