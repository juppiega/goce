function [  ] = fitIlModel(  )
rng(1, 'twister');

global numCoeffs;
global modelLbHeight;
modelLbHeight = 130;
numCoeffs = 5;

maxIter = 5;

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool();
end

% Check the existence of the data file.
if exist('ilData.mat', 'file')
    load ilData.mat
else
    error('File ilData.mat not found!')
end

[TempStruct, OStruct, rhoStruct] = removeAndFixData(TempStruct, OStruct, rhoStruct);
if ~exist('TexStruct', 'var')
    TexStruct = computeExosphericTemperatures(TempStruct);
end


opt = optimoptions('lsqnonlin', 'Jacobian', 'on', 'Algorithm', 'levenberg-marquardt');
ms = MultiStart('Display', 'off', 'UseParallel', true);
numStartPoints = 8;

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fprintf('\nFitting the IL Model\n\n')
for iter = 1:maxIter
    fprintf('Iteration %d:\n', iter)
    
    TexStruct = fitTex(TexStruct, opt, ms, numStartPoints);

    OStruct = fitSpecies(OStruct, TexStruct, opt, ms, numStartPoints);
    
    printProgressInfo(TexStruct, OStruct, rhoStruct);
    
    if iter < maxIter
        [TexStruct, OStruct] = estimateTexAndDensities(TexStruct, OStruct, rhoStruct);
    end

end
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end

function [TempStruct, OStruct, rhoStruct] = removeAndFixData(TempStruct, OStruct, rhoStruct)

% Make GOCE observations unbiased.
rhoStruct.data(rhoStruct.goce) = rhoStruct.data(rhoStruct.goce) * 1.23;
rhoStruct.numBiases = 0;
removeInd = ~ismember(1:length(rhoStruct.data), 1:10:length(rhoStruct.data)); % !!!!!!!!!! TESTAUS
rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, false, false);
goceWeight = 0.5*(1 - length(rhoStruct.goce) / length(rhoStruct.data));
champWeight = 1 - length(rhoStruct.champ)  / length(rhoStruct.data);
graceWeight = 1 - length(rhoStruct.grace) / length(rhoStruct.data);
rhoStruct.weights = zeros(length(rhoStruct.data),1);
rhoStruct.weights(rhoStruct.goce) = goceWeight;
rhoStruct.weights(rhoStruct.champ) = champWeight;
rhoStruct.weights(rhoStruct.grace) = graceWeight;

TempStruct.T0 = 507;
TempStruct.dT0 = 12.6;
% Remove bad temperature observations
removeInd = TempStruct.data <= TempStruct.T0 | TempStruct.data > 10000;
% For an unknown reason, there is one erroneous value in the local solar
% time present in the original ASCII data files.
removeInd(TempStruct.solarTime > 24 | TempStruct.solarTime < 0) = true(1);
% Remove observations with altitude > 900 km or < 250 or temperature > 3000 K.
removeInd(TempStruct.data > 3000 | TempStruct.altitude > 900 | TempStruct.altitude < 250) = true(1);
TempStruct = removeDataPoints(TempStruct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(TempStruct.de2) = 1;
satInd(TempStruct.aeC) = 2;
satInd(TempStruct.aeE) = 3;
satInd(removeInd) = [];
TempStruct.de2 = find(satInd == 1);
TempStruct.aeC = find(satInd == 2);
TempStruct.aeE = find(satInd == 3);
TempStruct.numBiases = 0;
TempStruct.dataEnd = length(TempStruct.data);

% Remove bad oxygen observations
removeInd = OStruct.data <= 1E6 | OStruct.data > 0.8E10;
OStruct = removeDataPoints(OStruct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(OStruct.de2) = 1;
satInd(OStruct.aeC) = 2;
satInd(OStruct.aeENace) = 3;
satInd(OStruct.aeEOss) = 4;
satInd(removeInd) = [];
OStruct.de2 = find(satInd == 1);
OStruct.aeC = find(satInd == 2);
OStruct.aeENace = find(satInd == 3);
OStruct.aeEOss = find(satInd == 4);
OStruct.numBiases = 4;
OStruct.biases = zeros(length(OStruct.data), OStruct.numBiases);
OStruct.biases(OStruct.de2, 1) = 1;
OStruct.biases(OStruct.aeC, 2) = 1;
OStruct.biases(OStruct.aeENace, 3) = 1;
OStruct.biases(OStruct.aeEOss, 4) = 1;
OStruct.dataEnd = length(OStruct.data);
OStruct.name = 'O';

TempStruct = computeGeopotentialHeight(TempStruct);
OStruct = computeGeopotentialHeight(OStruct);
rhoStruct = computeGeopotentialHeight(rhoStruct);

end

function [fixStruct] = removeDataPoints(fixStruct, removeInd, removeBiasMatRows, fixSatIndices, fixZ, fixWeights)

if (length(fixStruct.data) ~= length(removeInd))
    error('Lengths of struct data and remove indices (logical vector) must be the same!')
end

fixStruct.data(removeInd) = [];
fixStruct.timestamps(removeInd) = [];
fixStruct.latitude(removeInd) = [];
fixStruct.longitude(removeInd) = [];
fixStruct.solarTime(removeInd) = [];
fixStruct.altitude(removeInd) = [];
fixStruct.aeInt(removeInd,:) = [];
fixStruct.F(removeInd) = [];
fixStruct.FA(removeInd) = [];
fixStruct.apNow(removeInd) = [];
fixStruct.ap3h(removeInd) = [];
fixStruct.ap6h(removeInd) = [];
fixStruct.ap9h(removeInd) = [];
fixStruct.ap12To33h(removeInd) = [];
fixStruct.ap36To57h(removeInd) = [];
fixStruct.Ap(removeInd) = [];

if nargin > 2
    if removeBiasMatRows
        fixStruct.biases(removeInd,:) = [];
    end
    if fixSatIndices
        satInd = zeros(1, length(removeInd));
        satInd(fixStruct.goce) = 1;
        satInd(fixStruct.champ) = 2;
        satInd(fixStruct.grace) = 3;
        satInd(removeInd) = [];
        fixStruct.goce = find(satInd == 1);
        fixStruct.champ = find(satInd == 2);
        fixStruct.grace = find(satInd == 3);
    end
    if fixZ
        fixStruct.Z(removeInd) = [];
    end
    if fixWeights
        fixStruct.weights(removeInd) = [];
    end
end

end

function [addStruct] = computeGeopotentialHeight(addStruct)
% TODO: Lisaa leveyspiirin ja pyorimisen vaikutus.

R = 6371E3;
addStruct.Z = (R * addStruct.altitude) ./ (R + addStruct.altitude*1000);

end

function S = computeDensityRHS(S, Tex, dT0, T0)
global modelLbHeight
numBiasesOrig = S.numBiases; S.numBiases = 0;
S.numBiases = numBiasesOrig;
u2kg = 1.660538921E-27;
k = 1.38064852E-23;
g = 9.80665; % TODO: Onko oikein?
if strcmpi(S.name, 'O')
    molecMass = 16 * u2kg;
    alpha = 0;
else
    error('Incorrect name for gas species!')
end

sigma = dT0 ./ (Tex - T0);
T = Tex - (Tex - T0) .* exp(-sigma .* (S.Z - modelLbHeight));
gamma = molecMass * g ./ (k * sigma * 1E-3 .* Tex);
altTerm = (1 + gamma + alpha) .* log(T0 ./ T) - gamma .* sigma .* (S.Z - modelLbHeight);
S.rhs = log(S.data) - altTerm;

if any(~isfinite(S.rhs))
     a = 1;
end

end

function [TexStruct] = computeExosphericTemperatures(TempStruct)

TexVec = zeros(length(TempStruct.data),1);
% Set solution accuracy to 0.001.
options = optimoptions('fmincon', 'TolX', 1E-6, 'display', 'none');
% Loop over the temperature data.
Z = TempStruct.Z; data = TempStruct.data;
T0 = TempStruct.T0; dT0 = TempStruct.dT0;

targetCount = round(length(TexVec) / 1000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                'Computing Tex, ETA ', ...
                '. Now at ', ...
                'Completed in ' );
parfor i = 1:length(TexVec);
    % Because the Bates profile is nonlinear, use numerical root
    % finding.
    minFunc = @(Tex)TexFunc(Tex, Z(i), data(i), T0, dT0);
    TexVec(i) = fmincon(minFunc, 1000, [],[],[],[],data(i),10000,[], options);
    if mod(i,1000) == 0
        p.progress;
    end
end
p.stop;

TexStruct = TempStruct;
% Remove clear outliers.
removeInd = (TexVec - data) > 3000;
TexStruct = removeDataPoints(TexStruct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(TexStruct.de2) = 1;
satInd(TexStruct.aeC) = 2;
satInd(TexStruct.aeE) = 3;
satInd(removeInd) = [];
TexStruct.de2 = find(satInd == 1);
TexStruct.aeC = find(satInd == 2);
TexStruct.aeE = find(satInd == 3);
TexStruct.data = TexVec(~removeInd);
TexStruct.dataEnd = length(TexStruct.data);

save('ilData.mat', 'TexStruct', '-append')
end

% Function, whose root gives the Tex.
function y = TexFunc(Tex, Z, Tz, T0, dT0)
    %global modelLbHeight
    y = abs(Tex - (Tex - T0)*exp(-(Z-130) * dT0 / (Tex-T0)) - Tz);
end

function [addStruct] = computeVariablesForFit(addStruct)

x = cos(90 - addStruct.latitude);

% First degree functions.
P = legendre(1, x);
addStruct.P10 = P(1,:)';

% Second degree.
P = legendre(2, x);
addStruct.P20 = P(1,:)';

% Third degree.
P = legendre(3, x);
addStruct.P30 = P(1,:)';

% Fourth degree.
P = legendre(4, x);
addStruct.P40 = P(1,:)';

% Annual parameter.
[yr,~,~,~,~,~] = datevec(addStruct.timestamps);
yearVec = [yr, repmat([1,1,0,0,0], length(yr), 1)];
doy = addStruct.timestamps - datenum(yearVec) + 1;
addStruct.yv = 2*pi*(doy-1)/365;

end

function [removeStruct] = removeFitVariables(removeStruct)

removeStruct.P10 = [];
removeStruct.P20 = [];
removeStruct.P30 = [];
removeStruct.P40 = [];

removeStruct.yv = [];

end

function [TexStruct] = fitTex(TexStruct, options, multiStartSolver, startPointsForVar)
global numCoeffs

fprintf('Fitting Tex.\n')

initGuess = [1000, zeros(1,4)];
[G_lb, G_ub] = G_bounds();
lb = [500, 500*G_lb]; ub = [1500, 500*G_ub];

de = TexStruct.dataEnd;

TexStruct = computeVariablesForFit(TexStruct);
TexStruct.rhs = TexStruct.data;
Jacobian = ones(length(TexStruct.data), numCoeffs);
if de < length(TexStruct.data)
    numRhoTotDatPoints = length(TexStruct.data) - de;
    TexStruct.weights(1:de) = sqrt(norm(TexStruct.weights(de+1:end)) / numRhoTotDatPoints);
else
    TexStruct.weights = ones(length(TexStruct.data), 1);
end
fun = @(coeff)modelResidual(coeff, TexStruct, Jacobian);

problem = createOptimProblem('lsqnonlin', 'x0', initGuess, 'objective', fun, 'options', options);
numStartPoints = length(lb) * startPointsForVar;
startPoints = createRandomStartPoints(lb, ub, numStartPoints);

[optimalCoeff, fmin, flag, output, allmins] = run(multiStartSolver, problem, startPoints);

TexStruct.eval = @(S)evalTex(S, optimalCoeff);

TexStruct.optimalCoeff = optimalCoeff;
TexStruct.coeffErr = 1.96 * std(vertcat(allmins.X),1);

TexStruct = removeFitVariables(TexStruct);

end

function [F] = evalTex(S, coeff)

F = coeff(1) + G(coeff, S);

end

function [varStruct] = fitSpecies(varStruct, TexStruct, options, multiStartSolver, startPointsForVar)
global numCoeffs

fprintf('%s.\n', ['Fitting ', varStruct.name])
de = varStruct.dataEnd;
firstPass = de == length(varStruct.data);

varStruct.weights = ones(length(varStruct.data), 1);
varStruct = computeVariablesForFit(varStruct);
origNumBiases = varStruct.numBiases; varStruct.numBiases = 0;
Tex = TexStruct.eval(varStruct);
varStruct.numBiases = origNumBiases;
removeInd = Tex < TexStruct.T0 + 1; Tex(removeInd) = []; de = de - length(find(removeInd(1:de)));
fitStruct = removeDataPoints(varStruct, removeInd, true, ~firstPass, true, true);
fitStruct = computeVariablesForFit(fitStruct);
fitStruct = computeDensityRHS(fitStruct, Tex, TexStruct.dT0, TexStruct.T0);
if strcmpi(varStruct.name, 'O')
    ub = log(1E11);
    lb = log(0.5E10);
end
biasUb = ones(1, varStruct.numBiases) * log(1.0); ub = [ub, biasUb];
biasLb = ones(1, varStruct.numBiases) * log(0.5); lb = [lb, biasLb];
[G_lb, G_ub] = G_bounds();
lb = [lb, G_lb]; ub = [ub, G_ub];
initGuess = mean([lb;ub]);

Jacobian = ones(length(fitStruct.data), numCoeffs+fitStruct.numBiases);
Jacobian(:,2:fitStruct.numBiases+1) = fitStruct.biases;

if ~firstPass
    numRhoTotDatPoints = length(fitStruct.data) - de;
    fitStruct.weights(1:de) = sqrt(norm(fitStruct.weights(de+1:end)) / numRhoTotDatPoints);
end

fun = @(coeff)modelResidual(coeff, fitStruct, Jacobian);

problem = createOptimProblem('lsqnonlin', 'x0', initGuess, 'objective', fun, 'options', options);
numStartPoints = length(lb) * startPointsForVar;
startPoints = createRandomStartPoints(lb, ub, numStartPoints);

[optimalCoeff, fmin, flag, output, allmins] = run(multiStartSolver, problem, startPoints);

varStruct.evalLb = @(S)evalSpecies(S, optimalCoeff);

varStruct.optimalCoeff = optimalCoeff;
varStruct.coeffErr = 1.96 * std(vertcat(allmins.X),1);

varStruct = removeFitVariables(varStruct);

end

function [F] = evalSpecies(S, coeff)

F = exp(coeff(1) + G(coeff, S));

end

function [residual, Jacobian] = modelResidual(coeff, varStruct, Jacobian)

if nargout == 1
    Gvec = G(coeff, varStruct);
elseif nargout == 2
    [Gvec, Jacobian] = G(coeff, varStruct, Jacobian);
    Jacobian = bsxfun(@times, varStruct.weights, Jacobian);
else
    error('Requested wrong number of outputs!')
end

if varStruct.numBiases == 0
    residual = coeff(1) + Gvec - varStruct.rhs;
elseif varStruct.numBiases > 0
    residual = coeff(1) + sum(bsxfun(@times, coeff(2:varStruct.numBiases+1), varStruct.biases), 2) + Gvec - varStruct.rhs;
else
    error('Incorrect number of biases!')
end
residual = varStruct.weights .* residual;

end

function [result, Jacobian] = G(a, S, Jacobian)

k = S.numBiases + 1; % Counter, which helps adding terms.

% Latitude terms.
latitude = a(k+1)*S.P10 + a(k+2)*S.P20 + a(k+3)*S.P30 + a(k+4)*S.P40;
k = 4;

% % Solar activity terms.
% solarActivity = a(k+1)*S.F + a(k+2)*S.F.^2 + a(k+3)*S.FA + a(k+4)*S.FA.^2 + a(k+5)*S.F.*S.FA;
% k = k + 5;
% 
% % Annual terms.
% annual = (a(k+1) + a(k+2)*S.P20 + a(k+3)*S.P40) * (a(k+4)*sin(S.yv) + a(k+5)*cos(S.yv));

result = latitude;

if nargin == 3
    k = S.numBiases + 2;
    Jacobian(:,k:end) = [S.P10, S.P20, S.P30, S.P40]; % Latitude terms.
end

end

function [lb, ub] = G_bounds()

latitude = [1, 1, 1, 1];

ub = latitude;
lb = -ub;

end

function [startPoints] = createRandomStartPoints(lb, ub, numPoints)

if (any(lb > ub))
    error('Not all lower bound values were smaller than their corresponding upper bound!')
end

pd = makedist('Beta', 'a', 0.5, 'b', 0.5);

numVars = length(lb);
pointMat = zeros(numPoints, numVars);

boundMean = mean([lb;ub],1);
boundDiff = ub - lb;

for i = 1:numPoints
     for j = 1:numVars
         pointMat(i,j) = (boundDiff(j) * (random(pd)-0.5)) + boundMean(j);
     end
end

startPoints = CustomStartPointSet(pointMat);

end

function [TexStruct, OStruct] = estimateTexAndDensities(TexStruct, OStruct, rhoStruct)

fprintf('Estimating Tex and densities.\n')

numObs = length(rhoStruct.data);
fixedTex = ones(numObs, 1);
T0 = TexStruct.T0; dT0 = TexStruct.dT0;
Z = rhoStruct.Z;
rhoStruct = computeVariablesForFit(rhoStruct);
OlbDens = OStruct.evalLb(rhoStruct);
observations = rhoStruct.data;

options = optimset('TolX', 1E-6, 'Display', 'off');

targetCount = round(numObs / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                'Recomputing Tex, ETA ', ...
                '. Now at ', ...
                'Completed in ' );
parfor i = 1:numObs;
    % Because the Bates profile is nonlinear, use numerical root
    % finding.
    
    zeroFunc = @(Tex)(observations(i) - computeRho(T0, dT0, Tex, Z(i), OlbDens(i)))*1E15;
    fixedTex(i) = fzero(zeroFunc, T0+1, options);
    
    if mod(i,10000) == 0
        p.progress;
    end 
end
p.stop;

TexStruct = addAccelerometerMeasurements(TexStruct, fixedTex, rhoStruct);
removeInd = TexStruct.data < T0+1 | TexStruct.data > 9000 | isnan(TexStruct.data); removeInd(1:TexStruct.dataEnd) = false;
TexStruct = removeDataPoints(TexStruct, removeInd, false, true, true, true);

[~,newO] = computeRho(T0, dT0, fixedTex, Z, OlbDens);
OStruct = addAccelerometerMeasurements(OStruct, newO, rhoStruct, true);
OremInd = ismember(1:length(OStruct.data), find(removeInd == true) + (OStruct.dataEnd-TexStruct.dataEnd));
OStruct = removeDataPoints(OStruct, OremInd, true, true, true, true);

TexNans = find(~isfinite(TexStruct.data));
Onans = find(~isfinite(OStruct.data));

if length(TexNans)> 0 || length(Onans) > 0
    a = 1;
end

end

function addStruct = addAccelerometerMeasurements(addStruct, newData, rhoStruct, addZeroBias)

de = addStruct.dataEnd;

addStruct.data = [addStruct.data(1:de); newData];
addStruct.timestamps = [addStruct.timestamps(1:de); rhoStruct.timestamps];
addStruct.latitude = [addStruct.latitude(1:de); rhoStruct.latitude];
addStruct.longitude = [addStruct.longitude(1:de); rhoStruct.longitude];
addStruct.solarTime = [addStruct.solarTime(1:de); rhoStruct.solarTime];
addStruct.altitude = [addStruct.altitude(1:de); rhoStruct.altitude];
addStruct.aeInt = [addStruct.aeInt(1:de,:); rhoStruct.aeInt];
addStruct.F = [addStruct.F(1:de); rhoStruct.F];
addStruct.FA = [addStruct.FA(1:de); rhoStruct.FA];
addStruct.apNow = [addStruct.apNow(1:de); rhoStruct.apNow];
addStruct.ap3h = [addStruct.ap3h(1:de); rhoStruct.ap3h];
addStruct.ap6h = [addStruct.ap6h(1:de); rhoStruct.ap6h];
addStruct.ap9h = [addStruct.ap9h(1:de); rhoStruct.ap9h];
addStruct.ap12To33h = [addStruct.ap12To33h(1:de); rhoStruct.ap12To33h];
addStruct.ap36To57h = [addStruct.ap36To57h(1:de); rhoStruct.ap36To57h];
addStruct.Ap = [addStruct.Ap(1:de); rhoStruct.Ap];
addStruct.Z = [addStruct.Z(1:de); rhoStruct.Z];
addStruct.weights = [addStruct.weights(1:de); rhoStruct.weights];

addStruct.goce = rhoStruct.goce + de;
addStruct.champ = rhoStruct.champ + de;
addStruct.grace = rhoStruct.grace + de;

if nargin > 3 && addZeroBias == true
    addStruct.biases = [addStruct.biases(1:de,:); zeros(length(newData), addStruct.numBiases)];
end

end

function [rho, OnumDens] = computeRho(T0, dT0, Tex, Z, OlbDens)
%global modelLbHeight
sigma = dT0 ./ (Tex - T0);
g = 9.80665;
u2kg =  1.660538921E-27;
k = 1.38064852E-23;

T = Tex - (Tex - T0) .* exp(-sigma .* (Z - 130));

gamma_O = 16 * u2kg * g ./ (sigma*1E-3 .* k .* Tex);
f_O = (T0 ./ T).^(1+gamma_O) .* exp(-sigma .* (Z - 130) .* gamma_O);
OnumDens = OlbDens.*f_O; % [1/cm^3]


rho = (16*OnumDens) * u2kg * 1E6; % [kg/m^3]

end

function printProgressInfo(TexStruct, OStruct, rhoStruct)

T0 = TexStruct.T0; dT0 = TexStruct.dT0;
Z = rhoStruct.Z;
rhoStruct = computeVariablesForFit(rhoStruct);
OlbDens = OStruct.evalLb(rhoStruct);
Tex = TexStruct.eval(rhoStruct);
observations = rhoStruct.data;

rhoModel = computeRho(T0, dT0, Tex, Z, OlbDens);
ratio = rhoStruct.weights .* observations ./ rhoModel;

modelRms = rms(ratio);
modelStd = std(ratio);
obsToModel = mean(ratio);

fprintf('rms       std       O/C\n')
fprintf('%f        %f        %f\n\n', modelRms, modelStd, obsToModel);

end