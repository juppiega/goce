function [  ] = fitIlModel(  )
rng(1, 'twister');
%import java.lang.*
%r = Runtime.getRuntime;
%numThreads = r.availableProcessors;
numThreads = 64;

global numCoeffs;
global modelLbHeight;
modelLbHeight = 130;
numCoeffs = 108;

clear mex;

maxIter = 10;

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(24);
end

% Check the existence of the data file.
if exist('ilData.mat', 'file')
    load ilData.mat
else
    error('File ilData.mat not found!')
end

[TempStruct, OStruct, N2Struct, HeStruct, rhoStruct] = removeAndFixData(TempStruct, OStruct, N2Struct, HeStruct, rhoStruct);
if ~exist('TexStruct', 'var')
    TexStruct = computeExosphericTemperatures(TempStruct);
end


opt = optimoptions('lsqnonlin', 'Jacobian', 'on', 'Algorithm', 'trust-region-reflective', 'TolFun', 1E-3, ...
                 'TolX', 1E-7, 'Display', 'iter');
%opt = psoptimset('Display', 'iter', 'CompletePoll', 'on', 'InitialMeshSize', 1E-1, 'TolX', 1E-11, 'TolMesh', 1E-11, 'TolFun', 1E-5, 'UseParallel', true, 'Cache', 'on', 'CacheTol', 1E-12, 'MaxMeshSize', 1.0);


dataLen = length(TexStruct.data) + length(OStruct.data) + length(N2Struct.data) + length(HeStruct.data) + length(rhoStruct.data);
%opt = optimoptions('fmincon', 'Display', 'iter', 'FinDiffType', 'central', 'MaxFunEvals', 1000*numCoeffs, 'TolFun', 1E-3, 'Algorithm', 'sqp', ...
%    'TolX', 1E-11, 'UseParallel', true, 'Hessian', 'bfgs', 'ScaleProblem', 'obj-and-constr', 'OutputFcn', @(x,y,z)outfun(x,y,z,dataLen));

ms = MultiStart('Display', 'iter', 'UseParallel', true);
numStartPoints = 1;

% % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% fprintf('\nFitting the IL Model\n\n')
% for iter = 1:maxIter
%     fprintf('\nIteration %d:\n', iter)
%     
%     TexStruct = fitTex(TexStruct, opt, ms, numStartPoints);
% 
%     OStruct = fitSpecies(OStruct, TexStruct, opt, ms, numStartPoints);
%     
%     printProgressInfo(TexStruct, OStruct, rhoStruct);
%     
%     if iter < maxIter
%         [TexStruct, OStruct] = estimateTexAndDensities(TexStruct, OStruct, rhoStruct);
%     end
% 
% end
% % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fitModelVariables(TexStruct, OStruct, N2Struct, HeStruct, rhoStruct, opt, ms, numStartPoints, numThreads)

end

function [TempStruct, OStruct, N2Struct, HeStruct, rhoStruct] = removeAndFixData(TempStruct, OStruct, N2Struct, HeStruct, rhoStruct)

% Make GOCE observations unbiased.
rhoStruct.data(rhoStruct.goce) = rhoStruct.data(rhoStruct.goce) * 1.23;
rhoStruct.numBiases = 0;
removeInd = ~ismember(1:length(rhoStruct.data), 1:10:length(rhoStruct.data)) | rhoStruct.data' <= 0; % !!!!!!!!!! TESTAUS
rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, false, false);
goceWeight = 0.5*(1 - length(rhoStruct.goce) / length(rhoStruct.data));
champWeight = 1 - length(rhoStruct.champ)  / length(rhoStruct.data);
graceWeight = 1 - length(rhoStruct.grace) / length(rhoStruct.data);
rhoStruct.weights = zeros(length(rhoStruct.data),1);
rhoStruct.weights(rhoStruct.goce) = goceWeight;
rhoStruct.weights(rhoStruct.champ) = champWeight;
rhoStruct.weights(rhoStruct.grace) = graceWeight;
rhoStruct.numBiases = 0;

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
OStruct.numBiases = 0; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OStruct.biases = zeros(length(OStruct.data), OStruct.numBiases);
OStruct.biases(OStruct.de2, 1) = 1;
OStruct.biases(OStruct.aeC, 2) = 1;
OStruct.biases(OStruct.aeENace, 3) = 1;
OStruct.biases(OStruct.aeEOss, 4) = 1;
OStruct.dataEnd = length(OStruct.data);
OStruct.name = 'O';

% Remove N2 observations.
removeInd = N2Struct.data <= 1E2 | N2Struct.data > 0.5E12 | N2Struct.altitude > 600 | N2Struct.altitude < 140;
N2Struct = removeDataPoints(N2Struct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(N2Struct.de2) = 1;
satInd(N2Struct.aeC) = 2;
satInd(N2Struct.aeENace) = 3;
satInd(N2Struct.aeEOss) = 4;
satInd(N2Struct.aeros) = 5;
satInd(removeInd) = [];
N2Struct.de2 = find(satInd == 1);
N2Struct.aeC = find(satInd == 2);
N2Struct.aeENace = find(satInd == 3);
N2Struct.aeEOss = find(satInd == 4);
N2Struct.aeros = find(satInd == 5);
N2Struct.numBiases = 0;
N2Struct.dataEnd = length(N2Struct.data);
N2Struct.name = 'N2';

% Remove He observations
removeInd = HeStruct.data <= 1E5 | HeStruct.data > 1E10;
HeStruct = removeDataPoints(HeStruct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(HeStruct.de2) = 1;
satInd(HeStruct.aeC) = 2;
satInd(HeStruct.aeENace) = 3;
satInd(HeStruct.aeEOss) = 4;
satInd(HeStruct.aeros) = 5;
satInd(removeInd) = [];
HeStruct.de2 = find(satInd == 1);
HeStruct.aeC = find(satInd == 2);
HeStruct.aeENace = find(satInd == 3);
HeStruct.aeEOss = find(satInd == 4);
HeStruct.aeros = find(satInd == 5);
HeStruct.numBiases = 0;
HeStruct.dataEnd = length(HeStruct.data);
HeStruct.name = 'He';

TempStruct = computeGeopotentialHeight(TempStruct);
OStruct = computeGeopotentialHeight(OStruct);
N2Struct = computeGeopotentialHeight(N2Struct);
HeStruct = computeGeopotentialHeight(HeStruct);
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
%global modelLbHeight
numBiasesOrig = S.numBiases; S.numBiases = 0;
S.numBiases = numBiasesOrig;
u2kg = 1.660538921E-27;
k = 1.38064852E-23;
g = 9.80665; % TODO: Onko oikein?
if strcmpi(S.name, 'O')
    molecMass = 16 * u2kg;
    alpha = 0;
elseif strcmpi(S.name, 'N2')
    molecMass = 28 * u2kg;
    alpha = 0;
elseif strcmpi(S.name, 'He')
    molecMass = 4 * u2kg;
    alpha = -0.38;
else
    error('Incorrect name for gas species!')
end

sigma = dT0 ./ (Tex - T0);
T = Tex - (Tex - T0) .* exp(-sigma .* (S.Z - 130));
gamma = molecMass * g ./ (k * sigma * 1E-3 .* Tex);
altTerm = (1 + gamma + alpha) .* log(T0 ./ T) - gamma .* sigma .* (S.Z - 130);
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

function [removeStruct] = removeFitVariables(removeStruct)

removeStruct.P10 = [];
removeStruct.P20 = [];
removeStruct.P30 = [];
removeStruct.P40 = [];

removeStruct.yv = [];

end

function [lb, ub] = G_bounds()

latitude = ones(1,14);
solarActivity = ones(1,5);
annual = ones(1,15); annual(1) = 0.001; annual(8) = 0.0001;
diurnal = ones(1, 21); diurnal([8, 19]) = 0.001;
semidiurnal = ones(1,16);
terdiurnal = ones(1,8);
quaterdiurnal = ones(1,2);
geomagnetic = ones(1,26); geomagnetic([1,9,15,21,18,19,24,25]) = 0.001;

ub = [latitude, solarActivity, annual, diurnal, semidiurnal, terdiurnal, quaterdiurnal, geomagnetic];
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

function printModelStats(TexStruct, OStruct, rhoStruct)

T0 = TexStruct.T0; dT0 = TexStruct.dT0;
Z = rhoStruct.Z;
rhoStruct = computeVariablesForFit(rhoStruct);
OlbDens = OStruct.evalLb(rhoStruct);
Tex = TexStruct.eval(rhoStruct);
observations = rhoStruct.data;

rhoModel = computeRho(T0, dT0, Tex, Z, OlbDens);
ratio = rhoStruct.weights .* observations ./ rhoModel;
failedInd = (Tex <= T0 | ~isfinite(rhoModel));
ratio = ratio(~failedInd);

modelRms = rms(ratio);
modelStd = std(ratio);
obsToModel = mean(ratio);

numFailed = sum(failedInd);

fprintf('rms             std             O/C\n')
fprintf('%f        %f        %f          %d        %E        %E\n\n', modelRms, modelStd, obsToModel, numFailed, max(rhoModel(~failedInd)), min(rhoModel(~failedInd)));

end

function [residual] = computeSpeciesResidual(varStruct, Tex, dT0, T0, coeff)

varStruct = computeDensityRHS(varStruct, Tex, dT0, T0);
Gvec = G(coeff, varStruct);

if varStruct.numBiases == 0
    residual = (varStruct.rhs ./ max(coeff(1) + Gvec, 1)) - 1;
elseif varStruct.numBiases > 0
    residual = (varStruct.rhs ./ max(coeff(1) + sum(bsxfun(@times, coeff(2:varStruct.numBiases+1), varStruct.biases), 2) + Gvec, 1)) - 1;
else
    error('Incorrect number of biases!')
end

end

function [residual, Jacobian] = modelMinimizationFunction(TexStruct, OStruct, N2Struct, HeStruct, rhoStruct, weights, tolX, coeff)

dataLen = length(TexStruct.data) + length(OStruct.data) + length(N2Struct.data) + length(HeStruct.data) + length(rhoStruct.data);
residual = zeros(dataLen, 1);

TexMesuredEstimate = clamp(TexStruct.T0+1, evalTex(TexStruct, coeff(TexStruct.coeffInd)), 5000);
residInd = 1:length(TexStruct.data);
residual(residInd) = TexStruct.data./TexMesuredEstimate - 1;

[Tex, dT0, T0] = findTempsForFit(OStruct, TexStruct, coeff);
residInd = residInd(end) + (1:length(OStruct.data));
residual(residInd) = computeSpeciesResidual(OStruct, Tex, dT0, T0, coeff(OStruct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(N2Struct, TexStruct, coeff);
residInd = residInd(end) + (1:length(N2Struct.data));
residual(residInd) = computeSpeciesResidual(N2Struct, Tex, dT0, T0, coeff(N2Struct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(HeStruct, TexStruct, coeff);
residInd = residInd(end) + (1:length(HeStruct.data));
residual(residInd) = computeSpeciesResidual(HeStruct, Tex, dT0, T0, coeff(HeStruct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(rhoStruct, TexStruct, coeff);
residInd = residInd(end) + (1:length(rhoStruct.data));
OlbDens = clamp(10, evalSpecies(rhoStruct, coeff(OStruct.coeffInd)), 1E20);
N2lbDens = clamp(10, evalSpecies(rhoStruct, coeff(N2Struct.coeffInd)), 1E20);
HelbDens = clamp(10, evalSpecies(rhoStruct, coeff(HeStruct.coeffInd)), 1E20);
modelRho = clamp(1E-20, computeRho(T0, dT0, Tex, rhoStruct.Z, OlbDens, N2lbDens, HelbDens), 0.1);
residual(residInd) = (log(rhoStruct.data)./log(modelRho)) - 1;

if any(abs(residual) > 100)
    a=1;
end

residual = weights .* residual;

if nargout == 2
    fun = @(X)modelMinimizationFunction(TexStruct, OStruct, N2Struct, HeStruct, rhoStruct, weights, tolX, X);
    Jacobian = computeJAC(fun, coeff, dataLen, tolX);
end

end

function JAC = computeJAC(fun, x, dataLen, tolX)

JAC = zeros(dataLen, length(x));
dx = max(0.25*tolX*abs(x), 1E-10);

parfor i = 1:length(x)
    xForw = x; xBackw = x;
    xForw(i) = x(i) + dx(i);
    xBackw(i) = x(i) - dx(i);
    
    F_forw = feval(fun, xForw);
    F_backw = feval(fun, xBackw);
    
    result = (F_forw - F_backw) / (2*dx(i));
    
    infInd = ~isfinite(result);
    if any(infInd)
        result(infInd) = mean(result(~infInd));
    end
    
    JAC(:,i) = result;
end

end

function [] = fitModelVariables(TexStruct, OStruct, N2Struct, HeStruct, rhoStruct, options, multiStartSolver, numStartPoints, numThreads)
global numCoeffs;

dataLen = length(TexStruct.data) + length(OStruct.data) + length(N2Struct.data) + length(HeStruct.data) + length(rhoStruct.data);

TexStruct = computeVariablesForFit(TexStruct);
OStruct = computeVariablesForFit(OStruct);
N2Struct = computeVariablesForFit(N2Struct);
HeStruct = computeVariablesForFit(HeStruct);
rhoStruct = computeVariablesForFit(rhoStruct);

weights = computeWeights(TexStruct, OStruct, N2Struct, HeStruct, rhoStruct);

[G_lb, G_ub] = G_bounds();
G_lb = 4 * G_lb; G_ub = 4 * G_ub;

TexStruct.coeffInd = 1:numCoeffs;
lb = [500, 125*G_lb]; ub = [1500, 125*G_ub];

OStruct.coeffInd = TexStruct.coeffInd(end) + (1:numCoeffs+OStruct.numBiases);
lb = [lb, log(0.5E10), G_lb]; % MUISTA LISATA BIASET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ub = [ub, log(1E11), G_ub];

N2Struct.coeffInd = OStruct.coeffInd(end) + (1:numCoeffs+N2Struct.numBiases);
lb = [lb, log(0.5E11), G_lb]; % MUISTA LISATA BIASET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ub = [ub, log(1E12), G_ub];

HeStruct.coeffInd = N2Struct.coeffInd(end) + (1:numCoeffs+HeStruct.numBiases);
lb = [lb, log(0.5E7), G_lb]; % MUISTA LISATA BIASET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ub = [ub, log(1E8), G_ub];

startPoints = createRandomStartPoints(lb, ub, numStartPoints);
%initGuess = list(startPoints);
ind = ub < 0.5;
ind(1:numCoeffs) = ub(1:numCoeffs) < 100;
initGuess = mean([lb;ub]);% - 0.001;
initGuess(ind) = ub(ind);
ub(ind) = mode(ub(~ind));
TexInd = 2:numCoeffs; ub(TexInd) = mode(ub(TexInd));
lb = -ub;
tolX = options.TolX;
%fun = @(coeff)sum(modelMinimizationFunction(TexStruct, OStruct, N2Struct, HeStruct, rhoStruct, weights, tolX, coeff).^2);
fun = @(coeff)modelMinimizationFunction(TexStruct, OStruct, N2Struct, HeStruct, rhoStruct, weights, tolX, coeff);
%problem = createOptimProblem('lsqnonlin', 'x0', initGuess, 'objective', fun, 'options', options);

%tic;[optCoeff, fmin, flag, output, allmins] = run(multiStartSolver, problem, startPoints);toc;

%tic;[optCoeff, fval, flag, output] = patternsearch(fun, initGuess, [],[],[],[], lb, ub, [], options);toc
%tic;[optCoeff, fval, flag, output] = patternsearch(fun, initGuess, [],[],[],[], [], [], [], options);toc

%tic;[optCoeff, fval, flag, output] = fmincon(fun, initGuess, [],[],[],[], lb, ub, [], options);toc

setenv('OMP_NUM_THREADS', num2str(numThreads))
disp('Calling LM solver')
tic;optCoeff = levenbergMarquardt_mex(TexStruct, OStruct, N2Struct, HeStruct, rhoStruct, weights, initGuess);toc;
%[comp] = fun(initGuess); %disp([comp(1), optCoeff]);
%JTJ_diag = diag(J'*J);

%tic; [optCoeff] = lsqnonlin(fun, initGuess, lb, ub, options);toc;

save('optCoeff.mat', 'optCoeff', '-v7.3');
TexInd = TexStruct.coeffInd; save('optCoeff.mat', 'TexInd', '-append');
HeInd = HeStruct.coeffInd; save('optCoeff.mat', 'HeInd', '-append');
OInd = OStruct.coeffInd; save('optCoeff.mat', 'OInd', '-append');
N2Ind = N2Struct.coeffInd; save('optCoeff.mat', 'N2Ind', '-append');


end

function weights = computeWeights(TexStruct, OStruct, N2Struct, HeStruct, rhoStruct)

dataLen = length(TexStruct.data) + length(OStruct.data) + length(N2Struct.data) + length(HeStruct.data) + length(rhoStruct.data);
TempAndSpectrometerLen = length(TexStruct.data) + length(OStruct.data) + length(N2Struct.data) + length(HeStruct.data);
weights = ones(dataLen, 1);

meanRhoWeight = mean(rhoStruct.weights); numRho = length(rhoStruct.weights);

w = meanRhoWeight * (numRho / TempAndSpectrometerLen);
wInd = 1:TempAndSpectrometerLen;
weights(wInd) = w;
weights(1:length(TexStruct.data)) = 0.5 * w;

weights(wInd(end)+1:end) = rhoStruct.weights;

ae16h = [TexStruct.aeInt(:,4); OStruct.aeInt(:,4); N2Struct.aeInt(:,4); HeStruct.aeInt(:,4); rhoStruct.aeInt(:,4)];
aeNormalized = 1 + (2 * ae16h / max(ae16h));
weights = weights .* aeNormalized;

weights = sqrt(weights);

end

function stop = outfun(x, optimValues, state, dataLen)

t = sqrt(abs(optimValues.firstorderopt)/dataLen);
fprintf('%s\n', ['Normalized first-order opt.: ', num2str(t)])
stop = t < 0.005;

end