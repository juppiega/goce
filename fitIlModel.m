function [  ] = fitIlModel( recomputeTex, recomputeLbTemp, recomputeDT, recomputeQuietModel, recomputeStormModel, recomputeAlsoInsign, fitSimultaneous, optimizedMex, subsampPercent )

if optimizedMex
    mex -O FCFLAGS="\$FCFLAGS -std=f2008" -output levenbergMarquardt_mex lmSolver.F90 levenbergMarquardt_mex.F90 -llapack
    mex -O FCFLAGS="\$FCFLAGS -std=f2008" -output lbFit_mex lmSolver.F90 lbFit_mex.F90 -llapack
else
    mex -v -g FFLAGS='$FFLAGS -fcheck=all' -output levenbergMarquardt_mex lmSolver.F90 levenbergMarquardt_mex.F90 -llapack
    mex -v -g FFLAGS='$FFLAGS -fcheck=all' -output lbFit_mex lmSolver.F90 lbFit_mex.F90 -llapack
end

rng(1, 'twister');
%import java.lang.*
%r = Runtime.getRuntime;
%numThreads = r.availableProcessors;
numThreads = 64;
aeThreshold = 0;

global numCoeffs;
numCoeffs = 141;

clear mex;
% 
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

quietData = true;
[~, ~, ~, ~, ~, ~, ~, lbDTStruct, lbT0Struct] = ...
    removeAndFixData(rhoStruct, aeThreshold, TempStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct, quietData);


if ~exist('T0Coeffs', 'var') || recomputeLbTemp
    [T0Coeffs, JTWJ_T0] = fitLbTemerature(lbT0Struct, subsampPercent); 
    save('ilData.mat', 'T0Coeffs', '-append')
    save('ilData.mat', 'JTWJ_T0', '-append')
end

if ~exist('dTCoeffs', 'var') || recomputeDT
    [dTCoeffs, JTWJ_dT] = fitTemeratureGradient(lbDTStruct, subsampPercent); 
    save('ilData.mat', 'dTCoeffs', '-append')
    save('ilData.mat', 'JTWJ_dT', '-append')
end

if ~exist('TexStruct', 'var') || recomputeTex
    [~, TempStruct] = ...
    removeAndFixData(rhoStruct, aeThreshold, TempStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct);

    TexStruct = computeExosphericTemperatures(TempStruct, dTCoeffs, T0Coeffs);
    save('ilData.mat', 'TexStruct', '-append')
end


opt = optimoptions('lsqnonlin', 'Jacobian', 'on', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1E-5, ...
                 'TolX', 1E-6, 'Display', 'iter-detailed');
%opt = psoptimset('Display', 'iter', 'CompletePoll', 'on', 'InitialMeshSize', 1E-1, 'TolX', 1E-11, 'TolMesh', 1E-11, 'TolFun', 1E-5, 'UseParallel', true, 'Cache', 'on', 'CacheTol', 1E-12, 'MaxMeshSize', 1.0);


%dataLen = length(TexStruct.data) + length(OStruct.data) + length(N2Struct.data) + length(HeStruct.data) + length(rhoStruct.data);
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

if fitSimultaneous || recomputeQuietModel
    load ilData.mat
    quietData = true;
    [rhoStruct, TempStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct] = ...
    removeAndFixData(rhoStruct, aeThreshold, TempStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct, quietData);
    
    [rhoStruct, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct] = ...
    subsampleStructs(rhoStruct, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, subsampPercent);
    
    fitModelVariables(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, opt, ms, numStartPoints, numThreads, quietData, recomputeAlsoInsign, fitSimultaneous)
    
    fprintf('Quiet time fitted.\n')
end

if recomputeStormModel && ~fitSimultaneous
    load ilData.mat
    quietData = false;
    [rhoStruct, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct] = ...
    removeAndFixData(rhoStruct, aeThreshold, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct, quietData);

    [rhoStruct, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct] = ...
    subsampleStructs(rhoStruct, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, subsampPercent);

    load quietCoeffs.mat
    fitModelVariables(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, opt, ms, numStartPoints, numThreads, quietData, recomputeAlsoInsign, fitSimultaneous, optCoeff)
    fprintf('Storm time fitted.\n')
end

end

function S = computeDensityRHS(S, Tex, dT0, T0)
%global modelLbHeight
numBiasesOrig = S.numBiases; S.numBiases = 0;
S.numBiases = numBiasesOrig;
u2kg = 1.660538921E-27;
k = 1.38064852E-23;
g = 9.418; 
if strcmpi(S.name, 'O')
    molecMass = 16 * u2kg;
    alpha = 0;
elseif strcmpi(S.name, 'N2')
    molecMass = 28 * u2kg;
    alpha = 0;
elseif strcmpi(S.name, 'O2')
    molecMass = 32 * u2kg;
    alpha = 0;
elseif strcmpi(S.name, 'Ar')
    molecMass = 40 * u2kg;
    alpha = 0;
elseif strcmpi(S.name, 'He')
    molecMass = 4 * u2kg;
    alpha = -0.38;
else
    error('Incorrect name for gas species!')
end

sigma = dT0 ./ (Tex - T0);
T = Tex - (Tex - T0) .* exp(-sigma .* (S.Z));
gamma = molecMass * g ./ (k * sigma * 1E-3 .* Tex);
altTerm = (1 + gamma + alpha) .* log(T0 ./ T) - gamma .* sigma .* (S.Z);
S.rhs = log(S.data) - altTerm;

if any(~isfinite(S.rhs))
     a = 1;
end

end

function [TexStruct] = computeExosphericTemperatures(TempStruct, dTCoeffs, T0Coeffs)

fprintf('%s\n', 'Computing exospheric temperatures')

TexVec = zeros(length(TempStruct.data),1);
% Set solution accuracy to 0.001.
options = optimoptions('fmincon', 'TolX', 1E-6, 'display', 'none');
% Loop over the temperature data.
S = computeVariablesForFit(TempStruct);
Z = S.Z;
T0 = evalT0(S, T0Coeffs); 
data = S.data; data = clamp(T0+1, data, 10000);
dT0 = clamp(1, evalDT(S, dTCoeffs), 30);

targetCount = round(length(TexVec) / 1000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                'Computing Tex, ETA ', ...
                '. Now at ', ...
                'Completed in ' );
parfor i = 1:length(TexVec);
    % Because the Bates profile is nonlinear, use numerical root
    % finding.
    minFunc = @(Tex)TexFunc(Tex, Z(i), data(i), T0(i), dT0(i));
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
TexStruct.Z = TexStruct.Z(~removeInd);

end

function [dTCoeffs, JTWJ] = fitTemeratureGradient(lbDTStruct, subsampPercent)

fprintf('%s\n', 'Fitting Temperature gradient')

subsamp = round(1 / (subsampPercent / 100));

N = length(lbDTStruct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
lbDTStruct = removeDataPoints(lbDTStruct, removeInd);

lbDTStruct.data = [lbDTStruct.data(1:2:end); lbDTStruct.data(1:2:end)];
lbDTStruct.timestamps = [lbDTStruct.timestamps(1:2:end); lbDTStruct.timestamps(1:2:end) + 180];
lbDTStruct.latitude = [lbDTStruct.latitude(1:2:end); -lbDTStruct.latitude(1:2:end)];
lbDTStruct.longitude = [lbDTStruct.longitude(1:2:end); lbDTStruct.longitude(1:2:end)];
lbDTStruct.solarTime = [lbDTStruct.solarTime(1:2:end); lbDTStruct.solarTime(1:2:end)];
lbDTStruct.altitude = [lbDTStruct.altitude(1:2:end); lbDTStruct.altitude(1:2:end)];
lbDTStruct.aeInt = [lbDTStruct.aeInt(1:2:end); lbDTStruct.aeInt(1:2:end)];
lbDTStruct.F = [lbDTStruct.F(1:2:end); lbDTStruct.F(1:2:end)];
lbDTStruct.FA = [lbDTStruct.FA(1:2:end); lbDTStruct.FA(1:2:end)];
lbDTStruct.apNow = [lbDTStruct.apNow(1:2:end); lbDTStruct.apNow(1:2:end)];
lbDTStruct.ap3h = [lbDTStruct.ap3h(1:2:end); lbDTStruct.ap3h(1:2:end)];
lbDTStruct.ap6h = [lbDTStruct.ap6h(1:2:end); lbDTStruct.ap6h(1:2:end)];
lbDTStruct.ap9h = [lbDTStruct.ap9h(1:2:end); lbDTStruct.ap9h(1:2:end)];
lbDTStruct.ap12To33h = [lbDTStruct.ap12To33h(1:2:end); lbDTStruct.ap12To33h(1:2:end)];
lbDTStruct.ap36To57h = [lbDTStruct.ap36To57h(1:2:end); lbDTStruct.ap36To57h(1:2:end)];
lbDTStruct.Ap = [lbDTStruct.Ap(1:2:end); lbDTStruct.Ap(1:2:end)];

lbDTStruct = computeVariablesForFit(lbDTStruct);

latitude = zeros(1,14);
solarActivity = zeros(1,5);
annual = zeros(1,32); annual([1,7,12,16,19,25,30]) = 0.001;
diurnal = zeros(1, 21); diurnal([8, 19]) = 0.001;
semidiurnal = zeros(1,16);
terdiurnal = zeros(1,8);
quaterdiurnal = zeros(1,2);
longitudinal = zeros(1,13); longitudinal([2,5,9,12]) = 1E-4;
dTCoeffs = [zeros(1,1), latitude, solarActivity, annual, diurnal, semidiurnal, ...
    terdiurnal, quaterdiurnal, longitudinal];

dTCoeffs(1) = mean(lbDTStruct.data);


opt = optimoptions('lsqnonlin', 'Jacobian', 'on', 'Algorithm', 'Levenberg-Marquardt', 'TolFun', 1E-8, ...
                 'TolX', 1E-8, 'Display', 'iter', 'initDamping', 1E8, 'OutputFcn', @outfun);

fun = @(X) temperatureGradientMinimization(lbDTStruct, X);
%[fval, JAC] = fun(dTCoeffs);
%JTJ_diag = diag(JAC'*JAC);
%[dTCoeffs] = lsqnonlin(fun, dTCoeffs, [], [], opt);

tolX = 1E-8;
tolFun = 1E-5;
tolOpt = 1E-4;
lambda0 = 1E-2;
weights = ones(size(lbDTStruct.data));
endSemidiurnal = 89;
symmAnn = 21:38;
paramsToFit = setdiff(1:endSemidiurnal, symmAnn);
dTCoeffs(endSemidiurnal+1:end) = 0;
dTCoeffs(symmAnn) = 0;
quietInd = 1:length(dTCoeffs);
isGradient = 1;
significance = 2.0/3.0;

[dTCoeffs, JTWJ] = lbFit_mex(lbDTStruct, weights, dTCoeffs, paramsToFit, tolX, tolFun, tolOpt, lambda0, isGradient);
paramErrors = sqrt(abs(diag(inv(JTWJ)))); % POISTA ABS lopuillisessa.TESTAUS

lbDTStruct.coeffInd = 1:length(dTCoeffs);
lbDTStruct.numBiases = 0;
pe = ones(size(dTCoeffs));
pe(paramsToFit) = paramErrors;
paramsToFit = [];
[dTCoeffs, paramsToFit] = zeroOutInsignificantQuiet(dTCoeffs, paramsToFit, quietInd, pe, significance, lbDTStruct);

tolFun = 1E-10;
[dTCoeffs, JTWJ] = lbFit_mex(lbDTStruct, weights, dTCoeffs, paramsToFit, tolX, tolFun, tolOpt, lambda0, isGradient);

end

function [lbT0Coeffs, JTWJ] = fitLbTemerature(lbT0Struct, subsampPercent)

fprintf('%s\n', 'Fitting lower boundary temperature')

subsamp = round(1 / (subsampPercent / 100));

N = length(lbT0Struct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
lbT0Struct = removeDataPoints(lbT0Struct, removeInd);

Nobs = length(lbT0Struct.data);

lbT0Struct.data = [lbT0Struct.data; lbT0Struct.data];
lbT0Struct.timestamps = [lbT0Struct.timestamps; lbT0Struct.timestamps + 180];
lbT0Struct.latitude = [lbT0Struct.latitude; -lbT0Struct.latitude];
lbT0Struct.longitude = [lbT0Struct.longitude; lbT0Struct.longitude];
lbT0Struct.solarTime = [lbT0Struct.solarTime; lbT0Struct.solarTime];
lbT0Struct.altitude = [lbT0Struct.altitude; lbT0Struct.altitude];
lbT0Struct.aeInt = [lbT0Struct.aeInt; lbT0Struct.aeInt];
lbT0Struct.F = [lbT0Struct.F; lbT0Struct.F];
lbT0Struct.FA = [lbT0Struct.FA; lbT0Struct.FA];
lbT0Struct.apNow = [lbT0Struct.apNow; lbT0Struct.apNow];
lbT0Struct.ap3h = [lbT0Struct.ap3h; lbT0Struct.ap3h];
lbT0Struct.ap6h = [lbT0Struct.ap6h; lbT0Struct.ap6h];
lbT0Struct.ap9h = [lbT0Struct.ap9h; lbT0Struct.ap9h];
lbT0Struct.ap12To33h = [lbT0Struct.ap12To33h; lbT0Struct.ap12To33h];
lbT0Struct.ap36To57h = [lbT0Struct.ap36To57h; lbT0Struct.ap36To57h];
lbT0Struct.Ap = [lbT0Struct.Ap; lbT0Struct.Ap];
lbT0Struct.index = [lbT0Struct.index; lbT0Struct.index];

lbT0Struct = computeVariablesForFit(lbT0Struct);

latitude = zeros(1,14);
solarActivity = zeros(1,5);
annual = zeros(1,32); annual([1,7,12,16,19,25,30]) = 0.001;
diurnal = zeros(1, 21); diurnal([8, 19]) = -0.001;
semidiurnal = zeros(1,16);
terdiurnal = zeros(1,8);
quaterdiurnal = zeros(1,2);
longitudinal = zeros(1,13); longitudinal([2,5,9,12]) = 1E-4;
lbT0Coeffs = [zeros(1,1), latitude, solarActivity, annual, diurnal, semidiurnal, ...
    terdiurnal, quaterdiurnal, longitudinal];
lbT0Coeffs(1) = mean(lbT0Struct.data);

opt = optimoptions('lsqnonlin', 'Jacobian', 'on', 'Algorithm', 'trust-region-reflective', 'TolFun', 1E-8, ...
                 'TolX', 1E-7, 'Display', 'iter', 'MaxIter', 10000);

fun = @(X) lbTemperatureMinimization(lbT0Struct, X);

%[x,J] = fun(lbT0Coeffs);

N_santin = sum(lbT0Struct.index == 2);
weights = ones(size(lbT0Struct.data));
i = lbT0Struct.index == 1;
weights(i) = N_santin / sum(i);
i = lbT0Struct.index == 3;
weights(i) = N_santin / sum(i);
weights = sqrt(weights);

tolX = 1E-8;
tolFun = 1E-5;
tolOpt = 1E-5;
lambda0 = 1E-2;
endSemidiurnal = 89;
lat = 6:15;
solar = 16:20;
symmAnn = [25:26, 30:38];
assAnn = [43:44, 48:52];
diurnal = [56:59, 67:70];
semidiurnal = [76:79, 84:87];
paramsToFit = setdiff(1:endSemidiurnal, [lat, solar, symmAnn, assAnn, diurnal, semidiurnal]);
lbT0Coeffs(endSemidiurnal+1:end) = 0;
lbT0Coeffs([lat, solar, symmAnn, assAnn, diurnal, semidiurnal]) = 0;
quietInd = 1:length(lbT0Coeffs);
isGradient = 0;
significance = 2.0/3.0;

[lbT0Coeffs, JTWJ] = lbFit_mex(lbT0Struct, weights, lbT0Coeffs, paramsToFit, tolX, tolFun, tolOpt, lambda0, isGradient);
paramErrors = sqrt(abs(diag(inv(JTWJ))));

lbT0Struct.coeffInd = 1:length(lbT0Coeffs);
lbT0Struct.numBiases = 0;
pe = ones(size(lbT0Coeffs));
pe(paramsToFit) = paramErrors;
paramsToFit = [];
[lbT0Coeffs, paramsToFit] = zeroOutInsignificantQuiet(lbT0Coeffs, paramsToFit, quietInd, pe, significance, lbT0Struct, true);

tolFun = 1E-10;
[lbT0Coeffs, JTWJ] = lbFit_mex(lbT0Struct, weights, lbT0Coeffs, paramsToFit, tolX, tolFun, tolOpt, lambda0, isGradient);


end

% Function, whose root gives the Tex.
function y = TexFunc(Tex, Z, Tz, T0, dT0)
    %global modelLbHeight
    y = abs(Tex - (Tex - T0).*exp(-Z .* dT0 ./ (Tex-T0)) - Tz);
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
annual = ones(1,32); annual([1,7,12,16,19,25,30]) = 0.001;
diurnal = ones(1, 21); diurnal([8, 19]) = 0.003;
semidiurnal = ones(1,16);
terdiurnal = ones(1,8);
quaterdiurnal = ones(1,2);
longitudinal = ones(1,13); longitudinal([2,5,6,9,12,13]) = 1E-4;
geomagnetic = ones(1,29); geomagnetic([2,6,10,13,16,22,25]) = 1E-4;

ub = [latitude, solarActivity, annual, diurnal, semidiurnal, terdiurnal, quaterdiurnal, longitudinal, geomagnetic];
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

function [residual] = computeSpeciesResidual_major(varStruct, Tex, dT0, T0, coeff)

varStruct = computeDensityRHS(varStruct, Tex, dT0, T0);
Gvec = G_majorTex(coeff, varStruct, varStruct.numBiases);

if varStruct.numBiases == 0
    residual = ((varStruct.rhs) ./ (max(coeff(1) + Gvec, 1))) - 1;
elseif varStruct.numBiases > 0
    residual = ((varStruct.rhs) ./ (max(coeff(1) + sum(bsxfun(@times, coeff(2:varStruct.numBiases+1), varStruct.biases), 2) + Gvec, 1))) - 1;
else
    error('Incorrect number of biases!')
end

end

function [residual] = computeSpeciesResidual_minor(varStruct, Tex, dT0, T0, coeff)

varStruct = computeDensityRHS(varStruct, Tex, dT0, T0);
Gvec = G_minor(coeff, varStruct, varStruct.numBiases);

if varStruct.numBiases == 0
    residual = (varStruct.rhs ./ max(coeff(1) + Gvec, 1)) - 1;
elseif varStruct.numBiases > 0
    residual = (varStruct.rhs ./ max(coeff(1) + sum(bsxfun(@times, coeff(2:varStruct.numBiases+1), varStruct.biases), 2) + Gvec, 1)) - 1;
else
    error('Incorrect number of biases!')
end

end

function [residual] = computeSpeciesResidual_O2(varStruct, Tex, dT0, T0, coeff)

varStruct = computeDensityRHS(varStruct, Tex, dT0, T0);
averVal = coeff(1);

residual = ((varStruct.rhs) ./ (max(averVal, 1))) - 1;

end

function [residual, Jacobian] = modelMinimizationFunction(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, coeff, paramsToFit)

dataLen = length(TexStruct.data) + length(OStruct.data) + length(N2Struct.data) + length(HeStruct.data) + length(rhoStruct.data) ...
    + length(ArStruct.data) + length(O2Struct.data);
residual = zeros(dataLen, 1);

T0 = evalT0(TexStruct, T0Coeffs);
TexMesuredEstimate = clamp(T0+1, evalTex(TexStruct, coeff(TexStruct.coeffInd)), 5000);
residInd = 1:length(TexStruct.data);
residual(residInd) = TexStruct.data./TexMesuredEstimate - 1;

[Tex, dT0, T0] = findTempsForFit(OStruct, TexStruct, dTCoeffs, T0Coeffs, coeff);
residInd = residInd(end) + (1:length(OStruct.data));
residual(residInd) = computeSpeciesResidual_major(OStruct, Tex, dT0, T0, coeff(OStruct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(N2Struct, TexStruct, dTCoeffs, T0Coeffs, coeff);
residInd = residInd(end) + (1:length(N2Struct.data));
residual(residInd) = computeSpeciesResidual_major(N2Struct, Tex, dT0, T0, coeff(N2Struct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(HeStruct, TexStruct, dTCoeffs, T0Coeffs, coeff);
residInd = residInd(end) + (1:length(HeStruct.data));
residual(residInd) = computeSpeciesResidual_major(HeStruct, Tex, dT0, T0, coeff(HeStruct.coeffInd));

%[Tex, dT0, T0] = findTempsForFit(ArStruct, TexStruct, dTCoeffs, T0Coeffs, coeff);
%residInd = residInd(end) + (1:length(ArStruct.data));
%residual(residInd) = computeSpeciesResidual_major(ArStruct, Tex, dT0, T0, coeff(ArStruct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(O2Struct, TexStruct, dTCoeffs, T0Coeffs, coeff);
residInd = residInd(end) + (1:length(O2Struct.data));
residual(residInd) = computeSpeciesResidual_O2(O2Struct, Tex, dT0, T0, coeff(O2Struct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(rhoStruct, TexStruct, dTCoeffs, T0Coeffs, coeff);
residInd = residInd(end) + (1:length(rhoStruct.data));
OlbDens = clamp(10, evalMajorSpecies(rhoStruct, coeff(OStruct.coeffInd), OStruct.numBiases), 1E20);
N2lbDens = clamp(10, evalMajorSpecies(rhoStruct, coeff(N2Struct.coeffInd), N2Struct.numBiases), 1E20);
HelbDens = clamp(10, evalMajorSpecies(rhoStruct, coeff(HeStruct.coeffInd), HeStruct.numBiases), 1E20);
%ArlbDens = clamp(10, evalMajorSpecies(rhoStruct, coeff(ArStruct.coeffInd), ArStruct.numBiases), 1E20);
O2lbDens = clamp(10, exp(coeff(O2Struct.coeffInd)), 1E20);
modelRho = clamp(1E-20, computeRho(T0, dT0, Tex, rhoStruct.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens), 0.1);
residual(residInd) = ((rhoStruct.data)./(modelRho)) - 1;%(log(rhoStruct.data)./log(modelRho)) - 1;

if any(abs(residual) > 100)
    a=1;
end

residual = weights .* residual;

if nargout == 2
    fun = @(X)modelMinimizationFunction(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, X, paramsToFit);
    Jacobian = computeJAC(fun, coeff, dataLen, tolX, paramsToFit);
end

end

function [residual, Jacobian] = modelMinimizationFunction_log(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, allCoeff, paramsToFit, stormCoeff)

allCoeff(paramsToFit) = stormCoeff;

dataLen = length(TexStruct.data) + length(OStruct.data) + length(N2Struct.data) + length(HeStruct.data) + length(rhoStruct.data) ...
    + length(ArStruct.data) + length(O2Struct.data);
residual = zeros(dataLen, 1);

T0 = evalT0(TexStruct, T0Coeffs);
TexMesuredEstimate = clamp(T0+1, evalTex(TexStruct, allCoeff(TexStruct.coeffInd)), 5000);
residInd = 1:length(TexStruct.data);
residual(residInd) = TexStruct.data./TexMesuredEstimate - 1;

[Tex, dT0, T0] = findTempsForFit(OStruct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
residInd = residInd(end) + (1:length(OStruct.data));
residual(residInd) = computeSpeciesResidual_major(OStruct, Tex, dT0, T0, allCoeff(OStruct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(N2Struct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
residInd = residInd(end) + (1:length(N2Struct.data));
residual(residInd) = computeSpeciesResidual_major(N2Struct, Tex, dT0, T0, allCoeff(N2Struct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(HeStruct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
residInd = residInd(end) + (1:length(HeStruct.data));
residual(residInd) = computeSpeciesResidual_major(HeStruct, Tex, dT0, T0, allCoeff(HeStruct.coeffInd));

%[Tex, dT0, T0] = findTempsForFit(ArStruct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
%residInd = residInd(end) + (1:length(ArStruct.data));
%residual(residInd) = computeSpeciesResidual_major(ArStruct, Tex, dT0, T0, allCoeff(ArStruct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(O2Struct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
residInd = residInd(end) + (1:length(O2Struct.data));
residual(residInd) = computeSpeciesResidual_O2(O2Struct, Tex, dT0, T0, allCoeff(O2Struct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(rhoStruct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
residInd = residInd(end) + (1:length(rhoStruct.data));
OlbDens = clamp(10, evalMajorSpecies(rhoStruct, allCoeff(OStruct.coeffInd), OStruct.numBiases), 1E20);
N2lbDens = clamp(10, evalMajorSpecies(rhoStruct, allCoeff(N2Struct.coeffInd), N2Struct.numBiases), 1E20);
HelbDens = clamp(10, evalMajorSpecies(rhoStruct, allCoeff(HeStruct.coeffInd), HeStruct.numBiases), 1E20);
%ArlbDens = clamp(10, evalMajorSpecies(rhoStruct, allCoeff(ArStruct.coeffInd), ArStruct.numBiases), 1E20);
O2lbDens = clamp(10, exp(allCoeff(O2Struct.coeffInd)), 1E20);
modelRho = clamp(1E-20, computeRho(T0, dT0, Tex, rhoStruct.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens), 0.1);
residual(residInd) = (log(rhoStruct.data)./log(modelRho)) - 1;%(log(rhoStruct.data)./log(modelRho)) - 1;

if any(abs(residual) > 100)
    a=1;
end

residual = weights .* residual;

if nargout == 2
    fun = @(X)modelMinimizationFunction(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, X, paramsToFit);
    Jacobian = computeJAC(fun, allCoeff, dataLen, tolX, paramsToFit);
end

end

function [residual, Jacobian] = modelMinimizationFunction_lin(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, allCoeff, paramsToFit, stormCoeff)

allCoeff(paramsToFit) = stormCoeff;
if allCoeff(1) < 100
    allCoeff(1) = allCoeff(1)*100;
end

dataLen = length(TexStruct.data) + length(OStruct.data) + length(N2Struct.data) + length(HeStruct.data) + length(rhoStruct.data) ...
    + length(ArStruct.data) + length(O2Struct.data);
residual = zeros(dataLen, 1);

T0 = evalT0(TexStruct, T0Coeffs);
TexMesuredEstimate = clamp(T0+1, evalTex(TexStruct, allCoeff(TexStruct.coeffInd)), 5000);
residInd = 1:length(TexStruct.data);
residual(residInd) = TexStruct.data./TexMesuredEstimate - 1;

[Tex, dT0, T0] = findTempsForFit(OStruct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
residInd = residInd(end) + (1:length(OStruct.data));
residual(residInd) = computeSpeciesResidual_major(OStruct, Tex, dT0, T0, allCoeff(OStruct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(N2Struct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
residInd = residInd(end) + (1:length(N2Struct.data));
residual(residInd) = computeSpeciesResidual_major(N2Struct, Tex, dT0, T0, allCoeff(N2Struct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(HeStruct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
residInd = residInd(end) + (1:length(HeStruct.data));
residual(residInd) = computeSpeciesResidual_major(HeStruct, Tex, dT0, T0, allCoeff(HeStruct.coeffInd));

%[Tex, dT0, T0] = findTempsForFit(ArStruct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
%residInd = residInd(end) + (1:length(ArStruct.data));
%residual(residInd) = computeSpeciesResidual_major(ArStruct, Tex, dT0, T0, allCoeff(ArStruct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(O2Struct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
residInd = residInd(end) + (1:length(O2Struct.data));
residual(residInd) = computeSpeciesResidual_O2(O2Struct, Tex, dT0, T0, allCoeff(O2Struct.coeffInd));

[Tex, dT0, T0] = findTempsForFit(rhoStruct, TexStruct, dTCoeffs, T0Coeffs, allCoeff);
residInd = residInd(end) + (1:length(rhoStruct.data));
OlbDens = clamp(10, evalMajorSpecies(rhoStruct, allCoeff(OStruct.coeffInd), OStruct.numBiases), 1E20);
N2lbDens = clamp(10, evalMajorSpecies(rhoStruct, allCoeff(N2Struct.coeffInd), N2Struct.numBiases), 1E20);
HelbDens = clamp(10, evalMajorSpecies(rhoStruct, allCoeff(HeStruct.coeffInd), HeStruct.numBiases), 1E20);
%ArlbDens = clamp(10, evalMajorSpecies(rhoStruct, allCoeff(ArStruct.coeffInd), ArStruct.numBiases), 1E20);
O2lbDens = clamp(10, exp(allCoeff(O2Struct.coeffInd)), 1E20);
modelRho = clamp(1E-20, computeRho(T0, dT0, Tex, rhoStruct.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens), 0.1);
residual(residInd) = ((rhoStruct.data)./(modelRho)) - 1;%(log(rhoStruct.data)./log(modelRho)) - 1;

if any(abs(residual) > 100)
    a=1;
end

residual = weights .* residual;

residual = residual / 1;

if nargout == 2
    fun = @(X)modelMinimizationFunction_lin(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, allCoeff, paramsToFit, X);
    Jacobian = computeJAC(fun, stormCoeff, dataLen, tolX, 1:length(stormCoeff));
end

end



function residual = fun_efold(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, allCoeff, stormParams, numStorm, opt, efolds)

paramsToFit = stormParams;
paramsToFit(1:numStorm:end) = [];
allCoeff(stormParams(1:numStorm:end)) = efolds;
fun = @(coeff)modelMinimizationFunction_lin(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, allCoeff, paramsToFit, coeff);

lb = repmat([zeros(numStorm-1,1)-1E3],4,1); ub = repmat([zeros(numStorm-1,1)+1E3],4,1);
[optCoeff,~,funVec,~,output,~,JAC] = lsqnonlin(fun,allCoeff(paramsToFit),lb,ub,opt);

residual = sum(funVec.^2);

end

function residual = fun_bias(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, allCoeff, allInd, biasInd, opt, biases)

meanInd = setdiff(1:length(allInd), biasInd);
paramsToFit = allInd(meanInd);
paramsToFit((biasInd)) = [];
allCoeff(allInd(biasInd)) = biases;
fun = @(coeff)modelMinimizationFunction_lin(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, allCoeff, paramsToFit, coeff);

%lb = repmat([zeros(meanInd-1,1)-1E3],4,1); ub = repmat([zeros(meanInd-1,1)+1E3],4,1);
[optCoeff,~,funVec,~,output,~,JAC] = lsqnonlin(fun,allCoeff(paramsToFit),[],[],opt);

residual = sum(funVec.^2);

end

function [residual, Jacobian] = temperatureGradientMinimization(lbDTStruct, coeff)

modelDT = evalDT(lbDTStruct, coeff);
residual = lbDTStruct.data ./ modelDT - 1;

if nargout == 2
    fun = @(X)temperatureGradientMinimization(lbDTStruct, X);
    Jacobian = computeJAC(fun, coeff, length(lbDTStruct.data), 1E-5, 1:length(coeff));
end

end

function [residual, Jacobian] = lbTemperatureMinimization(lbT0Struct, coeff)

modelT0 = evalT0(lbT0Struct, coeff);
residual = (lbT0Struct.data ./ modelT0 - 1);

if nargout == 2
    fun = @(X)lbTemperatureMinimization(lbT0Struct, X);
    Jacobian = computeJAC(fun, coeff, length(lbT0Struct.data), 1E-5, 1:length(coeff));
end

end

function JAC = computeJAC(fun, x, dataLen, tolX, paramsToFit)

JAC = zeros(dataLen, length(paramsToFit));
dx = max(0.25*tolX*abs(x), 1E-10);

for i = 1:length(paramsToFit)
    k = paramsToFit(i);
    sumResult = 0;
    while sumResult == 0
        xForw = x; xBackw = x;        
        xForw(k) = x(k) + dx(k);
        xBackw(k) = x(k) - dx(k);

        F_forw = feval(fun, xForw);
        F_backw = feval(fun, xBackw);

        result = (F_forw - F_backw) / (2*dx(k));

        infInd = ~isfinite(result);
        if any(infInd)
            result(infInd) = mean(result(~infInd));
        end
        
        sumResult = 1;%sum(abs(result));
        if sumResult == 0
            dx(k) = dx(k) * 1E1;
            if (abs(x(k)) == 0 && dx(k) > 10) || (abs(x(k)) > 0 && dx(k) > abs(x(k)))
                error(['Unable to compute derivative for variable ',num2str(k)])
            end
        end
    end
    
    JAC(:,i) = result;
end

end

function [ind] = ignoreAnnual(varStruct, numQuietCoeffs)

ind = [1 : 20+varStruct.numBiases, (53 : numQuietCoeffs)+varStruct.numBiases];

end

function [ind] = ignoreAnnualAndLongitude(varStruct, numQuietCoeffs)

ind = [1 : 20+varStruct.numBiases, (53 : numQuietCoeffs-13)+varStruct.numBiases];

end

function [ind] = quietParams(var, numOtherCoeffs, trueInd)

if isnumeric(var)
    ind = 1:(length(var) - numOtherCoeffs);
else
    ind = 1 : numOtherCoeffs+var.numBiases;
end
if nargin == 3 && trueInd
    if isnumeric(var)
        ind = ind + var(1) - 1;
    else
        ind = ind + var.coeffInd(1) - 1;
    end
end

end

function [ind] = geomParams(varStruct, numQuietCoeffs)

ind = (numQuietCoeffs + varStruct.numBiases + 1) : length(varStruct.coeffInd);

end

function [ind] = ArParams(varStruct, numQuietCoeffs)

lat = 6:15;
symmAnn = [25:26, 30:38];
assAnn = [43:44, 48:52];
diur = [56:59, 67:70];
semiDiurn = [76:69, 84:87];
rest = 90:numQuietCoeffs;
removeInd = [lat, symmAnn, assAnn, diur, semiDiurn, rest] + varStruct.numBiases;
ind = setdiff(1 : numQuietCoeffs+varStruct.numBiases, removeInd);

end

function [ind] = HeParams(varStruct, numQuietCoeffs)

ind = ArParams(varStruct, numQuietCoeffs);

end

function [ind] = N2Params(varStruct, numQuietCoeffs)

lat = 10:15;
symmAnn = [32:38];
assAnn = [50:52];
rest = 90:numQuietCoeffs;
removeInd = [lat, symmAnn, assAnn, rest] + varStruct.numBiases;
ind = setdiff(1 : numQuietCoeffs+varStruct.numBiases, removeInd);

end

function [ind] = OParams(varStruct, numQuietCoeffs)

ind = quietParams(varStruct, numQuietCoeffs);

end

function [ind] = TexParams(varStruct, numQuietCoeffs)

ind = quietParams(varStruct, numQuietCoeffs);

end

function [] = fitModelVariables(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, options, multiStartSolver, numStartPoints, numThreads, quietData, fitBaseAgain, fitSimultaneously, quietCoeffs)
global numCoeffs;
numMinorCoeffs = 50;
numQuietCoeffs = 112;
numStormPrevious = numCoeffs - numQuietCoeffs;

removeInd = rhoStruct.swarm;
rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, false, true);
%if ~quietData
    removeInd = true(size(ArStruct.data)); %removeInd(1) = false;
    ArStruct = removeDataPoints(ArStruct, removeInd, true, true, true, true);
%end

% DEBUG
%removeInd = true(size(rhoStruct.data)); removeInd(1) = false;
%rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, false, true);

fprintf('%s\n', 'Computing final fit')

dataLen = length(TexStruct.data) + length(OStruct.data) + ...
length(N2Struct.data) + length(HeStruct.data) + length(ArStruct.data) +...
length(O2Struct.data) + length(rhoStruct.data);

TexStruct = computeVariablesForFit(TexStruct);
OStruct = computeVariablesForFit(OStruct);
N2Struct = computeVariablesForFit(N2Struct);
HeStruct = computeVariablesForFit(HeStruct);
%ArStruct = computeVariablesForFit(ArStruct);
O2Struct = computeVariablesForFit(O2Struct);
rhoStruct = computeVariablesForFit(rhoStruct);

if quietData
    tempSpecRelWeight = 0.25;
else
    tempSpecRelWeight = 0.05;
end

weights = computeWeights(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, tempSpecRelWeight); % 

%weights = ones(dataLen,1);

[G_lb, G_ub] = G_bounds();
if length(G_lb) + 1 ~= numCoeffs;
    error('length(G_lb) + 1 ~= numCoeffs');
end
%G_lb = 4 * G_lb; G_ub = 4 * G_ub;

TexStruct.coeffInd = 1:numCoeffs;
lb = [500, G_lb]; ub = [1500, G_ub];

OStruct.coeffInd = TexStruct.coeffInd(end) + (1:numCoeffs+OStruct.numBiases);
lb = [lb, log(0.5E10), zeros(1, OStruct.numBiases)-0, G_lb]; % MUISTA LISATA BIASET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ub = [ub, log(1E11), zeros(1, OStruct.numBiases)-0, G_ub];

N2Struct.coeffInd = OStruct.coeffInd(end) + (1:numCoeffs+N2Struct.numBiases);
lb = [lb, log(0.5E11), zeros(1, N2Struct.numBiases)-0, G_lb]; % MUISTA LISATA BIASET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ub = [ub, log(1E12), zeros(1, N2Struct.numBiases)-0, G_ub];

HeStruct.coeffInd = N2Struct.coeffInd(end) + (1:numCoeffs+HeStruct.numBiases);
lb = [lb, log(0.5E7), zeros(1, HeStruct.numBiases)-0, G_lb]; % MUISTA LISATA BIASET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ub = [ub, log(1E8), zeros(1, HeStruct.numBiases)-0, G_ub];

%ArStruct.coeffInd = HeStruct.coeffInd(end) + (1:numCoeffs+ArStruct.numBiases);
%lb = [lb, log(0.5E9), zeros(1, ArStruct.numBiases)-0, G_lb]; % MUISTA LISATA BIASET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%ub = [ub, log(2E9), zeros(1, ArStruct.numBiases)-0, G_ub];

%ArStruct.coeffInd = HeStruct.coeffInd(end) + (1:numCoeffs+ArStruct.numBiases);
%O2Struct.coeffInd = ArStruct.coeffInd(end) + 1;
O2Struct.coeffInd = HeStruct.coeffInd(end) + 1;

%startPoints = createRandomStartPoints(lb, ub, numStartPoints);
%initGuess = list(startPoints);
ind = ub < 0.5 | ub > 2;
%ind(1:numCoeffs) = ub(1:numCoeffs) < 100;
initGuess = mean([lb;ub]);% - 0.001;
initGuess(ind) = ub(ind);
%ArCoeffs = zeros(numCoeffs+ArStruct.numBiases, 1);
%ArCoeffs([17, 24, 28, 32, 35, 39, 41, 45, 46, 50] + ArStruct.numBiases) = 0.001;
%initGuess(ArStruct.coeffInd) = ArCoeffs;

initGuess(TexStruct.coeffInd(1)) = 950;
initGuess(OStruct.coeffInd(1)) = log(4E10);
initGuess(N2Struct.coeffInd(1)) = log(1.4E11);
initGuess(HeStruct.coeffInd(1)) = log(2.2E7);
%initGuess(ArStruct.coeffInd(1)) = log(0.4E9);
initGuess(O2Struct.coeffInd) = log(1.4E10);

quietInd = [TexStruct.coeffInd(TexParams(TexStruct, numQuietCoeffs)),...
            OStruct.coeffInd(OParams(OStruct, numQuietCoeffs)),...
            N2Struct.coeffInd(N2Params(N2Struct, numQuietCoeffs)),...
            HeStruct.coeffInd(HeParams(HeStruct, numQuietCoeffs)),...
            O2Struct.coeffInd(1)]; 

stormInd = [TexStruct.coeffInd(geomParams(TexStruct, numQuietCoeffs)),...
            OStruct.coeffInd(geomParams(OStruct, numQuietCoeffs)),...
            N2Struct.coeffInd(geomParams(N2Struct, numQuietCoeffs)),...
            HeStruct.coeffInd(geomParams(HeStruct, numQuietCoeffs))]; 
%initGuess(ArStruct.coeffInd(geomParams(ArStruct, numQuietCoeffs))) = 0;

tolX = 1E-8;
tolFun = 1E-6;
tolOpt = 1E-4;
lambda0 = 1E0;
if quietData
    minLambda = 1E-10;
else
    minLambda = 1E3;
end
% if quietData && ~fitSimultaneously
%     paramsToFit = [TexStruct.coeffInd(1),...
%             OStruct.coeffInd(1:1+OStruct.numBiases),...
%             N2Struct.coeffInd(1:1+N2Struct.numBiases),...
%             HeStruct.coeffInd(1:1+HeStruct.numBiases),...
%             ArStruct.coeffInd(1:1+ArStruct.numBiases),...
%             O2Struct.coeffInd(1)];
%     setenv('OMP_NUM_THREADS', num2str(numThreads))
%     disp('Calling LM solver')
%     clear mex;
%     meanGuess = initGuess;
%     meanGuess(setdiff(1:length(meanGuess), paramsToFit)) = 0;
%     tic;[optCoeff, JTWJ] = levenbergMarquardt_mex(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, meanGuess, paramsToFit, tolX, tolFun, tolOpt, lambda0);toc;
%     fprintf('Mean parameters refitted.\n');
%     
%     initGuess(paramsToFit) = optCoeff(paramsToFit);
% end

if ~fitSimultaneously
    if quietData
        paramsToFit = quietInd;
        otherInd = setdiff(1:length(initGuess), quietInd);
        initGuess(otherInd) = 0;
%         load quietCoeffsAll.F30.mat
%         numStormPrevious = TexInd(length(TexInd)) - numQuietCoeffs;
%         initGuess(quietParams(TexStruct,numQuietCoeffs,true)) = optCoeff(quietParams(TexInd,numStormPrevious,true));
%         initGuess(quietParams(OStruct,numQuietCoeffs,true)) = optCoeff(quietParams(OInd,numStormPrevious,true));
%         initGuess(quietParams(N2Struct,numQuietCoeffs,true)) = optCoeff(quietParams(N2Ind,numStormPrevious,true));
%         initGuess(quietParams(HeStruct,numQuietCoeffs,true)) = optCoeff(quietParams(HeInd,numStormPrevious,true));
%         initGuess(quietParams(ArStruct,numQuietCoeffs,true)) = optCoeff(quietParams(ArInd,numStormPrevious,true));
%         initGuess(O2Struct.coeffInd(1)) = optCoeff(O2Ind);
    else
        paramsToFit = stormInd;
        load quietCoeffs.mat
        %initGuess(quietInd) = quietCoeffs(quietInd);
        numStormPrevious = TexInd(length(TexInd)) - numQuietCoeffs;
        initGuess(quietParams(TexStruct,numQuietCoeffs,true)) = quietCoeffs(quietParams(TexInd,numStormPrevious,true));
        initGuess(quietParams(OStruct,numQuietCoeffs,true)) = quietCoeffs(quietParams(OInd,numStormPrevious,true));
        initGuess(quietParams(N2Struct,numQuietCoeffs,true)) = quietCoeffs(quietParams(N2Ind,numStormPrevious,true));
        initGuess(quietParams(HeStruct,numQuietCoeffs,true)) = quietCoeffs(quietParams(HeInd,numStormPrevious,true));
        %initGuess(quietParams(ArStruct,numQuietCoeffs,true)) = quietCoeffs(quietParams(ArInd,numStormPrevious,true));
        initGuess(O2Struct.coeffInd(1)) = quietCoeffs(O2Ind);
    end
else
    paramsToFit = 1:length(initGuess);
end

fun = @(coeff)modelMinimizationFunction(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, coeff, paramsToFit);
[comp] = fun(initGuess);
if ~quietData
    fun_log = @(coeff)modelMinimizationFunction_log(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, initGuess, paramsToFit, coeff);
    %[comp_log,JAC] = fun_log(initGuess(paramsToFit));
    
    fun_lin = @(coeff)modelMinimizationFunction_lin(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, initGuess, paramsToFit, coeff);
    %[comp_lin,JAC] = fun_lin(initGuess(paramsToFit));
end
%error('Not complete');
[comp] = fun(initGuess);
%[derivNorms, indSort] = sort(rms(JAC)); indSort = paramsToFit(indSort);
%[minDiff, indMin] = min(diff(derivNorms));
%fprintf('Min difference between param effects: %e for indices (%d, %d)\n', minDiff, indSort(indMin), indSort(indMin+1))
%fprintf('Median difference: %e\n', median(diff(derivNorms)));

%JTJ_diag = diag(JAC'*JAC);

if ~fitSimultaneously
    if quietData
        filename = 'quietCoeffsAll.mat';
        tolFun = 1E-4;
        tolOpt = 1E0;
    else
        filename = 'stormCoeffsAll.mat';
        tolFun = 1E-5;
        tolOpt = 1E0;
    end
else
    filename = 'coeffsAll.mat';
end

if fitSimultaneously || fitBaseAgain
    setenv('OMP_NUM_THREADS', num2str(numThreads))
    disp('Calling LM solver')
    clear mex;
    if quietData
        
          Obiases = OStruct.coeffInd(2:1+OStruct.numBiases);
          N2biases = N2Struct.coeffInd(2:1+N2Struct.numBiases);
          HeBiases = HeStruct.coeffInd(2:1+HeStruct.numBiases);
          %ArBiases =  ArStruct.coeffInd(2:1+ArStruct.numBiases);
          initGuess(Obiases) = [0.0680	0.2379	0.0665	0.1053	-0.1162];
          initGuess(N2biases) = [-0.1043	0.1267	0.0401	-0.0289	0.1312	-0.1166];
          initGuess(HeBiases) = [-0.0223	-0.0086	0.0474	-0.0968	0.1847];
          %initGuess(ArBiases) = [-0.0043	0.0880];
          paramsToFit = setdiff(paramsToFit,[Obiases, N2biases, HeBiases]);
%         
%         rmInd = setdiff(1:length(initGuess), paramsToFit);
%         initGuess(rmInd) = 0;
%         biasInd = find(initGuess(paramsToFit) == 0); 
%         meanInd = setdiff(1:length(paramsToFit), biasInd);
%         
%         opt = optimoptions('lsqnonlin', 'Jacobian', 'off', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1E-5, ...
%                   'TolX', 1E-4, 'Display', 'off');
%         
%         %fun = @(X) fun_bias(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, allCoeff, paramsToFit, biasInd, opt, biases);
%         fun = @(coeff)sum(modelMinimizationFunction_lin(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, initGuess, paramsToFit, coeff).^2);
%         
%         initGuess(paramsToFit(1)) = 9.5;
%         lb = initGuess(paramsToFit); ub = lb;
%         lb(biasInd) = log(0.3); ub(biasInd) = log(1.5);
%         lb(meanInd) = [8.0, initGuess(paramsToFit(meanInd(2:end)))-3];
%         ub(meanInd) = [11.0, initGuess(paramsToFit(meanInd(2:end)))+3];
%         options = psoptimset('Display','iter','tolfun',1E-4,'tolmesh',1E-2,'useparallel',false,'CompletePoll','off','Vectorized','off',...
%             'maxiter',1000,'MaxMeshSize',1.0,'cache','on','cacheTol',1E-4);
%         p = gcp('nocreate');
%         numStartPoints = p.NumWorkers;
%         initPoints = repmat(initGuess(paramsToFit), numStartPoints, 1);
%         final_points = zeros(numStartPoints, length(paramsToFit));
%         final_fvals = zeros(numStartPoints, 1);
%         for i = 1:numStartPoints
%             initPoints(i,biasInd) = (lb(biasInd)) + rand(size(lb(biasInd))).*(ub(biasInd)-lb(biasInd));
%         end
%         
%         parfor i = 1:numStartPoints
%             [final_points(i,:),final_fvals(i)] = patternsearch(fun,initPoints(i,:),[],[],[],[],lb,ub,[],options);
%         end
%         save('biases.mat','final_points','final_fvals','initPoints');

        tic;[optCoeff, JTWJ] = levenbergMarquardt_mex(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, initGuess, paramsToFit, tolX, tolFun, tolOpt, lambda0, minLambda);toc;
    else
         numStorm = numCoeffs - numQuietCoeffs;
%         efolds_init = [7, 11.4, 13.0, 7.3];
%         opt = optimoptions('lsqnonlin', 'Jacobian', 'on', 'Algorithm', 'trust-region-reflective', 'TolFun', 1E-5, ...
%                  'TolX', 1E-4, 'Display', 'off');
%         fun = @(X) fun_efold(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, initGuess, paramsToFit, numStorm, opt, X);
%         %comp = fun(efolds_init);
%         %lb = efolds_init - 6; ub = efolds_init + 6; lb(3) = 4;
%         lb = 2*ones(size(efolds_init)); ub = 18*ones(size(efolds_init));
%         options = psoptimset('Display','iter','tolfun',1E-4,'tolmesh',1E-2,'useparallel',false,'CompletePoll','off','Vectorized','off','maxiter',1000);
%         %prevThreads = maxNumCompThreads(1);
%         
%         p = gcp('nocreate');
%         numStartPoints = p.NumWorkers;
%         initPoints = zeros(numStartPoints, length(efolds_init));
%         final_efolds = zeros(numStartPoints, length(efolds_init));
%         final_fvals = zeros(numStartPoints, 1);
%         for i = 1:numStartPoints
%             initPoints(i,:) = (lb+1) + rand(size(lb)).*(ub-lb-2);
%         end
%         parfor i = 1:numStartPoints
%             [final_efolds(i,:),final_fvals(i)] = patternsearch(fun,initPoints(i,:),[],[],[],[],lb,ub,[],options);
%         end
%         save('efolds.mat','final_efolds','final_fvals','initPoints');
        
        %maxNumCompThreads(prevThreads);

        initGuess(paramsToFit(1:numStorm:end)) = [6.0546   17.9177    5.2326    3.9629];
        paramsToFit(1:numStorm:end) = [];
        %lb = repmat([2; zeros(1,1)-1E3],4,1); ub = repmat([22; zeros(1,1)+1E3],4,1);
        fun_lin = @(coeff)modelMinimizationFunction_lin(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, initGuess, paramsToFit, coeff);
        %[comp,JAC] = fun_lin(initGuess(paramsToFit));
        tic;[optCoeff,~,funVec,~,output,~,JAC] = lsqnonlin(fun_lin,initGuess(paramsToFit),[],[],options);toc;
               
%         weights = computeWeights(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, 0.1);
%         fun_lin = @(coeff)modelMinimizationFunction_lin(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, tolX, initGuess, paramsToFit, coeff);
%         [optCoeff,~,funVec,~,output,~,JAC] = lsqnonlin(fun_lin,optCoeff,lb,ub,options);
        
        initGuess(paramsToFit) = optCoeff;
        optCoeff = initGuess;
        
        stdFit = funVec' * funVec / (length(funVec) - length(paramsToFit) + 1);
        JTWJ = JAC' * JAC / stdFit;
    end
    if quietData && all(optCoeff == initGuess')
        error('Cholesky failed?');
    end
    %[comp] = fun(initGuess); %disp([comp(1), optCoeff]);
    %JTJ_diag_matlab = diag(J'*J);

    saveToFile(filename, optCoeff, JTWJ, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, dTCoeffs, T0Coeffs)
    fprintf('All parameters refitted.\n');
else
    load(filename)
%     load onePercent_err % TESTAUS
%     load onePercent_coeff % TESTAUS
end

if quietData
    ind = paramsToFit;
else
    ind = paramsToFit;
end
i=find(ismember(ind,TexStruct.coeffInd)); paramErrors_Tex = sqrt(abs(diag(inv(JTWJ(i,i)))));
i=find(ismember(ind,OStruct.coeffInd)); paramErrors_O = sqrt(abs(diag(inv(JTWJ(i,i)))));
i=find(ismember(ind,N2Struct.coeffInd)); paramErrors_N2 = sqrt(abs(diag(inv(JTWJ(i,i)))));
i=find(ismember(ind,HeStruct.coeffInd)); paramErrors_He = sqrt(abs(diag(inv(JTWJ(i,i)))));
%i=find(ismember(ind,ArStruct.coeffInd)); paramErrors_Ar = sqrt(abs(diag(inv(JTWJ(i,i)))));
if quietData
    paramErrors_O2 = sqrt(1/abs(JTWJ(end,end)));
    paramErrors = [paramErrors_Tex; paramErrors_O; paramErrors_N2; paramErrors_He; ...
                paramErrors_O2];
else
    paramErrors = [paramErrors_Tex; paramErrors_O; paramErrors_N2; paramErrors_He;];
end

allInd = 1:length(optCoeff);
pe = ones(size(allInd)); 
if quietData
    pe(paramsToFit) = paramErrors;
else
    pe(paramsToFit) = paramErrors;
end
significance = 2.0/3.0;
if ~fitSimultaneously
    if quietData
        paramsToFit = [];
        [optCoeff, paramsToFit] = zeroOutInsignificantQuiet(optCoeff, paramsToFit, allInd, pe, significance, TexStruct);
        [optCoeff, paramsToFit] = zeroOutInsignificantQuiet(optCoeff, paramsToFit, allInd, pe, significance, OStruct);
        [optCoeff, paramsToFit] = zeroOutInsignificantQuiet(optCoeff, paramsToFit, allInd, pe, significance, N2Struct);
        [optCoeff, paramsToFit] = zeroOutInsignificantQuiet(optCoeff, paramsToFit, allInd, pe, significance, HeStruct);
        %[optCoeff, paramsToFit] = zeroOutInsignificantQuiet(optCoeff, paramsToFit, allInd, pe, significance, ArStruct);
        paramsToFit = [paramsToFit, O2Struct.coeffInd]; %optCoeff(O2Struct.coeffInd) = 20.0;
    else
        %[optCoeff, paramsToFit] = zeroOutInsignificantStorm(optCoeff, paramsToFit, stormInd, paramErrors, significance);TESTAUS
    end
else
    paramsToFit = [];
        [optCoeff, paramsToFit] = zeroOutInsignificantQuiet(optCoeff, paramsToFit, allInd, pe, significance, TexStruct);
        [optCoeff, paramsToFit] = zeroOutInsignificantQuiet(optCoeff, paramsToFit, allInd, pe, significance, OStruct);
        [optCoeff, paramsToFit] = zeroOutInsignificantQuiet(optCoeff, paramsToFit, allInd, pe, significance, N2Struct);
        [optCoeff, paramsToFit] = zeroOutInsignificantQuiet(optCoeff, paramsToFit, allInd, pe, significance, HeStruct);
        %[optCoeff, paramsToFit] = zeroOutInsignificantQuiet(optCoeff, paramsToFit, allInd, pe, significance, ArStruct);
        paramsToFit = [paramsToFit, O2Struct.coeffInd]; %optCoeff(O2Struct.coeffInd) = 20.0;
    %[optCoeff, paramsToFit] = zeroOutInsignificantStorm(optCoeff, paramsToFit, stormInd, paramErrors, significance);TESTAUS
    paramsToFit = sort([paramsToFit, stormInd]); % Testaus
end

tolFun = 1E-5;
tolOpt = 1E0;
lambda0 = 1E-2;

if fitSimultaneously || quietData % TESTAUS. Kunnes Myrsky-yhtalo saavuttanut loppulisen muotonsa ja zeroOutInsignificantStorm on koodattu
    setenv('OMP_NUM_THREADS', num2str(numThreads))
    disp('Calling LM solver')
    clear mex;
    tic;[optCoeff, JTWJ] = levenbergMarquardt_mex(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, dTCoeffs, T0Coeffs, weights, optCoeff, paramsToFit, tolX, tolFun, tolOpt, lambda0, minLambda);toc;
    fprintf('Significant parameters refitted.\n');
end

if iscolumn(optCoeff) optCoeff = optCoeff'; end
if quietData
    filename = 'quietCoeffs.mat';
else
    filename = 'optCoeff.mat';
end
saveToFile(filename, optCoeff, JTWJ, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, dTCoeffs, T0Coeffs)

end

function [] = saveToFile(filename, optCoeff, JTWJ, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, dTCoeffs, T0Coeffs)

save(filename, 'optCoeff', '-v7.3');
TexInd = TexStruct.coeffInd; save(filename, 'TexInd', '-append');
HeInd = HeStruct.coeffInd; save(filename, 'HeInd', '-append');
OInd = OStruct.coeffInd; save(filename, 'OInd', '-append');
N2Ind = N2Struct.coeffInd; save(filename, 'N2Ind', '-append');
%ArInd = ArStruct.coeffInd; save(filename, 'ArInd', '-append');
O2Ind = O2Struct.coeffInd; save(filename, 'O2Ind', '-append');
save(filename, 'JTWJ', '-append')
save(filename, 'dTCoeffs', '-append');
save(filename, 'T0Coeffs', '-append');

end

function weights = computeWeights(TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct, tempSpecRelWeight)

wTex = ones(size(TexStruct.data)); wTex(TexStruct.de2) = length(TexStruct.aeE) / length(TexStruct.de2);
wO = ones(size(OStruct.data)); wO(OStruct.de2) = length(OStruct.aeENace) / length(OStruct.de2); wO(OStruct.guvi) = length(OStruct.aeENace) / length(OStruct.guvi);
wN2 = ones(size(N2Struct.data)); wN2(N2Struct.de2) = length(N2Struct.aeENace) / length(N2Struct.de2); wN2(N2Struct.aeros) = 0.5*length(N2Struct.aeENace) / length(N2Struct.aeros); wN2(N2Struct.guvi) = length(N2Struct.aeENace) / length(N2Struct.guvi);
wHe = ones(size(HeStruct.data)); wHe(HeStruct.de2) = length(HeStruct.aeENace) / length(HeStruct.de2); wHe(HeStruct.aeros) = 0.5*length(HeStruct.aeENace) / length(HeStruct.aeros); 
%wAr = ones(size(ArStruct.data)); wAr(ArStruct.de2) = 2*length(ArStruct.aeros) / length(ArStruct.de2);
wO2 = ones(size(O2Struct.data));

%wAr = wAr / 4; 
wTex = wTex * 2; wN2 = wN2 / 2;
tempSpecWeight = [wTex; wO; wN2; wHe; wO2];

dataLen = length(TexStruct.data) + length(OStruct.data) + length(N2Struct.data) + length(HeStruct.data) + length(rhoStruct.data)...
    + length(ArStruct.data) + length(O2Struct.data);
TempAndSpectrometerLen = dataLen - length(rhoStruct.data);
weights = ones(dataLen, 1);

w = (tempSpecRelWeight / (1 - tempSpecRelWeight)) * sum(rhoStruct.weights) / sum(tempSpecWeight);
wInd = 1:TempAndSpectrometerLen;
weights(wInd) = tempSpecWeight * w;

weights(wInd(end)+1:end) = rhoStruct.weights;

% ae16h = [TexStruct.aeInt(:,4); OStruct.aeInt(:,4); N2Struct.aeInt(:,4); HeStruct.aeInt(:,4); ...
%     ArStruct.aeInt(:,4); O2Struct.aeInt(:,4); rhoStruct.aeInt(:,4)];
% aeThreshold = 250;
% ind = ae16h >= aeThreshold;
% w = sum(weights(~ind)) / sum(weights(ind));
% weights(ind) = w * weights(ind);

% goceInd = TempAndSpectrometerLen + rhoStruct.goce;
% graceInd = TempAndSpectrometerLen + rhoStruct.grace;
% swarmInd = TempAndSpectrometerLen + rhoStruct.swarm;
% wGoce = 0.125 * sum(weights(graceInd)) / sum(weights(goceInd));
% weights(goceInd) = wGoce * weights(goceInd);
% wSwarm = 0.5 * sum(weights(goceInd)) / length(swarmInd);
% weights(swarmInd) = wSwarm;

% aeNormalized = 1 + (2 * ae16h / max(ae16h));
% weights = weights .* aeNormalized;

wTempSpec = wInd;
save('weights.mat','weights')
save('weights.mat','wTempSpec','-append')

weights = sqrt(weights);

end

function stop = outfun(x, optimValues, state)

fprintf('Coefficient(1): %12.6f \n', x(1))
stop = false;

end
