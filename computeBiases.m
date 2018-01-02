function [] = computeBiases (name, ind)

if strcmpi(name,'O')
    load('ilData.mat','OStruct');
    S = OStruct;
elseif strcmpi(name,'N2')
    load('ilData.mat','N2Struct');
    S = N2Struct;
elseif strcmpi(name,'He')
    load('ilData.mat','HeStruct');
    S = HeStruct;
elseif strcmpi(name,'Tex')
    load('ilData.mat','TexStruct');
    S = TexStruct;
end

rmInd = setdiff(1:length(S.data),ind);
S = removeDataPoints(S, rmInd,true,true,true,true);
S = computeVariablesForFit(S);

load optCoeff

if ~strcmpi(name, 'Tex') && ~strcmpi(name, 'T0') && ~strcmpi(name, 'dT')
    [Tex, dT0, T0] = findTempsForFit_this(S, optCoeff(TexInd), dTCoeffs, T0Coeffs);
    OlbDens = evalMajorSpecies(S, optCoeff(OInd), 5);
    N2lbDens = evalMajorSpecies(S, optCoeff(N2Ind), 6);
    HelbDens = evalMajorSpecies(S, optCoeff(HeInd), 5);
    O2lbDens = exp(optCoeff(O2Ind));
    S = computeDensityRHS(S, Tex, dT0, T0);
    data_il = S.rhs;
    
    Tex_msis = computeMsis(S);
    Tex_dtm = computeDtm(S);
    [T0_msis, dT0_msis, T0_dtm, dT0_dtm] = computeMsisDtmLb(S);
    S = computeDensityRHS(S, Tex_msis, dT0_msis, T0_msis);
    data_msis = S.rhs;
    
    S = computeDensityRHS(S, Tex_dtm, dT0_dtm, T0_dtm);
    data_dtm = S.rhs;
end

S.altitude(:) = 130; S.Z(:) = 0;

if strcmpi(name,'O')
    [~,~,model_dtm] = computeDtm(S);
    [~,~,model_msis] = computeMsis(S);
    [~,model_il] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);
elseif strcmpi(name,'N2')
    [~,~,~,model_dtm] = computeDtm(S);
    [~,~,~,model_msis] = computeMsis(S);
    [~,~,model_il] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);
elseif strcmpi(name,'He')
    [~,~,~,~,model_dtm] = computeDtm(S);
    [~,~,~,~,model_msis] = computeMsis(S);
    [~,~,~,model_il] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);
elseif strcmpi(name,'Tex')
    [model_dtm] = computeDtm(S);
    [model_msis] = computeMsis(S);
    model_il = evalTex(S, optCoeff(TexInd));
end

logOM_dtm = findBestBias(log(model_dtm), data_dtm)
logOM_msis = findBestBias(log(model_msis), data_msis)
logOM_il = findBestBias(log(model_il), data_il)


end

function bias = findBestBias(model, observed)

fun = @(x)mean(observed ./ (x + model)) - 1;
bias = fzero(fun, 0);

end

function [Tex, dT0, T0] = findTempsForFit_this(varStruct, TexCoeffs, dTCoeffs, T0Coeffs)

Tex_est = evalTex(varStruct, TexCoeffs);

T0 = clamp(200, evalT0(varStruct, T0Coeffs), 1000);
dT0 = clamp(1, evalDT(varStruct, dTCoeffs), 30);
Tex = clamp(T0+1, Tex_est, 5000);

end

function S = computeDensityRHS(S, Tex, dT0, T0)

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
