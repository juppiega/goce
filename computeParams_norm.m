function [param, msisParam, dtmParam, obs_130] = computeParams_norm(S, coeffStruct, paramName, numBiasesStruct, norm_130)

TexCoeffs = coeffStruct.TexCoeff; dTCoeffs = coeffStruct.dTCoeff;
T0Coeffs = coeffStruct.T0Coeff;


if ~strcmpi(paramName, 'Tex') && ~strcmpi(paramName, 'T0') && ~strcmpi(paramName, 'dT')
    [Tex, dT0, T0] = findTempsForFit_this(S, TexCoeffs, dTCoeffs, T0Coeffs);
    OlbDens = evalMajorSpecies(S, coeffStruct.OCoeff, numBiasesStruct.O);
    N2lbDens = evalMajorSpecies(S, coeffStruct.N2Coeff, numBiasesStruct.N2);
    HelbDens = evalMajorSpecies(S, coeffStruct.HeCoeff, numBiasesStruct.He);
    %ArlbDens = evalMajorSpecies(S, coeffStruct.ArCoeff, numBiasesStruct.Ar);
    O2lbDens = exp(coeffStruct.O2Coeff);
end

if nargin > 4 && norm_130
    S.name = paramName;
    S = computeDensityRHS(S, Tex, dT0, T0);
    obs_130 = exp(S.rhs);
    S.Z(:) = 0;
    S.altitude(:) = 130;
end

if strcmpi(paramName, 'Tex')
    param = evalTex(S, TexCoeffs);
    msisParam = computeMsis(S);
    dtmParam = computeDtm(S);
elseif strcmpi(paramName, 'T0')
    param = evalT0(S, T0Coeffs);
    [msisParam, ~, dtmParam, ~] = computeMsisDtmLb(S);
elseif strcmpi(paramName, 'dT')
    param = evalDT(S, dTCoeffs);
    [~, msisParam, ~, dtmParam] = computeMsisDtmLb(S);
elseif strcmpi(paramName, 'rho')
    param = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);
    [~,msisParam] = computeMsis(S);
    [~,dtmParam] = computeDtm(S);
elseif strcmpi(paramName, 'O')
    [~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);
    [~,~,msisParam] = computeMsis(S);
    [~,~,dtmParam] = computeDtm(S);
elseif strcmpi(paramName, 'N2')
    [~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);
    [~,~,~,msisParam] = computeMsis(S);
    [~,~,~,dtmParam] = computeDtm(S);
elseif strcmpi(paramName, 'O_N2')
    [~,~,il_N2] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);
    [~,~,~,msis_N2] = computeMsis(S);
    [~,~,~,dtm_N2] = computeDtm(S);
    [~,il_O] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);
    [~,~,msis_O] = computeMsis(S);
    [~,~,dtm_O] = computeDtm(S);
    param = il_O./il_N2;
    msisParam = msis_O./msis_N2;
    dtmParam = dtm_O./dtm_N2;
elseif strcmpi(paramName, 'He')
    [~,~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);   
    [~,~,~,~,msisParam] = computeMsis(S);
    [~,~,~,~,dtmParam] = computeDtm(S);
elseif strcmpi(paramName, 'Ar')
    [~,~,~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);   
    [~,~,~,~,~,msisParam] = computeMsis(S);
    dtmParam = zeros(size(msisParam));
elseif strcmpi(paramName, 'O2')
    [~,~,~,~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);   
    [~,~,~,~,~,~,msisParam] = computeMsis(S);
    dtmParam = zeros(size(msisParam));
else
    error(['Unknown variable: ', paramName])
end



end

function [Tex, dT0, T0] = findTempsForFit_this(varStruct, TexCoeffs, dTCoeffs, T0Coeffs, coeff)

Tex_est = evalTex(varStruct, TexCoeffs);

T0 = clamp(200, evalT0(varStruct, T0Coeffs), 1000);
dT0 = clamp(1, evalDT(varStruct, dTCoeffs), 30);
Tex = clamp(T0+1, Tex_est, 5000);

end

function S = computeDensityRHS(S, Tex, dT0, T0)
%global modelLbHeight

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