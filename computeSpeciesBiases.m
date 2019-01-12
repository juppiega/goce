function [param, msisParam, dtmParam] = computeParams(S, coeffStruct, paramName, numBiasesStruct)

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