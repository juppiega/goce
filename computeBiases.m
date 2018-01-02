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
    %ArlbDens = evalMajorSpecies(S, coeffStruct.ArCoeff, numBiasesStruct.Ar);
    O2lbDens = exp(optCoeff(O2Ind));
end

if strcmpi(name,'O')
    [~,~,obs_dtm] = computeDtm(S);
    [~,~,obs_msis] = computeMsis(S);
    [~,obs_il] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);
elseif strcmpi(name,'N2')
    [~,~,~,obs_dtm] = computeDtm(S);
    [~,~,~,obs_msis] = computeMsis(S);
    [~,~,obs_il] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);
elseif strcmpi(name,'He')
    [~,~,~,~,obs_dtm] = computeDtm(S);
    [~,~,~,~,obs_msis] = computeMsis(S);
    [~,~,~,obs_il] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);
elseif strcmpi(name,'Tex')
    [obs_dtm] = computeDtm(S);
    [obs_msis] = computeMsis(S);
    obs_il = evalTex(S, optCoeff(TexInd));
end

logOM_dtm = mean(log(S.data./obs_dtm))
logOM_msis = mean(log(S.data./obs_msis))
logOM_il = mean(log(S.data./obs_il))

end

function [Tex, dT0, T0] = findTempsForFit_this(varStruct, TexCoeffs, dTCoeffs, T0Coeffs)

Tex_est = evalTex(varStruct, TexCoeffs);

T0 = clamp(200, evalT0(varStruct, T0Coeffs), 1000);
dT0 = clamp(1, evalDT(varStruct, dTCoeffs), 30);
Tex = clamp(T0+1, Tex_est, 5000);

end
