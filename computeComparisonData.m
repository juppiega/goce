function [ilRho, msisRho, dtmRho] = computeComparisonData(rhoStruct, coeffStruct, numBiasesStruct)

rhoStruct = computeVariablesForFit(rhoStruct);
TexCoeffs = coeffStruct.TexCoeff; dTCoeffs = coeffStruct.dTCoeff;
T0Coeffs = coeffStruct.T0Coeff;
[Tex, dT0, T0] = findTempsForFit_coeff(rhoStruct, TexCoeffs, dTCoeffs, T0Coeffs);

OlbDens = evalMajorSpecies(rhoStruct, coeffStruct.OCoeff, numBiasesStruct.O);
N2lbDens = evalMajorSpecies(rhoStruct, coeffStruct.N2Coeff, numBiasesStruct.N2);
HelbDens = evalMajorSpecies(rhoStruct, coeffStruct.HeCoeff, numBiasesStruct.He);
ArlbDens = evalMajorSpecies(rhoStruct, coeffStruct.ArCoeff, numBiasesStruct.Ar);
O2lbDens = exp(coeffStruct.O2Coeff);

ilRho = computeRho(T0, dT0, Tex, rhoStruct.Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);

if nargout > 1
    [~, msisRho] = computeMsis(rhoStruct);
    [~, dtmRho] = computeDtm(rhoStruct);
end

end

function [Tex, dT0, T0] = findTempsForFit_coeff(varStruct, TexCoeffs, dTCoeffs, T0Coeffs)

Tex_est = evalTex(varStruct, TexCoeffs);

T0 = clamp(200, evalT0(varStruct, T0Coeffs), 1000);
dT0 = clamp(1, evalDT(varStruct, dTCoeffs), 30);
Tex = clamp(T0+1, Tex_est, 5000);

end