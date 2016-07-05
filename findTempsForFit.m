function [Tex, dT0, T0] = findTempsForFit(varStruct, TexStruct, dTCoeffs, T0Coeffs, coeff)

origNumBiases = varStruct.numBiases; varStruct.numBiases = 0;
Tex_est = evalTex(varStruct, coeff(TexStruct.coeffInd));

T0 = clamp(200, evalT0(varStruct, T0Coeffs), 1000);
dT0 = clamp(1, evalDT(varStruct, dTCoeffs), 30);
Tex = clamp(T0+1, Tex_est, 5000);

varStruct.numBiases = origNumBiases;

end