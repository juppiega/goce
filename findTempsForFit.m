function [Tex, dT0, T0] = findTempsForFit(varStruct, TexStruct, coeff)

origNumBiases = varStruct.numBiases; varStruct.numBiases = 0;
Tex_est = evalTex(varStruct, coeff(TexStruct.coeffInd));
varStruct.numBiases = origNumBiases;

T0 = TexStruct.T0;
dT0 = TexStruct.dT0;
Tex = clamp(T0+1, Tex_est, 5000);

end