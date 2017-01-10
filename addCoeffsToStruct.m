function S = addCoeffsToStruct(S, OStruct, HeStruct, N2Struct, ArStruct, O2Struct)

load optCoeff

S.TexCoeff = optCoeff(TexInd);
S.dTCoeff = dTCoeffs;
S.T0Coeff = T0Coeffs;
S.OCoeff = optCoeff(OInd);
S.HeCoeff = optCoeff(HeInd);
S.N2Coeff = optCoeff(N2Ind);
S.ArCoeff = optCoeff(ArInd);
S.O2Coeff = optCoeff(O2Ind);

S.O_numBiases = OStruct.numBiases;
S.He_numBiases = HeStruct.numBiases;
S.N2_numBiases = N2Struct.numBiases;
S.Ar_numBiases = ArStruct.numBiases;
S.O2_numBiases = O2Struct.numBiases;

end