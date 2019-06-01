function testIL()

S.timestamps = datenum('2019-04-15 16:00:00');
S.altitude = 400;
S.latitude = 10;
S.longitude = -60;
S.solarTime = 12;
S.F = 100;
S.FA = 100;
S.aeInt = zeros(1, 24) + 600;

S = computeVariablesForFit(S);
load optCoeff.mat

[Tex, dT0, T0] = findTempsForFit_this(S, optCoeff(TexInd), dTCoeffs, T0Coeffs);
OlbDens = evalMajorSpecies(S, optCoeff(OInd), 5);
N2lbDens = evalMajorSpecies(S, optCoeff(N2Ind), 6);
HelbDens = evalMajorSpecies(S, optCoeff(HeInd), 5);
O2lbDens = exp(optCoeff(O2Ind));

[rho, ~, ~, ~, ~, ~, T] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens)

end

function [Tex, dT0, T0] = findTempsForFit_this(varStruct, TexCoeffs, dTCoeffs, T0Coeffs, coeff)

Tex_est = evalTex(varStruct, TexCoeffs);

T0 = clamp(200, evalT0(varStruct, T0Coeffs), 1000);
dT0 = clamp(1, evalDT(varStruct, dTCoeffs), 30);
Tex = clamp(T0+1, Tex_est, 5000);

end