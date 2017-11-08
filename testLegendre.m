S.latitude = 0; S.longitude = 25; S.altitude = 400;
S.doy = 100;S.solarTime = 12;
S.timestamps = now;
S = computeVariablesForFit(S);
S.FA = 100; S.F = 100; S.aeInt = zeros(1,24);

load optCoeff.mat

Tex = evalTex(S, optCoeff(TexInd));
T0 = evalT0(S, T0Coeffs);
dT = evalDT(S, dTCoeffs);

O_lb = evalMajorSpecies(S, optCoeff(OInd), 5);
N2_lb = evalMajorSpecies(S, optCoeff(N2Ind), 6);
He_lb = evalMajorSpecies(S, optCoeff(HeInd), 5);
Ar_lb = evalMajorSpecies(S, optCoeff(ArInd), 2);
O2_lb = exp(optCoeff(O2Ind));

[rho, O, N2, He, Ar, O2, T] = ...
    computeRho(T0, dT, Tex, S.Z, O_lb, N2_lb, He_lb, Ar_lb, O2_lb);

fprintf('%16.12e\n',rho)
fprintf('%16.12e\n',T)
fprintf('%e\n',[O,N2,He,Ar,O2,Tex,T0,dT])