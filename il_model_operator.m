function [rho, outputStruct] =...
    il_model_operator(state, S, index)

if state(1) > 0
    S.F = state(1);
    S.FA = state(1);
end

T0 = clamp(200, evalT0(S, S.T0Coeff), 1000);
dT0 = clamp(1, evalDT(S, S.dTCoeff), 30);
Tex = clamp(T0+1, evalTex(S, S.TexCoeff), 5000);

% Add scalar corrections to dT0 and T0
T0 = clamp(300, T0 + state(2), 700);
dT0 = clamp(1, dT0 + state(3), 20);

% Add T2 correction to Tex
Tex = clamp(T0+1, Tex + bulge(state(4:11), S), 5000);

OlbDens = evalMajorSpecies(S, S.OCoeff, S.O_numBiases);
N2lbDens = evalMajorSpecies(S, S.N2Coeff, S.N2_numBiases);
HelbDens = evalMajorSpecies(S, S.HeCoeff, S.He_numBiases);
ArlbDens = evalMajorSpecies(S, S.ArCoeff, S.Ar_numBiases);
O2lbDens = exp(S.O2Coeff);

Z = computeGeopotentialHeight(S.altitude);

[rho, O, N2, He, Ar, O2, T] = ...
    computeRho(T0, dT0, Tex, Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);
rho = log(rho);

if nargout > 1
    outputStruct = struct('O', O, 'N2', N2, 'He', He, 'Ar', Ar, 'O2', O2, 'Tex', Tex,...
                            'dT', dT0, 'T', T, 'T0', T0);
end

end