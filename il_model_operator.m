function [rho, outputStruct] =...
    il_model_operator(state, S)

[Tex, dT0, T0] = findTempsForFit(S, S.TexCoeff, S.dTCoeff, S.T0Coeff);

% Add scalar corrections to dT0 and T0
T0 = clamp(300, T0 + state(1), 700);
dT0 = clamp(1, dT0 + state(2), 20);

% Add T2 correction to Tex
Tex = clamp(T0+1, Tex + bulge(state(3:11), S), 5000);

OlbDens = evalMajorSpecies(S, S.OCoeff, S.O_numBiases);
N2lbDens = evalMajorSpecies(S, S.N2Coeff, S.N2_numBiases);
HelbDens = evalMajorSpecies(S, S.HeCoeff, S.He_numBiases);
ArlbDens = evalMajorSpecies(S, S.ArCoeff, S.Ar_numBiases);
O2lbDens = exp(S.O2Coeff);

Z = computeGeopotentialHeight(S.altitude);

[rho, O, N2, He, Ar, O2, T] = ...
    computeRho(T0, dT0, Tex, Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);

if nargout > 1
    outputStruct = struct('O', O, 'N2', N2, 'He', He, 'Ar', Ar, 'O2', O2, 'Tex', Tex,...
                            'dT', dT0, 'T', T, 'T0', T0);
end

end