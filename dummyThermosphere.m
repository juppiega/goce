function [rho, outputStruct] =...
    dummyThermosphere(state, observationStruct, index)

T0 = clamp(300, 507 + state(1), 700);
dT0 = clamp(1, 12 + state(2), 20);
Tex = clamp(T0+1, 1030 + bulge(state(3:11), observationStruct), 5000);
OlbDens = log(8.47E10);
N2lbDens = log(3.2E11);
HelbDens = log(2.5E7);
ArlbDens = log(8.6E8);
O2lbDens = log(1.3E5);
Z = computeGeopotentialHeight(observationStruct.altitude);

[rho, O, N2, He, Ar, O2, T] = computeRho(T0, dT0, Tex, Z, exp(OlbDens), exp(N2lbDens), ...
        exp(HelbDens), exp(ArlbDens), exp(O2lbDens));
rho = log(rho);
    
if nargout > 1
    outputStruct = struct('O', O, 'N2', N2, 'He', He, 'Ar', Ar, 'O2', O2, 'Tex', Tex,...
                            'dT', dT0, 'T', T, 'T0', T0);
end

end