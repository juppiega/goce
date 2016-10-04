function [simulatedObservations, outputStruct] =...
    dummyThermosphere(state, observationStruct)

T0 = 507;
dT0 = 12;
Tex = max(state(1), T0+1);
OlbDens = max(state(2), log(1E8));
N2lbDens = max(state(3), log(1E9));
HelbDens = log(2.5E7);
ArlbDens = log(8.6E8);
O2lbDens = log(1.3E5);
Z = computeGeopotentialHeight(observationStruct.altitude);

[simulatedObservations, O, N2, He, Ar, O2, T] = computeRho(T0, dT0, Tex, Z, exp(OlbDens), exp(N2lbDens), ...
        exp(HelbDens), exp(ArlbDens), exp(O2lbDens));
    
if nargout > 1
    outputStruct = struct('O', O, 'N2', N2, 'He', He, 'Ar', Ar, 'O2', O2, 'Tex', Tex,...
                            'dT', dT0, 'T', T, 'T0', T0);
end

end