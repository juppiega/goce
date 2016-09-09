function [simulatedObservations, outputStruct] =...
    dummyThermosphere(state, observationStruct)

T0 = clamp(300,state(1),700);
dT0 = clamp(1,state(2),30);
Tex = max(state(3), T0+1);
OlbDens = max(state(4), log(1E8));
N2lbDens = max(state(5), log(1E9));
HelbDens = max(state(6), log(1E6));
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