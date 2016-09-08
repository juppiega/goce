function [simulatedObservations] = dummyThermosphere(state, observationStruct)

T0 = state(1);
dT0 = state(2);
Tex = state(3);
OlbDens = state(4);
N2lbDens = state(5);
HelbDens = state(6);
ArlbDens = 8.6E8;
O2lbDens = ;

computeRho(T0, dT0, Tex, Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);

end