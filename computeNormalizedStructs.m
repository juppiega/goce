function [normalizedData] = computeNormalizedStructs()

load('ilData.mat','OStruct');

OStruct = computeVariablesForFit(OStruct);
[T0, dT] = computeMsisDtmLb(OStruct);
Tex = computeMsis(OStruct);

normalizedData = computeNormalizedDensity(OStruct.data, OStruct.Z, 'O', Tex, dT, T0);

save('ONormalized.mat','normalizedData');

end