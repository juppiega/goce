function [rhoStruct, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct] = ...
    subsampleStructs(rhoStruct, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, subsampPercent)
% subsampPercent ]0,100]

subsamp = round(1 / (subsampPercent / 100));

N = length(rhoStruct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, false, false);

N = length(TexStruct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
TexStruct = removeDataPoints(TexStruct, removeInd);

N = length(OStruct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
OStruct = removeDataPoints(OStruct, removeInd);

N = length(N2Struct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
N2Struct = removeDataPoints(N2Struct, removeInd);

N = length(HeStruct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
HeStruct = removeDataPoints(HeStruct, removeInd);

N = length(ArStruct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
ArStruct = removeDataPoints(ArStruct, removeInd);

N = length(O2Struct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
O2Struct = removeDataPoints(O2Struct, removeInd);

end