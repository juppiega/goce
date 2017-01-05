function [rhoStruct, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct] = ...
    subsampleStructs(rhoStruct, TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, subsampPercent)
% subsampPercent ]0,100]

subsamp = round(1 / (subsampPercent / 100));

N = length(rhoStruct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
rhoStruct = removeDataPoints(rhoStruct, removeInd, true, true, true, true);

N = length(TexStruct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
TexStruct = removeDataPoints(TexStruct, removeInd, true, true, true, true);

N = length(OStruct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
OStruct = removeDataPoints(OStruct, removeInd, true, true, true, true);

N = length(N2Struct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
N2Struct = removeDataPoints(N2Struct, removeInd, true, true, true, true);

N = length(HeStruct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
HeStruct = removeDataPoints(HeStruct, removeInd, true, true, true, true);

N = length(ArStruct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
ArStruct = removeDataPoints(ArStruct, removeInd, true, true, true, true);

N = length(O2Struct.data);
removeInd = ~ismember(1:N, 1:subsamp:N);
O2Struct = removeDataPoints(O2Struct, removeInd, true, true, true, true);

end