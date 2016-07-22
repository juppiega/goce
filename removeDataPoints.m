function [fixStruct] = removeDataPoints(fixStruct, removeInd, removeBiasMatRows, fixSatIndices, fixZ, fixWeights)

if (islogical(removeInd) && length(fixStruct.data) ~= length(removeInd))
    error('Lengths of struct data and remove indices (logical vector) must be the same!')
end

fixStruct.data(removeInd) = [];
fixStruct.timestamps(removeInd) = [];
fixStruct.latitude(removeInd) = [];
fixStruct.longitude(removeInd) = [];
fixStruct.solarTime(removeInd) = [];
fixStruct.altitude(removeInd) = [];
fixStruct.aeInt(removeInd,:) = [];
fixStruct.F(removeInd) = [];
fixStruct.FA(removeInd) = [];
fixStruct.apNow(removeInd) = [];
fixStruct.ap3h(removeInd) = [];
fixStruct.ap6h(removeInd) = [];
fixStruct.ap9h(removeInd) = [];
fixStruct.ap12To33h(removeInd) = [];
fixStruct.ap36To57h(removeInd) = [];
fixStruct.Ap(removeInd) = [];

if nargin > 2
    if removeBiasMatRows
        fixStruct.biases(removeInd,:) = [];
    end
    if fixSatIndices
        satInd = zeros(1, length(removeInd));
        satInd(fixStruct.goce) = 1;
        satInd(fixStruct.champ) = 2;
        satInd(fixStruct.grace) = 3;
        satInd(removeInd) = [];
        fixStruct.goce = find(satInd == 1);
        fixStruct.champ = find(satInd == 2);
        fixStruct.grace = find(satInd == 3);
    end
    if fixZ
        fixStruct.Z(removeInd) = [];
    end
    if fixWeights
        fixStruct.weights(removeInd) = [];
    end
end

end