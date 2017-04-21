function [fixStruct] = removeDataPoints(fixStruct, removeInd, removeBiasMatRows, fixSatIndices, fixZ, fixWeights)

if (islogical(removeInd) && length(fixStruct.timestamps) ~= length(removeInd))
    error('Lengths of struct data and remove indices (logical vector) must be the same!')
end

if isfield(fixStruct, 'data'); fixStruct.data(removeInd) = []; end
fixStruct.timestamps(removeInd) = [];
fixStruct.latitude(removeInd) = [];
fixStruct.longitude(removeInd) = [];
fixStruct.solarTime(removeInd) = [];
fixStruct.altitude(removeInd) = [];
if isfield(fixStruct, 'aeInt');fixStruct.aeInt(removeInd,:) = [];end
if isfield(fixStruct, 'F');fixStruct.F(removeInd) = [];end
if isfield(fixStruct, 'FA');fixStruct.FA(removeInd) = [];end
if isfield(fixStruct, 'apNow'); fixStruct.apNow(removeInd) = [];end
if isfield(fixStruct, 'ap3h');fixStruct.ap3h(removeInd) = [];end
if isfield(fixStruct, 'ap6h');fixStruct.ap6h(removeInd) = [];end
if isfield(fixStruct, 'ap9h');fixStruct.ap9h(removeInd) = [];end
if isfield(fixStruct, 'ap12To33h');fixStruct.ap12To33h(removeInd) = [];end
if isfield(fixStruct, 'ap36To57h');fixStruct.ap36To57h(removeInd) = [];end
if isfield(fixStruct, 'Ap');fixStruct.Ap(removeInd) = [];end
if isfield(fixStruct, 'dst')
    fixStruct.dst(removeInd) = [];
end
if isfield(fixStruct, 'sigma')
    fixStruct.sigma(removeInd) = [];
end

if nargin > 2
    if removeBiasMatRows
        if isfield(fixStruct, 'biases') fixStruct.biases(removeInd,:) = []; end
    end
    if fixSatIndices
        satInd = zeros(1, length(removeInd));
        
        if isfield(fixStruct,'goce') satInd(fixStruct.goce) = 1; end
        if isfield(fixStruct,'champ') satInd(fixStruct.champ) = 2; end
        if isfield(fixStruct,'grace') satInd(fixStruct.grace) = 3; end
        if isfield(fixStruct,'swarm') satInd(fixStruct.swarm) = 4; end
        if isfield(fixStruct,'de2') satInd(fixStruct.de2) = 5; end
        if isfield(fixStruct,'aeC') satInd(fixStruct.aeC) = 6; end
        if isfield(fixStruct,'aeE') satInd(fixStruct.aeE) = 7; end
        if isfield(fixStruct,'aeENace') satInd(fixStruct.aeENace) = 8; end
        if isfield(fixStruct,'aeEOss') satInd(fixStruct.aeEOss) = 9; end
        if isfield(fixStruct,'guvi') satInd(fixStruct.guvi) = 10; end
        if isfield(fixStruct,'aeros') satInd(fixStruct.aeros) = 11; end
        if isfield(fixStruct,'aeCOss') satInd(fixStruct.aeCOss) = 12; end
        satInd(removeInd) = [];
        if isfield(fixStruct,'goce') fixStruct.goce = find(satInd == 1); end
        if isfield(fixStruct,'champ') fixStruct.champ = find(satInd == 2); end
        if isfield(fixStruct,'grace') fixStruct.grace = find(satInd == 3); end
        if isfield(fixStruct,'swarm') fixStruct.swarm = find(satInd == 4); end
        if isfield(fixStruct,'de2') fixStruct.de2 = find(satInd == 5); end
        if isfield(fixStruct,'aeC') fixStruct.aeC = find(satInd == 6); end
        if isfield(fixStruct,'aeE') fixStruct.aeE = find(satInd == 7); end
        if isfield(fixStruct,'aeENace') fixStruct.aeENace = find(satInd == 8); end
        if isfield(fixStruct,'aeEOss') fixStruct.aeEOss = find(satInd == 9); end
        if isfield(fixStruct,'guvi') fixStruct.guvi = find(satInd == 10); end
        if isfield(fixStruct,'aeros') fixStruct.aeros = find(satInd == 11); end
        if isfield(fixStruct,'aeCOss') fixStruct.aeCOss = find(satInd == 12); end
    end
    if fixZ
        if isfield(fixStruct, 'Z') fixStruct.Z(removeInd) = []; end
        if isfield(fixStruct, 'doy') fixStruct.doy(removeInd) = []; end
    end
    if fixWeights
        if isfield(fixStruct, 'weights') fixStruct.weights(removeInd) = []; end
    end
end

fixStruct.dataEnd = length(fixStruct.timestamps);

end