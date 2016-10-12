function newTleMap = selectTLEs(tleMap, oldestOrNewestOrCustom, getIndexVector)

objects = keys(tleMap);
if strcmpi(oldestOrNewestOrCustom, 'oldest')
    oldest = true;
    custom = false;
elseif strcmpi('index', oldestOrNewestOrCustom)
    oldest = false;
    custom = true;
else
    oldest = false;
    custom = false;
end

newTleMap = containers.Map('KeyType', 'double', 'ValueType', 'any');

for i = 1:length(objects)
    sgp4Info = tleMap(objects{i}).sgp4info;
    if oldest
        sgp4InfoReduced = sgp4Info(1);
    elseif custom
        sgp4InfoReduced = sgp4Info(getIndexVector(i));
    else
        sgp4InfoReduced = sgp4Info(end);
    end
    newTleMap(objects{i}) = struct('sgp4info', sgp4InfoReduced, 'Btrue', tleMap(objects{i}).Btrue);
end

end