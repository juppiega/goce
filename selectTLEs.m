function tleMap = selectTLEs(tleMap, oldestOrNewest)

objects = keys(tleMap);
if strcmpi(oldestOrNewest, 'oldest')
    oldest = true;
else
    oldest = false;
end

for i = 1:lenth(objects)
    sgp4Info = tleMap(objects(i));
    if oldest
        sgp4InfoReduced = sgp4Info(1);
    else
        sgp4InfoReduced = sgp4Info(end);
    end
    tleMap(objects(i)) = sgp4InfoReduced;
end

end