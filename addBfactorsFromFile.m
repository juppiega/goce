function tleMap = addBfactorsFromFile(tleMap)

load Bfactors.dat

satellites = Bfactors(:,1);
B = Bfactors(:,2);

obj = keys(tleMap);

for i = 1:length(obj)
    ind = satellites == obj{i};
    if sum(ind) == 0
        error(['Could not find Btrue for object: ', obj{i}]);
    end
    if sum(ind) > 1
        error(['File contained Btrue twice (or more) for object: ', obj{i}, '. A satellite must have only one Btrue value!']);
    end
    tleMap(obj(i)).Btrue = B(ind);
end

end