function goingSouth = satelliteIsGoingSouth(magneticLatitude)
%

orbitEndIndices = find(abs(diff(magneticLatitude)) > 45);
testIndex1 = round(mean([orbitEndIndices(3) orbitEndIndices(2)]));
testIndex2 = testIndex1 + 1;

if magneticLatitude(testIndex1) > magneticLatitude(testIndex2)
   goingSouth = 1; 
else
   goingSouth = 0;
end

end