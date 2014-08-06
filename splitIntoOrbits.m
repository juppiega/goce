function [orbits, exactOrbitIndices] = splitIntoOrbits(latitude, timestamps10s)
%
oneHinSec = 60 * 60;
if satelliteIsGoingSouth(latitude)
    orbitEndIndices = [find(latitude(1:end-1) < latitude(2:end) | timestamps10s(1:end -1) - timestamps10s(2:end) < -1 * oneHinSec); length(latitude)];
else
    orbitEndIndices = [find(latitude(1:end-1) > latitude(2:end) | timestamps10s(1:end -1) - timestamps10s(2:end) < -1 * oneHinSec); length(latitude)];
end
orbitBeginIndices = [1; (orbitEndIndices(1:end-1) + 1)];

orbits = [orbitBeginIndices orbitEndIndices];

orbitsToDelete = orbits(:,2) == orbits(:,1);
orbits = orbits(~orbitsToDelete, :);

exactOrbitIndices = [];
for i = 1:length(orbits(:,1))
    exactOrbitIndices = [exactOrbitIndices (orbits(i,1) : orbits(i,2))];
end

end