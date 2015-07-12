function [valAtSatLoc] = interpSatellite(lon, lat, altGrid, modelTime, modelField, satLon, satLat, satAlt, satTime, interpOption)

if satTime < modelTime(1) || satTime > modelTime(end)
    valAtSatLoc = nan(1);
    return
end

dt = modelTime(2) - modelTime(1);

i = find(lon <= satLon, 1, 'last');
j = find(lat <= satLat, 1, 'last');
if isempty(i)
    i = 1;
end
if isempty(j)
    j = 1;
end

if satLat < min(lat)
    jNext = 0;
else jNext = j + 1; 
end
if i == length(lon)
    iNext = 1;
elseif satLon < min(lon)
    iNext = length(lon);
else iNext = i + 1;
end

if jNext == 0
    interpLon = [lon; lon; lon; lon];
    interpLat = ones(size(interpLon)) * min(lat);
    lonInd = (1:length(lon))'; latInd = ones(size(lonInd))*1;
elseif jNext == length(lat) + 1
    interpLon = [lon; lon; lon; lon];
    interpLat = ones(size(interpLon)) * max(lat);
    lonInd = (1:length(lon))'; latInd = ones(size(lonInd))*length(lat);
else
    interpLon = repmat([lon(i); lon(i); lon(iNext); lon(iNext)], 4, 1);
    interpLat = repmat([lat(j); lat(jNext); lat(jNext); lat(j)], 4, 1);
    lonInd = [i;i;iNext;iNext];
    latInd = [j;jNext;jNext;j];
end


nearestTime = find(modelTime <= satTime, 1, 'last');
if nearestTime == length(modelTime)
    nextTime = nearestTime;
else nextTime = nearestTime + 1;
end

normalizedTime = (modelTime - modelTime(nearestTime)) / dt;

vals = zeros(2,1);
tiegcmCoord = [];
interpField = [];
for n = 1:2
    if n == 1
        t = nearestTime;
    else t = nextTime;
    end
    
    if satAlt < altGrid(i,j,1,t) || satAlt > altGrid(i,j,end-1,t)
        valAtSatLoc = nan(1);
        return
    end
    
    altCol = altGrid(i,j,:,t);
    k = find(altCol(1:end-1) <= satAlt, 1, 'last');
    
    tInd = ones(size(lonInd)) * t;
    nlev = size(altGrid, 3);
    interpAlt = [];
    %interpField = [];
    for m = -1:2
        if k + m < 1 || k + m > nlev-1
            k_Ind = ones(size(lonInd)) * k;
            km_Ind = ones(size(lonInd)) * (k-m);
            linInd1 = sub2ind(size(altGrid), lonInd, latInd, k_Ind, tInd);
            linInd2 = sub2ind(size(altGrid), lonInd, latInd, km_Ind, tInd);
            extrapAlt = 2 * altGrid(linInd1) - altGrid(linInd2);
            extrapField = modelField(linInd1).^2 ./ modelField(linInd2);
            interpAlt = [interpAlt; extrapAlt];
            interpField = [interpField; extrapField];
        else
            km_Ind = ones(size(lonInd)) * (k+m);
            linInd = sub2ind(size(altGrid), lonInd, latInd, km_Ind, tInd);
            interpAlt = [interpAlt; altGrid(linInd)];
            interpField = [interpField; modelField(linInd)];
        end
    end

    tiegcmCoordThisT = geod2ecef(interpLat, interpLon, interpAlt) / 1E6;
    if nearestTime ~= nextTime || n == 1
        tiegcmCoordThisT = [tiegcmCoordThisT, ones(size(tiegcmCoordThisT,1),1)*normalizedTime(t)];
    else
        ghostTime = normalizedTime(t) + 1;
        tiegcmCoordThisT = [tiegcmCoordThisT, ones(size(tiegcmCoordThisT,1),1)*ghostTime];
    end
    tiegcmCoord = [tiegcmCoord; tiegcmCoordThisT];
    %satCoord = geod2ecef(satLat, satLon, satAlt) / 1E6;
    %vals(n) = griddatan(double(tiegcmCoord), double(interpField), satCoord, 'linear');
end

satCoord = geod2ecef(satLat, satLon, satAlt) / 1E6;
normSatTime = (satTime - modelTime(nearestTime)) / dt;
satCoord = [satCoord, normSatTime];

valAtSatLoc = griddatan(double(tiegcmCoord), double(interpField), satCoord); %interp1(tInterp, vals, satTime);

end