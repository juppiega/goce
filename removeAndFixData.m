function [rhoStruct, TempStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct, removeIndRho] = ...
    removeAndFixData(rhoStruct, aeThreshold, TempStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, lbDTStruct, lbT0Struct, quietData, aeQuietLimit)
if nargin <= 11
    aeQuietLimit = 800;
end

rhoStruct.numBiases = 0;
if aeThreshold <= 0
    removeInd = ~ismember(1:length(rhoStruct.data), 1:1:length(rhoStruct.data)) | rhoStruct.data' <= 0; % !!!!!!!!!! TESTAUS
else
    removeInd = rhoStruct.aeInt(:,4) < aeThreshold;
end
if exist('quietData','var')
    [~,~,stormInd] = findStormsForSat(rhoStruct,'ae',aeQuietLimit,0,2,true);
    if quietData
        removeInd(stormInd) = true;
    else
        removeInd(~stormInd) = true;
    end

%     quietInd = all(rhoStruct.aeInt(:,4:18) < aeQuietLimit, 2);
%     if quietData
%         removeInd(~quietInd) = true;
%     else
%         removeInd(quietInd) = true;
%     end
end
% removeChamp = removeInd(rhoStruct.champ);
% removeChamp(rhoStruct.data(rhoStruct.champ) > 1E-10) = true;
% removeInd(rhoStruct.champ) = removeChamp;
% removeGrace = removeInd(rhoStruct.grace);
% removeGrace(rhoStruct.data(rhoStruct.grace) > 2E-11) = true;
% removeInd(rhoStruct.grace) = removeGrace;

relErr = rhoStruct.sigma ./ rhoStruct.data;
removeInd(relErr > 0.15) = true;
removeIndRho = removeInd;
rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, true, true);
% % GOCE has bad lst, F10.7 coverage.
% goceWeight = 0.5*0.5*(1 - length(rhoStruct.goce) / length(rhoStruct.data));
% champWeight = 1 - length(rhoStruct.champ)  / length(rhoStruct.data);
% graceWeight = 1 - length(rhoStruct.grace) / length(rhoStruct.data);
goceWeight = 0.5*0.5 * length(rhoStruct.goce) / length(rhoStruct.grace);
champWeight = 1;
graceWeight = 1;
rhoStruct.weights = zeros(length(rhoStruct.data),1);
rhoStruct.weights(rhoStruct.goce) = goceWeight;
rhoStruct.weights(rhoStruct.champ) = champWeight;
rhoStruct.weights(rhoStruct.grace) = graceWeight;
rhoStruct.numBiases = 0;

[yr,~,~,~,~,~] = datevec(rhoStruct.timestamps);
yearVec = [yr, repmat([1,1,0,0,0], length(yr), 1)];
rhoStruct.doy = rhoStruct.timestamps - datenum(yearVec) + 1;
rhoStruct = computeGeopotentialHeight(rhoStruct);

%if nargin == 2 && nargout == 1
if nargout == 1 || nargout == 10;
    return
end

% Remove bad temperature observations
removeInd = TempStruct.data <= 300 | TempStruct.data > 10000;
% For an unknown reason, there is one erroneous value in the local solar
% time present in the original ASCII data files.
removeInd(TempStruct.solarTime > 24 | TempStruct.solarTime < 0) = true(1);
% Remove observations with altitude > 900 km or < 250 or temperature > 3000 K.
removeInd(TempStruct.data > 3000 | TempStruct.altitude > 900 | TempStruct.altitude < 250) = true(1);
if nargin >= 11
    [~,~,stormInd] = findStormsForSat(TempStruct,'ae',aeQuietLimit,0,2,true);
    if quietData
        removeInd(stormInd) = true;
    else
        removeInd(~stormInd) = true;
    end
%     quietInd = all(TempStruct.aeInt < aeQuietLimit, 2);
%     if quietData
%         removeInd(~quietInd) = true;
%     else
%         removeInd(quietInd) = true;
%     end
end
TempStruct = removeDataPoints(TempStruct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(TempStruct.de2) = 1;
satInd(TempStruct.aeC) = 2;
satInd(TempStruct.aeE) = 3;
satInd(removeInd) = [];
TempStruct.de2 = find(satInd == 1);
TempStruct.aeC = find(satInd == 2);
TempStruct.aeE = find(satInd == 3);
TempStruct.numBiases = 0;
TempStruct.dataEnd = length(TempStruct.data);

% Remove DT observations
lbDTStruct.F = lbDTStruct.F10;
lbDTStruct.FA = lbDTStruct.F81A;
lbDTStruct.ap12To33h = lbDTStruct.apAver12To33h;
lbDTStruct.ap36To57h = lbDTStruct.apAver36To57h;
lbDTStruct.Ap = lbDTStruct.ApDaily;
lbDTStruct.altitude = lbDTStruct.altitude(1) * ones(size(lbDTStruct.data));

removeInd = lbDTStruct.solarTime < 0;
removeInd(lbDTStruct.data <= 1 | lbDTStruct.data > 20) = true;
%removeInd(85 < lbDTStruct.zenithAngle & lbDTStruct.zenithAngle < 95) = true;
removeInd(lbDTStruct.apNow > 15) = true;
%removeInd(lbDTStruct.nightObservation == 2) = true; % Remove twilight observations.
%removeInd(lbDTStruct.nightObservation == 1 & abs(lbDTStruct.latitude) > 52) = true;
if nargin >= 11
    [~,~,stormInd] = findStormsForSat(lbDTStruct,'ae',aeQuietLimit,0,2,true);
    if quietData
        removeInd(stormInd) = true;
    else
        removeInd(~stormInd) = true;
    end
end
lbDTStruct = removeDataPoints(lbDTStruct, removeInd);
lbDTStruct.auroraFlag(removeInd) = [];
lbDTStruct.zenithAngle(removeInd) = [];
lbDTStruct.nightObservation(removeInd) = [];
lbDTStruct.numBiases = 0;

% Remove T0 observations
lbT0Struct.F = lbT0Struct.F10;
lbT0Struct.FA = lbT0Struct.F81A;
lbT0Struct.ap12To33h = lbT0Struct.apAver12To33h;
lbT0Struct.ap36To57h = lbT0Struct.apAver36To57h;
lbT0Struct.Ap = lbT0Struct.ApDaily;
lbT0Struct.fundPulseLen = lbT0Struct.fundPulsLen;

removeInd = lbT0Struct.apNow > 15;
spuriousLocationChanges = find(abs(diff(lbT0Struct.latitude)) > 1.0) + 1;
removeInd(spuriousLocationChanges) = true;
removeInd(lbT0Struct.data <= 250 | lbT0Struct.data >= 850) = true;
removeInd(lbT0Struct.Ti_err > 1000) = true(1);
removeInd(lbT0Struct.fundPulseLen > 100*1E-6) = true; % Remove > 100 us (>~15 km) pulses.
if nargin >= 11
    [~,~,stormInd] = findStormsForSat(lbT0Struct,'ae',aeQuietLimit,0,2,true);
    if quietData
        removeInd(stormInd) = true;
    else
        removeInd(~stormInd) = true;
    end
end
lbT0Struct = removeDataPoints(lbT0Struct, removeInd);
lbT0Struct.Ti(removeInd) = [];
lbT0Struct.Ti_err(removeInd) = [];
lbT0Struct.Tn_err(removeInd) = [];
lbT0Struct.pulseLen(removeInd) = [];
lbT0Struct.fundPulseLen(removeInd) = [];
lbT0Struct.vertRes(removeInd) = [];
lbT0Struct.index(removeInd) = [];
lbT0Struct.numBiases = 0;

% Remove bad oxygen observations
removeInd = OStruct.data <= 1E6 | OStruct.data > 1E11;
if nargin >= 11
    [~,~,stormInd] = findStormsForSat(OStruct,'ae',aeQuietLimit,0,2,true);
    if quietData
        removeInd(stormInd) = true;
    else
        removeInd(~stormInd) = true;
    end
end
OStruct = removeDataPoints(OStruct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(OStruct.de2) = 1;
satInd(OStruct.aeC) = 2;
satInd(OStruct.aeENace) = 3;
satInd(OStruct.aeEOss) = 4;
satInd(OStruct.guvi) = 5;
satInd(removeInd) = [];
OStruct.de2 = find(satInd == 1);
OStruct.aeC = find(satInd == 2);
OStruct.aeENace = find(satInd == 3);
OStruct.aeEOss = find(satInd == 4);
OStruct.guvi = find(satInd == 5);
OStruct.numBiases = 5; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OStruct.biases = zeros(length(OStruct.data), OStruct.numBiases);
OStruct.biases(OStruct.de2, 1) = 1;
OStruct.biases(OStruct.aeC, 2) = 1;
OStruct.biases(OStruct.aeENace, 3) = 1;
OStruct.biases(OStruct.aeEOss, 4) = 1;
OStruct.biases(OStruct.guvi, 5) = 1;
OStruct.dataEnd = length(OStruct.data);
OStruct.name = 'O';

% Remove N2 observations.
removeInd = N2Struct.data <= 1E2 | N2Struct.data > 0.5E12 | N2Struct.altitude > 600 | N2Struct.altitude < 140;
if nargin >= 11
    [~,~,stormInd] = findStormsForSat(N2Struct,'ae',aeQuietLimit,0,2,true);
    if quietData
        removeInd(stormInd) = true;
    else
        removeInd(~stormInd) = true;
    end
end
N2Struct = removeDataPoints(N2Struct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(N2Struct.de2) = 1;
satInd(N2Struct.aeC) = 2;
satInd(N2Struct.aeENace) = 3;
satInd(N2Struct.aeEOss) = 4;
satInd(N2Struct.aeros) = 5;
satInd(N2Struct.guvi) = 6;
satInd(removeInd) = [];
N2Struct.de2 = find(satInd == 1);
N2Struct.aeC = find(satInd == 2);
N2Struct.aeENace = find(satInd == 3);
N2Struct.aeEOss = find(satInd == 4);
N2Struct.aeros = find(satInd == 5);
N2Struct.guvi = find(satInd == 6);
N2Struct.numBiases = 6;
N2Struct.biases = zeros(length(N2Struct.data), N2Struct.numBiases);
N2Struct.biases(N2Struct.de2, 1) = 1;
N2Struct.biases(N2Struct.aeC, 2) = 1;
N2Struct.biases(N2Struct.aeENace, 3) = 1;
N2Struct.biases(N2Struct.aeEOss, 4) = 1;
N2Struct.biases(N2Struct.aeros, 5) = 1;
N2Struct.biases(N2Struct.guvi, 6) = 1;
N2Struct.dataEnd = length(N2Struct.data);
N2Struct.name = 'N2';

% Remove He observations
removeInd = HeStruct.data <= 1E5 | HeStruct.data > 1E10;
if nargin >= 11
    [~,~,stormInd] = findStormsForSat(HeStruct,'ae',aeQuietLimit,0,2,true);
    if quietData
        removeInd(stormInd) = true;
    else
        removeInd(~stormInd) = true;
    end
end
HeStruct = removeDataPoints(HeStruct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(HeStruct.de2) = 1;
satInd(HeStruct.aeC) = 2;
satInd(HeStruct.aeENace) = 3;
satInd(HeStruct.aeEOss) = 4;
satInd(HeStruct.aeros) = 5;
satInd(removeInd) = [];
HeStruct.de2 = find(satInd == 1);
HeStruct.aeC = find(satInd == 2);
HeStruct.aeENace = find(satInd == 3);
HeStruct.aeEOss = find(satInd == 4);
HeStruct.aeros = find(satInd == 5);
HeStruct.numBiases = 5;
HeStruct.biases = zeros(length(HeStruct.data), HeStruct.numBiases);
HeStruct.biases(HeStruct.de2, 1) = 1;
HeStruct.biases(HeStruct.aeC, 2) = 1;
HeStruct.biases(HeStruct.aeENace, 3) = 1;
HeStruct.biases(HeStruct.aeEOss, 4) = 1;
HeStruct.biases(HeStruct.aeros, 5) = 1;
HeStruct.dataEnd = length(HeStruct.data);
HeStruct.name = 'He';

% Remove Ar observations
removeInd = ArStruct.data <= 100 | ArStruct.data > 1E16; % KORJAA!!!!
if nargin >= 11
    [~,~,stormInd] = findStormsForSat(ArStruct,'ae',aeQuietLimit,0,2,true);
    if quietData
        removeInd(stormInd) = true;
    else
        removeInd(~stormInd) = true;
    end
end
ArStruct = removeDataPoints(ArStruct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(ArStruct.de2) = 1;
satInd(ArStruct.aeros) = 2;
satInd(removeInd) = [];
ArStruct.de2 = find(satInd == 1);
ArStruct.aeros = find(satInd == 2);
ArStruct.numBiases = 2;
ArStruct.biases = zeros(length(ArStruct.data), ArStruct.numBiases);
ArStruct.biases(ArStruct.de2, 1) = 1;
ArStruct.biases(ArStruct.aeros, 2) = 1;
ArStruct.dataEnd = length(ArStruct.data);
ArStruct.name = 'Ar';

% Remove O2 observations
removeInd = O2Struct.data <= 0 | O2Struct.data > 1E17; % KORJAA!!!!
if nargin >= 11
    quietInd = all(O2Struct.aeInt < aeQuietLimit, 2);
    if quietData
        removeInd(~quietInd) = true;
    else
        removeInd(quietInd) = true;
        if sum(~removeInd) == 0
            removeInd(1) = false;
        end
    end
end
removeInd(O2Struct.aeEOss) = true;
removeInd(O2Struct.apNow > 15) = true;
removeInd(O2Struct.altitude > 220) = true;
O2Struct = removeDataPoints(O2Struct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(O2Struct.aeCOss) = 1;
satInd(removeInd) = [];
O2Struct.aeCOss = find(satInd == 1);
O2Struct.numBiases = 0;
%O2Struct.biases = zeros(length(O2Struct.data), O2Struct.numBiases);
%O2Struct.biases(O2Struct.aeCOss, 1) = 1;
O2Struct.dataEnd = length(O2Struct.data);
O2Struct.name = 'O2';

TempStruct = computeGeopotentialHeight(TempStruct);
OStruct = computeGeopotentialHeight(OStruct);
N2Struct = computeGeopotentialHeight(N2Struct);
HeStruct = computeGeopotentialHeight(HeStruct);
ArStruct = computeGeopotentialHeight(ArStruct);
O2Struct = computeGeopotentialHeight(O2Struct);


end
