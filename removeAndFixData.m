function [TempStruct, OStruct, N2Struct, HeStruct, rhoStruct] = removeAndFixData(TempStruct, OStruct, N2Struct, HeStruct, rhoStruct)

% Make GOCE observations unbiased.
rhoStruct.data(rhoStruct.goce) = rhoStruct.data(rhoStruct.goce) * 1.23;
rhoStruct.numBiases = 0;
removeInd = ~ismember(1:length(rhoStruct.data), 1:1:length(rhoStruct.data)) | rhoStruct.data' <= 0; % !!!!!!!!!! TESTAUS
rhoStruct = removeDataPoints(rhoStruct, removeInd, false, true, false, false);
goceWeight = 0.5*(1 - length(rhoStruct.goce) / length(rhoStruct.data));
champWeight = 1 - length(rhoStruct.champ)  / length(rhoStruct.data);
graceWeight = 1 - length(rhoStruct.grace) / length(rhoStruct.data);
rhoStruct.weights = zeros(length(rhoStruct.data),1);
rhoStruct.weights(rhoStruct.goce) = goceWeight;
rhoStruct.weights(rhoStruct.champ) = champWeight;
rhoStruct.weights(rhoStruct.grace) = graceWeight;
rhoStruct.numBiases = 0;

TempStruct.T0 = 507;
TempStruct.dT0 = 12.6;
% Remove bad temperature observations
removeInd = TempStruct.data <= TempStruct.T0 | TempStruct.data > 10000;
% For an unknown reason, there is one erroneous value in the local solar
% time present in the original ASCII data files.
removeInd(TempStruct.solarTime > 24 | TempStruct.solarTime < 0) = true(1);
% Remove observations with altitude > 900 km or < 250 or temperature > 3000 K.
removeInd(TempStruct.data > 3000 | TempStruct.altitude > 900 | TempStruct.altitude < 250) = true(1);
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

% Remove bad oxygen observations
removeInd = OStruct.data <= 1E6 | OStruct.data > 0.8E10;
OStruct = removeDataPoints(OStruct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(OStruct.de2) = 1;
satInd(OStruct.aeC) = 2;
satInd(OStruct.aeENace) = 3;
satInd(OStruct.aeEOss) = 4;
satInd(removeInd) = [];
OStruct.de2 = find(satInd == 1);
OStruct.aeC = find(satInd == 2);
OStruct.aeENace = find(satInd == 3);
OStruct.aeEOss = find(satInd == 4);
OStruct.numBiases = 0; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OStruct.biases = zeros(length(OStruct.data), OStruct.numBiases);
OStruct.biases(OStruct.de2, 1) = 1;
OStruct.biases(OStruct.aeC, 2) = 1;
OStruct.biases(OStruct.aeENace, 3) = 1;
OStruct.biases(OStruct.aeEOss, 4) = 1;
OStruct.dataEnd = length(OStruct.data);
OStruct.name = 'O';

% Remove N2 observations.
removeInd = N2Struct.data <= 1E2 | N2Struct.data > 0.5E12 | N2Struct.altitude > 600 | N2Struct.altitude < 140;
N2Struct = removeDataPoints(N2Struct, removeInd);
satInd = zeros(1, length(removeInd));
satInd(N2Struct.de2) = 1;
satInd(N2Struct.aeC) = 2;
satInd(N2Struct.aeENace) = 3;
satInd(N2Struct.aeEOss) = 4;
satInd(N2Struct.aeros) = 5;
satInd(removeInd) = [];
N2Struct.de2 = find(satInd == 1);
N2Struct.aeC = find(satInd == 2);
N2Struct.aeENace = find(satInd == 3);
N2Struct.aeEOss = find(satInd == 4);
N2Struct.aeros = find(satInd == 5);
N2Struct.numBiases = 0;
N2Struct.dataEnd = length(N2Struct.data);
N2Struct.name = 'N2';

% Remove He observations
removeInd = HeStruct.data <= 1E5 | HeStruct.data > 1E10;
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
HeStruct.numBiases = 0;
HeStruct.dataEnd = length(HeStruct.data);
HeStruct.name = 'He';

TempStruct = computeGeopotentialHeight(TempStruct);
OStruct = computeGeopotentialHeight(OStruct);
N2Struct = computeGeopotentialHeight(N2Struct);
HeStruct = computeGeopotentialHeight(HeStruct);
rhoStruct = computeGeopotentialHeight(rhoStruct);

end

function [addStruct] = computeGeopotentialHeight(addStruct)
% TODO: Lisaa leveyspiirin ja pyorimisen vaikutus.

R = 6371E3;
addStruct.Z = (R * addStruct.altitude) ./ (R + addStruct.altitude*1000);

end