function readFiles()
%

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool();
end

% if(matlabpool('size')==0)
%     matlabpool;
% end


tic;
[ae, timestampsAeDatenum] = readAeFiles();

[F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, F10datenum, datenumToJulian] = readSolarIndexFiles();

[absB, timestampsAbsBDatenum, akasofuEpsilon, epsilonQualityFlag, timestampsEpsilonDatenum, vBz] ...
    = readSolarWindFiles(timestampsAeDatenum);

[density, longitude, latitude, altitude, solarTime, magneticLatitude, magneticLongitude, magneticLocalTime, crwindEast, crwindNorth, crwindUp, ...
 timestamps10sFixed, timestamps1min, timestamps1minFixed, timestampsDensityDatenum, doy, timestampsAbsB, timestampsEpsilon, firstDatenum, absDensityError, ...
 absWindError, noiseAffected, eclipseAffected, isMorningPass, ionThrusterActive] = readDensityFile(timestampsAeDatenum, timestampsAbsBDatenum, timestampsEpsilonDatenum);

goceData = struct('density', density * 1E-11, 'densityError', absDensityError * 1E-11, 'timestamps', timestampsDensityDatenum, 'latitude', latitude, 'longitude', longitude,...
                  'altitude', altitude/1000, 'solarTime', solarTime);


[champData, graceData, de2Data, aeCData, aeEData, aerosData] = computeOtherDensities();

load temperatureGradient.mat
load T0.mat

T0 = vertcatStructFields(T0);

[tiegcmDensityVariableAlt, tiegcmDensity270km] = readTiegcmFile(timestampsDensityDatenum);

[aeIntegrals, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0] = ...
    computeAeIntegrals(ae, timestampsAeDatenum, timestampsDensityDatenum, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0);

[apGoce, ap, timestampsAp, amGoce, dtc] = readApAndDtcFiles(timestampsDensityDatenum);

[F10, F81A, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0] = ...
    giveSolarInputForModels(F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, F10datenum, timestampsDensityDatenum, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0);

[apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, ApDaily, am3h, amAver24h, ap, timestamps3h, timestamps3hFixed, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0] =...
    giveApValsForMSIS(apGoce, amGoce, ap, timestampsAp, timestamps10sFixed, timestamps1minFixed, timestamps1min, timestampsDensityDatenum, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0);

if ~exist('goceVariables.mat', 'file')
[densityNoBg, msisDensityVariableAlt, msisDensity270km, msisDensity270kmNoAp, jb2008DensityVariableAlt, jb2008Density270km, jb2008Density270kmNoDtc, ...
    dtm2013Density270km, dtm2013DensityVariableAlt, dtm2013Density270kmNoAm, hwmU, hwmV, densityIndex, densityIndex1min, densityIndexNoBg, averagedDensity, averagedDensityNoBg, density3h]  = ...
    relateMsisToDensity(density, altitude, datenumToJulian, timestampsDensityDatenum, doy,...
    timestampsAeDatenum, timestamps1minFixed, timestamps3hFixed, solarTime, latitude, longitude, F10, F81A, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, am3h, amAver24h, dtc, ...
    ApDaily, apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h);

    measuredDensity = density;
    [morningTimestamps10s, eveningTimestamps10s, morningLatitude, eveningLatitude, morningDoy, eveningDoy] = ...
        splitBySolarTime(timestamps10sFixed, latitude, doy, solarTime);

    [morningGrid, morningBins] = computeTimeCells(morningTimestamps10s, firstDatenum, morningLatitude);
    [eveningGrid, eveningBins] = computeTimeCells(eveningTimestamps10s, firstDatenum, eveningLatitude);

    %msisSimulatedResidue = msisDensity270km - msisDensity270kmNoAp;

    [morningFourierGrid, eveningFourierGrid, morningLatDoy, eveningLatDoy] = computeFourierGrids(morningGrid, morningBins, eveningGrid, ...
        eveningBins, doy, morningLatitude, eveningLatitude, morningDoy, eveningDoy, timestamps10sFixed, densityIndex, aeIntegrals(:,9), median(F10));

    [hwmU, hwmV] = computeCrossTrackWind(hwmU, hwmV, crwindEast, crwindNorth);
end

if fclose('all') ~= 0
    display('File close unsuccesful - check if some of the files are reserved by another editor.')
end

if ~exist('ilData.mat', 'file') 
    fprintf('%s\n', 'Arranging satellite data by variable')
    [TempStruct, lbDTStruct, lbT0Struct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct] = arrangeByVar(goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0);
    fprintf('%s\n', 'Saving results to "ilData.mat" file')
    save('ilData.mat', 'TempStruct', '-v7.3')
    save('ilData.mat', 'OStruct', '-append')
    save('ilData.mat', 'N2Struct', '-append')
    save('ilData.mat', 'HeStruct', '-append')
    save('ilData.mat', 'ArStruct', '-append')
    save('ilData.mat', 'O2Struct', '-append')
    save('ilData.mat', 'rhoStruct', '-append')
    save('ilData.mat', 'lbDTStruct', '-append')
    save('ilData.mat', 'lbT0Struct', '-append')
end

if ~exist('goceVariables.mat', 'file')
    fprintf('%s\n', 'Saving results to "goceVariables.mat" file')
    save('goceVariables.mat', 'ae', '-v7.3')
    save('goceVariables.mat', 'aeIntegrals', '-append')
    save('goceVariables.mat', 'ap', '-append')
    save('goceVariables.mat', 'absB', '-append')
    save('goceVariables.mat', 'vBz', '-append')
    save('goceVariables.mat', 'densityIndex', '-append')
    save('goceVariables.mat', 'densityIndex1min', '-append')
    save('goceVariables.mat', 'densityIndexNoBg', '-append')
    save('goceVariables.mat', 'averagedDensity', '-append')
    save('goceVariables.mat', 'averagedDensityNoBg', '-append')
    save('goceVariables.mat', 'density3h', '-append')
    save('goceVariables.mat', 'akasofuEpsilon', '-append')
    save('goceVariables.mat', 'epsilonQualityFlag', '-append')
    save('goceVariables.mat', 'longitude', '-append')
    save('goceVariables.mat', 'latitude', '-append')
    save('goceVariables.mat', 'altitude', '-append')
    save('goceVariables.mat', 'solarTime', '-append')
    save('goceVariables.mat', 'crwindEast', '-append')
    save('goceVariables.mat', 'crwindNorth', '-append')
    save('goceVariables.mat', 'crwindUp', '-append')
    save('goceVariables.mat', 'magneticLatitude', '-append')
    save('goceVariables.mat', 'magneticLongitude', '-append')
    save('goceVariables.mat', 'magneticLocalTime', '-append')
    save('goceVariables.mat', 'densityNoBg', '-append')
    save('goceVariables.mat', 'msisDensityVariableAlt', '-append')
    save('goceVariables.mat', 'msisDensity270km', '-append')
    save('goceVariables.mat', 'msisDensity270kmNoAp', '-append')
    save('goceVariables.mat', 'dtm2013Density270km', '-append')
    save('goceVariables.mat', 'dtm2013Density270kmNoAm', '-append')
    save('goceVariables.mat', 'dtm2013DensityVariableAlt', '-append')
    save('goceVariables.mat', 'jb2008Density270km', '-append')
    save('goceVariables.mat', 'jb2008Density270kmNoDtc', '-append')
    save('goceVariables.mat', 'jb2008DensityVariableAlt', '-append')
    save('goceVariables.mat', 'tiegcmDensity270km', '-append')
    save('goceVariables.mat', 'tiegcmDensityVariableAlt', '-append')
    save('goceVariables.mat', 'hwmU', '-append')
    save('goceVariables.mat', 'hwmV', '-append')
    save('goceVariables.mat', 'morningFourierGrid', '-append')
    save('goceVariables.mat', 'eveningFourierGrid', '-append')
    save('goceVariables.mat', 'morningLatDoy', '-append')
    save('goceVariables.mat', 'eveningLatDoy', '-append')
    save('goceVariables.mat', 'morningBins', '-append')
    save('goceVariables.mat', 'eveningBins', '-append')
    save('goceVariables.mat', 'measuredDensity', '-append')
    save('goceVariables.mat', 'absDensityError', '-append')
    save('goceVariables.mat', 'absWindError', '-append')
    save('goceVariables.mat', 'noiseAffected', '-append')
    save('goceVariables.mat', 'eclipseAffected', '-append')
    save('goceVariables.mat', 'isMorningPass', '-append')
    save('goceVariables.mat', 'ionThrusterActive', '-append')
    save('goceVariables.mat', 'timestamps10sFixed', '-append')
    save('goceVariables.mat', 'timestamps1min', '-append')
    save('goceVariables.mat', 'timestamps1minFixed', '-append')
    save('goceVariables.mat', 'timestamps3h', '-append')
    save('goceVariables.mat', 'timestamps3hFixed', '-append')
    save('goceVariables.mat', 'timestampsAbsB', '-append')
    save('goceVariables.mat', 'timestampsEpsilon', '-append')
    save('goceVariables.mat', 'firstDatenum', '-append')
    save('goceVariables.mat', 'timestampsDensityDatenum', '-append')
    save('goceVariables.mat', 'doy', '-append')
    save('goceVariables.mat', 'timestampsAeDatenum', '-append')
    save('goceVariables.mat', 'timestampsEpsilonDatenum', '-append')
end

toc;
end

function [TempStruct, lbDTStruct, lbT0Struct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct] = arrangeByVar(goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberDT, T0)

TempStruct = struct('data', [], 'timestamps', [], 'longitude', [], 'latitude', [], 'altitude', [], 'solarTime', [], 'aeInt', zeros(0,9),...
                    'F', [], 'FA', [], 'apNow', [], 'ap3h', [], 'ap6h', [],'ap9h', [], 'ap12To33h', [], 'ap36To57h', [], 'Ap', []);

[tempTime, order] = unique(aeCData.temp.TTimes); data = aeCData.temp.T(order); aeC = 1:length(data); TempStruct = fillPosAndInd(TempStruct, aeCData, tempTime(aeC));
[t, order] = unique(aeEData.temp.TTimes); tempTime = [tempTime; t]; data = [data; aeEData.temp.T(order)]; aeE = aeC(end)+1 : length(data); TempStruct = fillPosAndInd(TempStruct, aeEData, tempTime(aeE));
[t, order] = unique(de2Data.temp.TTimes); tempTime = [tempTime; t]; data = [data; de2Data.temp.T(order)]; de2 = aeE(end)+1 : length(data); TempStruct = fillPosAndInd(TempStruct, de2Data, tempTime(de2));
TempStruct.data = data; TempStruct.timestamps = tempTime; TempStruct.aeC = aeC; TempStruct.aeE = aeE; TempStruct.de2 = de2;

OStruct = struct('data', [], 'timestamps', [], 'longitude', [], 'latitude', [], 'altitude', [], 'solarTime', [], 'aeInt', zeros(0,9),...
                    'F', [], 'FA', [], 'apNow', [], 'ap3h', [], 'ap6h', [],'ap9h', [], 'ap12To33h', [], 'ap36To57h', [], 'Ap', []);

[Otime, order] = unique(aeCData.oss.OTimes); data = aeCData.oss.O(order); aeC = 1:length(data); OStruct = fillPosAndInd(OStruct, aeCData, Otime(aeC));
[t, order] = unique(aeEData.nace.OTimes); Otime = [Otime; t]; data = [data; aeEData.nace.O(order)]; aeENace = aeC(end)+1 : length(data); OStruct = fillPosAndInd(OStruct, aeEData, Otime(aeENace));
[t, order] = unique(aeEData.oss.OTimes); Otime = [Otime; t]; data = [data; aeEData.oss.O(order)]; aeEOss = aeENace(end)+1 : length(data); OStruct = fillPosAndInd(OStruct, aeEData, Otime(aeEOss));
[t, order] = unique(de2Data.dens.OTimes); Otime = [Otime; t]; data = [data; de2Data.dens.O(order)]; de2 = aeEOss(end)+1 : length(data); OStruct = fillPosAndInd(OStruct, de2Data, Otime(de2));
OStruct.data = data; OStruct.timestamps = Otime; OStruct.aeC = aeC; OStruct.aeENace = aeENace; OStruct.aeEOss = aeEOss; OStruct.de2 = de2;

N2Struct = struct('data', [], 'timestamps', [], 'longitude', [], 'latitude', [], 'altitude', [], 'solarTime', [], 'aeInt', zeros(0,9),...
                    'F', [], 'FA', [], 'apNow', [], 'ap3h', [], 'ap6h', [],'ap9h', [], 'ap12To33h', [], 'ap36To57h', [], 'Ap', []);

[N2time, order] = unique(aeCData.oss.N2Times); data = aeCData.oss.N2(order); aeC = 1:length(data); N2Struct = fillPosAndInd(N2Struct, aeCData, N2time(aeC));
[t, order] = unique(aeEData.nace.N2Times); N2time = [N2time; t]; data = [data; aeEData.nace.N2(order)]; aeENace = aeC(end)+1 : length(data); N2Struct = fillPosAndInd(N2Struct, aeEData, N2time(aeENace));
[t, order] = unique(aeEData.oss.N2Times); N2time = [N2time; t]; data = [data; aeEData.oss.N2(order)]; aeEOss = aeENace(end)+1 : length(data); N2Struct = fillPosAndInd(N2Struct, aeEData, N2time(aeEOss));
[t, order] = unique(de2Data.dens.N2Times); N2time = [N2time; t]; data = [data; de2Data.dens.N2(order)]; de2 = aeEOss(end)+1 : length(data); N2Struct = fillPosAndInd(N2Struct, de2Data, N2time(de2));
[t, order] = unique(aerosData.dens.N2Times); N2time = [N2time; t]; data = [data; aerosData.dens.N2(order)]; aeros = de2(end)+1 : length(data); N2Struct = fillPosAndInd(N2Struct, aerosData, N2time(aeros));
N2Struct.data = data; N2Struct.timestamps = N2time; N2Struct.aeC = aeC; N2Struct.aeENace = aeENace; N2Struct.aeEOss = aeEOss; N2Struct.de2 = de2; N2Struct.aeros = aeros;

HeStruct = struct('data', [], 'timestamps', [], 'longitude', [], 'latitude', [], 'altitude', [], 'solarTime', [], 'aeInt', zeros(0,9),...
                    'F', [], 'FA', [], 'apNow', [], 'ap3h', [], 'ap6h', [],'ap9h', [], 'ap12To33h', [], 'ap36To57h', [], 'Ap', []);

[Hetime, order] = unique(aeCData.oss.HeTimes); data = aeCData.oss.He(order); aeC = 1:length(data); HeStruct = fillPosAndInd(HeStruct, aeCData, Hetime(aeC));
[t, order] = unique(aeEData.nace.HeTimes); Hetime = [Hetime; t]; data = [data; aeEData.nace.He(order)]; aeENace = aeC(end)+1 : length(data); HeStruct = fillPosAndInd(HeStruct, aeEData, Hetime(aeENace));
[t, order] = unique(aeEData.oss.HeTimes); Hetime = [Hetime; t]; data = [data; aeEData.oss.He(order)]; aeEOss = aeENace(end)+1 : length(data); HeStruct = fillPosAndInd(HeStruct, aeEData, Hetime(aeEOss));
[t, order] = unique(de2Data.dens.HeTimes); Hetime = [Hetime; t]; data = [data; de2Data.dens.He(order)]; de2 = aeEOss(end)+1 : length(data); HeStruct = fillPosAndInd(HeStruct, de2Data, Hetime(de2));
[t, order] = unique(aerosData.dens.HeTimes); Hetime = [Hetime; t]; data = [data; aerosData.dens.He(order)]; aeros = de2(end)+1 : length(data); HeStruct = fillPosAndInd(HeStruct, aerosData, Hetime(aeros));
HeStruct.data = data; HeStruct.timestamps = Hetime; HeStruct.aeC = aeC; HeStruct.aeENace = aeENace; HeStruct.aeEOss = aeEOss; HeStruct.de2 = de2; HeStruct.aeros = aeros;

ArStruct = struct('data', [], 'timestamps', [], 'longitude', [], 'latitude', [], 'altitude', [], 'solarTime', [], 'aeInt', zeros(0,9),...
                    'F', [], 'FA', [], 'apNow', [], 'ap3h', [], 'ap6h', [],'ap9h', [], 'ap12To33h', [], 'ap36To57h', [], 'Ap', []);

[Artime, order] = unique(aeCData.oss.ArTimes); data = aeCData.oss.Ar(order); aeCOss = 1:length(data); ArStruct = fillPosAndInd(ArStruct, aeCData, Artime(aeCOss));
[t, order] = unique(aeEData.oss.ArTimes); Artime = [Artime; t]; data = [data; aeEData.oss.Ar(order)]; aeEOss = aeCOss(end)+1 : length(data); ArStruct = fillPosAndInd(ArStruct, aeEData, Artime(aeEOss));
[t, order] = unique(de2Data.dens.ArTimes); Artime = [Artime; t]; data = [data; de2Data.dens.Ar(order)]; de2 = aeEOss(end)+1 : length(data); ArStruct = fillPosAndInd(ArStruct, de2Data, Artime(de2));
[t, order] = unique(aerosData.dens.ArTimes); Artime = [Artime; t]; data = [data; aerosData.dens.Ar(order)]; aeros = de2(end)+1 : length(data); ArStruct = fillPosAndInd(ArStruct, aerosData, Artime(aeros));
ArStruct.data = data; ArStruct.timestamps = Artime; ArStruct.aeCOss = aeCOss; ArStruct.aeEOss = aeEOss; ArStruct.de2 = de2; ArStruct.aeros = aeros;

O2Struct = struct('data', [], 'timestamps', [], 'longitude', [], 'latitude', [], 'altitude', [], 'solarTime', [], 'aeInt', zeros(0,9),...
                    'F', [], 'FA', [], 'apNow', [], 'ap3h', [], 'ap6h', [],'ap9h', [], 'ap12To33h', [], 'ap36To57h', [], 'Ap', []);

[O2time, order] = unique(aeCData.oss.O2Times); data = aeCData.oss.O2(order); aeCOss = 1:length(data); O2Struct = fillPosAndInd(O2Struct, aeCData, O2time(aeCOss));
[t, order] = unique(aeEData.oss.O2Times); O2time = [O2time; t]; data = [data; aeEData.oss.O2(order)]; aeEOss = aeCOss(end)+1 : length(data); O2Struct = fillPosAndInd(O2Struct, aeEData, O2time(aeEOss));
O2Struct.data = data; O2Struct.timestamps = O2time; O2Struct.aeCOss = aeCOss; O2Struct.aeEOss = aeEOss; O2Struct.de2 = de2; O2Struct.aeros = aeros;

lbDTStruct = saberDT;
lbT0Struct = T0;

rhoStruct = struct('data', [], 'timestamps', [], 'longitude', [], 'latitude', [], 'altitude', [], 'solarTime', [], 'aeInt', zeros(0,9),...
                    'F', [], 'FA', [], 'apNow', [], 'ap3h', [], 'ap6h', [],'ap9h', [], 'ap12To33h', [], 'ap36To57h', [], 'Ap', []);

combStruct = [goceData; champData; graceData];
rhoStruct.data = vertcat(combStruct.density);
rhoStruct.timestamps = vertcat(combStruct.timestamps);
rhoStruct.longitude = vertcat(combStruct.longitude);
rhoStruct.latitude = vertcat(combStruct.latitude);
rhoStruct.altitude = vertcat(combStruct.altitude);
rhoStruct.solarTime = vertcat(combStruct.solarTime);
rhoStruct.aeInt = vertcat(combStruct.aeInt);
rhoStruct.F = vertcat(combStruct.F10);
rhoStruct.FA = vertcat(combStruct.F81A);
rhoStruct.apNow = vertcat(combStruct.apNow);
rhoStruct.ap3h = vertcat(combStruct.ap3h);
rhoStruct.ap6h = vertcat(combStruct.ap6h);
rhoStruct.ap9h = vertcat(combStruct.ap9h);
rhoStruct.ap12To33h = vertcat(combStruct.apAver12To33h);
rhoStruct.ap36To57h = vertcat(combStruct.apAver36To57h);
rhoStruct.Ap = vertcat(combStruct.ApDaily);
k = length(goceData.timestamps); rhoStruct.goce = 1:k;
rhoStruct.champ = k+1 : k + length(champData.timestamps); k = k + length(champData.timestamps);
rhoStruct.grace = k+1 : k + length(graceData.timestamps);

end

function [fillStruct] = fillPosAndInd(fillStruct, satDataStruct, obsTimes)
    [satTimes, order] = unique(satDataStruct.timestamps);
    satLat = satDataStruct.latitude(order);
    satLon = satDataStruct.longitude(order);
    satLst = satDataStruct.solarTime(order);
    satAlt = satDataStruct.altitude(order);
    satAeInt = satDataStruct.aeInt(order,:);
    satF = satDataStruct.F10(order);
    satFA = satDataStruct.F81A(order);
    satApNow = satDataStruct.apNow(order);
    satAp3h = satDataStruct.ap3h(order);
    satAp6h = satDataStruct.ap6h(order);
    satAp9h = satDataStruct.ap9h(order);
    satAp12To33h = satDataStruct.apAver12To33h(order);
    satAp36To57h = satDataStruct.apAver36To57h(order);
    satAp = satDataStruct.ApDaily(order);
    
    ind = ismember(satTimes, obsTimes);
    fillStruct.latitude = [fillStruct.latitude; satLat(ind)];
    fillStruct.longitude = [fillStruct.longitude; satLon(ind)];
    fillStruct.solarTime = [fillStruct.solarTime; satLst(ind)];
    fillStruct.altitude = [fillStruct.altitude; satAlt(ind)];
    fillStruct.aeInt = [fillStruct.aeInt; satAeInt(ind, :)];
    fillStruct.F = [fillStruct.F; satF(ind)];
    fillStruct.FA = [fillStruct.FA; satFA(ind)];
    fillStruct.apNow = [fillStruct.apNow; satApNow(ind)];
    fillStruct.ap3h = [fillStruct.ap3h; satAp3h(ind)];
    fillStruct.ap6h = [fillStruct.ap6h; satAp6h(ind)];
    fillStruct.ap9h = [fillStruct.ap9h; satAp9h(ind)];
    fillStruct.ap12To33h = [fillStruct.ap12To33h; satAp12To33h(ind)];
    fillStruct.ap36To57h = [fillStruct.ap36To57h; satAp36To57h(ind)];
    fillStruct.Ap = [fillStruct.Ap; satAp(ind)];
end

function [ae, timestampsAeDatenum] = readAeFiles()
%
fprintf('%s\n', 'Began reading AE files')

ae = [];
timestampsAeDatenum = [];

aeFiles = [dir('ae2*'); dir('ae198*'); dir('ae1975'); dir('ae1978'); dir('ae1979')];
parfor i = 1:length(aeFiles)
    aeFile = fopen(aeFiles(i).name);
    if aeFile == -1
        error('aeFile open unsuccesful')
    end

    aeData = textscan(aeFile, '%s %s %f %f %f %f %f', 'MultipleDelimsAsOne',1, 'HeaderLines',15);

    timestampsAeDatenum = [timestampsAeDatenum; datenum(strcat(aeData{1}, aeData{2}), 'yyyy-mm-ddHH:MM:SS.FFF')];
    ae = [ae; aeData{4}];
end

aeFiles = [dir('ae1972*'); dir('ae1973*'); dir('ae1974*')];
parfor i = 1:length(aeFiles)
    aeFile = fopen(aeFiles(i).name);
    if aeFile == -1
        error('aeFile open unsuccesful')
    end
    
    aeData = textscan(aeFile, '%f %f %f %f %f %f %f %f', 'MultipleDelimsAsOne',1, 'HeaderLines',15);
    yr = aeData{1}; mon = aeData{2}; d = aeData{3};
    hh = floor(aeData{4} / 10000); minute = floor((aeData{4} - hh*10000) / 100); sec = aeData{4} - hh*10000 - minute*100;
    timestampsAeDatenum = [timestampsAeDatenum; datenum([yr, mon, d, hh, minute, sec])];
    ae = [ae; aeData{5}];
end

[timestampsAeDatenum, indicesToConserve, ~] = unique(timestampsAeDatenum);
ae = ae(indicesToConserve);

interpIndices = ae > 50000;
aeInterp = ae(~interpIndices);
tInterp = timestampsAeDatenum(~interpIndices);
ae = interp1(tInterp, aeInterp, timestampsAeDatenum, 'linear', 'extrap');

end

function [apGoce, ap, timestampsAp, am, dtc] = readApAndDtcFiles(timestampsDensityDatenum)
%
fprintf('%s\n', 'Began reading ap file')

apFile = fopen('apData');
if apFile == -1
    error('ap File open unsuccesful! Check that you have a file "apData" in your WORKING DIRECTORY.')
end

apData = textscan(apFile, '%s %s %s %s %s %s %s %s %s %s %f %f %f %f %f %f %f %f %f', ...
    'delimiter', {' ', '+', '-'}, 'MultipleDelimsAsOne',1, 'HeaderLines',1);
apValues = horzcat(apData{11:18});

[row, ~] = find(isnan(apValues) | apValues >= 1000);
row = unique(row);
for i = 1:length(row)
    frewind(apFile);
    nanLine = textscan(apFile, '%s', 1, 'delimiter', '\n', 'HeaderLines', row(i));
    nanLine = nanLine{1}(1);
    correctValue = nan(52, 1);
    for k = 31:3:52
     correctValue(k) = str2double(horzcat(nanLine{1}(k-2), nanLine{1}(k-1), nanLine{1}(k)));
    end
    correctValue = correctValue(~isnan(correctValue));
    apValues(row(i), :) = correctValue;
end

timestampsAp = datenum(apData{1}, 'yyyymmdd');
apRows = find(timestampsAp > timestampsDensityDatenum(1) - 4 & timestampsAp <= timestampsDensityDatenum(end));
apGoce = reshape(apValues(apRows, :)', numel(apValues(apRows, :)), 1);

timestampsAp = repmat(timestampsAp', 8, 1); 
timestampsAp = timestampsAp + repmat((3:3:24)'/24, 1, size(timestampsAp, 2));
timestampsAp = timestampsAp(:);
apValues = apValues'; ap = apValues(:);

firstHour = round(str2double(datestr(timestampsDensityDatenum(1), 'HH')));
apGoce(1:firstHour/3) = [];


amFile = fopen('am_file_spider.dat');
if amFile == -1
    error('am File open unsuccesful! Check that you have a file "am_file_spider.data" in your WORKING DIRECTORY.')
end

amData = textscan(amFile, '%f %f %f %f %f %f %f %f %f %f %f %f %f', 'MultipleDelimsAsOne',1);
amValues = horzcat(amData{end-8:end});
[row, ~] = find(isnan(amValues));
amValues(row,2:9) = amValues(row,1:8);
amValues = amValues(:,2:9);

year = mod(amData{2}, 10000);
month = amData{1};
day = floor(amData{2}/10000);
timestampsAm = datenum(year, month, day);
amValues(amValues < 0) = 0.0;
amRows = find(timestampsAm > timestampsDensityDatenum(1) - 4 & timestampsAm <= timestampsDensityDatenum(end));
am = reshape(amValues(amRows, :)', numel(amValues(amRows, :)), 1);
am(1:firstHour/3) = [];


dtcFile = fopen('DTCFILE.TXT');
if dtcFile == -1
    error('Dtc file open unsuccesful! Check that you have a file "DTCFILE.TXT" in your WORKING DIRECTORY.')
end

dtcData = textscan(dtcFile, '%s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d', 'MultipleDelimsAsOne',1, 'EmptyValue', 9999);
dtcData = horzcat(dtcData{2:end});
row = find(dtcData(:,end) > 999);
for i = 1:length(row)
    frewind(dtcFile);
    nanLine = textscan(dtcFile, '%s', 1, 'delimiter', '\n', 'HeaderLines', row(i) - 1);
    nanLine = nanLine{1}(1);
    yy = str2double(nanLine{1}(5:6));
    ddd = str2double(nanLine{1}(7:9));
    dtcData(row(i), :) = [yy ddd dtcData(row(i), 2:end-1)];
end

year = dtcData(:,1);
lastCentury = year >= 97;
year(lastCentury) = year(lastCentury) + 1900;
year(~lastCentury) = year(~lastCentury) + 2000;
dayOfYear = dtcData(:,2);
dtcDatenum = datenum(double(year), 1, 1) + double(dayOfYear) - 1;
dtcRows = find(dtcDatenum > timestampsDensityDatenum(1) - 1 & dtcDatenum <= timestampsDensityDatenum(end));
dtcValues = dtcData(:,3:end);
dtc = reshape(dtcValues(dtcRows, :)', numel(dtcValues(dtcRows, :)), 1);

numOfDays = length(dtcRows);
firstDtcDay = dtcDatenum(dtcRows(1));
timestampsDtc = (0 : 3600 : numOfDays * 86400 - 1)';
timestampsFixed = round((timestampsDensityDatenum - firstDtcDay) * 86400);
dtc = interp1(timestampsDtc, double(dtc), timestampsFixed, 'linear', 'extrap');

end

function [F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, F10datenum, datenumToJulian] = readSolarIndexFiles
%
fprintf('%s\n', 'Began reading solar index files')

solarFile = fopen('SOLFSMY.TXT');
if solarFile == -1
    error('Solar index file open unsuccesful! Check that you have a file "SOLFSMY.TXT" in your WORKING DIRECTORY.')
end

solarData = textscan(solarFile, '%d %d %f %f %f %f %f %f %f %f %f %s', 'MultipleDelimsAsOne',1, 'CommentStyle','#');
%F10 = solarData{4};
%F81A = solarData{5};
S10 = solarData{6};
S81A = solarData{7};
M10 = solarData{8};
M81A = solarData{9};
Y10 = solarData{10};
Y81A = solarData{11};

SnoNans = find(S10 > 1);
S81noNans = find(S81A > 1);
S10 = interp1(SnoNans, S10(SnoNans), 1:length(S10), 'linear', 'extrap')';
S81A = interp1(S81noNans, S81A(S81noNans), 1:length(S81A), 'linear', 'extrap')';

year = solarData{1};
dayOfYear = solarData{2};
indexDatenums = datenum(double(year), 1, 1) + double(dayOfYear) - 1;
julianDate = solarData{3};
datenumToJulian = julianDate(1) - indexDatenums(1) - 0.5;


solarFile = fopen('proxies_unadjusted.dat');
if solarFile == -1
    error('Solar index file open unsuccesful! Check that you have a file "proxies_unadjusted.dat" in your WORKING DIRECTORY.')
end

solarData = textscan(solarFile, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'MultipleDelimsAsOne',1, 'CommentStyle','//');
F30 = solarData{14};
F30A = solarData{15};
F10 = solarData{5};
F81A = solarData{6};
F30toF10 = mean(F10 ./ F30);

year = solarData{1};
month = solarData{2};
day = solarData{3};
F10datenum = datenum(year, month, day);
F30 = F30toF10 * F30; %interp1(F10datenum, F30, indexDatenums, 'nearest', 'extrap');
F30A = F30toF10 * F30A;% interp1(F10datenum, F30A, indexDatenums, 'nearest', 'extrap');

end

function [absB, timestampsAbsBDatenum, akasofuEpsilon, epsilonQualityFlag, timestampsEpsilonDatenum, vBz] = ...
    readSolarWindFiles(timestampsAeDatenum)
%
fprintf('%s\n', 'Began reading IMF and SW files')

imfFiles = dir('ac_h0s_mfi*');

if isempty(imfFiles)
    error('No IMF B files found! Make sure there are files beginning "ac_h0s_mfi_*" in the WORKING DIRECTORY')
end

timestampsBDatenum = [];
Bx = [];
By = [];
Bz = [];
for i = 1:length(imfFiles)
     data = cdfread(imfFiles(i).name, 'CombineRecords',true, 'ConvertEpochToDatenum',true);
     timestampsBDatenum = [timestampsBDatenum; data{1}];
     Bx = [Bx; data{2}(1:3:end)'];
     By = [By; data{2}(2:3:end)'];
     Bz = [Bz; data{2}(3:3:end)'];
end
[timestampsBDatenum, indicesToConserve, ~] = unique(timestampsBDatenum);
Bx = Bx(indicesToConserve);
By = By(indicesToConserve);
Bz = Bz(indicesToConserve);
B2 = Bx.^2 + By.^2 + Bz.^2;
theta = atan2(By, Bz);

indicesToConserve = find(Bx > -1e30 & By > -1e30 & Bz > -1e30);

tBNoNans = timestampsBDatenum(indicesToConserve);
thetaNoNans = theta(indicesToConserve);
B2NoNans = B2(indicesToConserve);
BzNoNans = Bz(indicesToConserve);

tBToInterpolate = timestampsAeDatenum(timestampsAeDatenum>=min(timestampsBDatenum) & timestampsAeDatenum<=max(timestampsBDatenum));
tThetaToInterpolate = tBToInterpolate;
B2Interpolated = interp1(tBNoNans, B2NoNans, tBToInterpolate, 'linear', 'extrap');
BzInterpolated = interp1(tBNoNans, BzNoNans, tBToInterpolate, 'linear', 'extrap');
thetaInterpolated = interp1(tBNoNans, thetaNoNans, tThetaToInterpolate, 'linear', 'extrap');

secondsInDay = 24 * 60 * 60;
absB = sqrt(B2Interpolated);
timestampsAbsBDatenum = tBToInterpolate;

solarWindFiles = dir('ac_h0s_swe*');

if isempty(solarWindFiles)
    error('No Solar Wind files found! Make sure there are files beginning "ac_h0s_swe_*" in the WORKING DIRECTORY')
end

timestampsVDatenum = [];
Vx = [];
Vy = [];
Vz = [];
for i = 1:length(solarWindFiles)
     data = cdfread(solarWindFiles(i).name, 'CombineRecords',true, 'ConvertEpochToDatenum',true);
     timestampsVDatenum = [timestampsVDatenum; data{1}];
     Vx = [Vx; data{2}(1:3:end)'];
     Vy = [Vy; data{2}(2:3:end)'];
     Vz = [Vz; data{2}(3:3:end)'];
end
[timestampsVDatenum, indicesToConserve, ~] = unique(timestampsVDatenum);
Vx = Vx(indicesToConserve);
Vy = Vy(indicesToConserve);
Vz = Vz(indicesToConserve);
V = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
indicesToConserve = find(Vx > -1e30 & Vy > -1e30 & Vz > -1e30);
tVNoNans = timestampsVDatenum(indicesToConserve);
VNoNans = V(indicesToConserve);

tVToInterpolate = timestampsAeDatenum(timestampsAeDatenum>=min(timestampsVDatenum) & timestampsAeDatenum<=max(timestampsVDatenum));
VInterpolated = interp1(tVNoNans, VNoNans, tVToInterpolate, 'linear', 'extrap');
% timesToNeareastNonNan = arrayfun(@(x) min(abs(tVNoNans - x)), tVToInterpolate);
% epsilonQualityFlag = ones(length(tVToInterpolate), 1);
% epsilonQualityFlag(timesToNeareastNonNan > oneMinuteInDays) = 0;  
commonIndices = timestampsVDatenum >= min(tVToInterpolate) & timestampsVDatenum <= max(tVToInterpolate);
badDataIndices = ~ismember(timestampsVDatenum(commonIndices), tVNoNans);
epsilonQualityFlag = ones(length(tVToInterpolate), 1);
epsilonQualityFlag(badDataIndices) = 0; 

VInterpolated = VInterpolated(ismember(tVToInterpolate, tBToInterpolate)) * 1000;
epsilonQualityFlag = epsilonQualityFlag(ismember(tVToInterpolate, tBToInterpolate));
B2Interpolated = B2Interpolated(ismember(tBToInterpolate, tVToInterpolate)) * 1e-9;
BzInterpolated = BzInterpolated(ismember(tBToInterpolate, tVToInterpolate)) * 1e-9;
thetaInterpolated = thetaInterpolated(ismember(tThetaToInterpolate, tVToInterpolate));

akasofuEpsilon = VInterpolated .* B2Interpolated .* sin(thetaInterpolated * 0.5).^4;
vBz = VInterpolated .* BzInterpolated;

tVToInterpolate = tVToInterpolate(ismember(tVToInterpolate, tBToInterpolate));
timestampsEpsilonDatenum = tVToInterpolate;

end

function [density, longitude, latitude, altitude, solarTime, magneticLatitude, magneticLongitude, magneticLocalTime, crwindEast, crwindNorth, crwindUp,...
    timestamps10sFixed, timestamps1min, timestamps1minFixed, timestampsDensityDatenum, doy, timestampsAbsB, timestampsEpsilon, firstDatenum,...
    absDensityError, absWindError, noiseAffected, eclipseAffected, isMorningPass, ionThrusterActive] = ...
    readDensityFile(timestampsAeDatenum, timestampsAbsBDatenum, timestampsEpsilonDatenum)
%
fprintf('%s\n', 'Began reading GOCE files')

densityFiles = dir('goce_dens*');

density = [];
longitude = [];
latitude = [];
altitude = [];
solarTime = [];
crwindEast = [];
crwindNorth = [];
crwindUp = [];
absDensityError = [];
absWindError = [];
noiseAffected = [];
eclipseAffected = [];
isMorningPass = [];
ionThrusterActive = [];

timestampsDensityDatenum = [];
doy = [];

secondsInDay = 60 * 60 * 24;
parfor i = 1:length(densityFiles)
    densityFile = fopen(densityFiles(i).name);

    if densityFile == -1
        error('densityFile open unsuccesful')
    end

    densityData = textscan(densityFile, '%s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'commentStyle','#');

    density = [density; densityData{9} * power(10, 11)];
    longitude = [longitude; densityData{5}];
    latitude = [latitude; densityData{6}];
    altitude = [altitude; densityData{4}];
    solarTime = [solarTime; densityData{7}];
    crwindEast = [crwindEast; densityData{10}];
    crwindNorth = [crwindNorth; densityData{11}];
    crwindUp = [crwindUp; densityData{12}];
    absDensityError = [absDensityError; densityData{13} * power(10, 11)];
    absWindError = [absWindError; densityData{14}];
    noiseAffected = [noiseAffected; densityData{15}];
    eclipseAffected = [eclipseAffected; densityData{16}];
    isMorningPass = [isMorningPass; densityData{17}];
    ionThrusterActive = [ionThrusterActive; densityData{18}];

    timestampsThisFile = datenum(strcat(densityData{1}, densityData{2}), 'yyyy-mm-ddHH:MM:SS.FFF');
    timestampsDensityDatenum = [timestampsDensityDatenum; timestampsThisFile];
    
    doyThisFile = timestampsThisFile - datenum(datestr(timestampsThisFile, 'yyyy'), 'yyyy');
    doy = [doy; doyThisFile];
end

indicesToRemove = roundTimestampsToBeginFromNearestThreeHourMark(timestampsDensityDatenum);

[~, uniqueIndices, ~] = unique(timestampsDensityDatenum);
uniqueIndices = ismember(1:length(timestampsDensityDatenum), uniqueIndices)';
indicesToConserve = ~indicesToRemove & uniqueIndices;

timestampsDensityDatenum = timestampsDensityDatenum(indicesToConserve);
doy = doy(indicesToConserve);
density = density(indicesToConserve);
longitude = longitude(indicesToConserve);
latitude = latitude(indicesToConserve);
altitude = altitude(indicesToConserve);
solarTime = solarTime(indicesToConserve);
crwindEast = crwindEast(indicesToConserve);
crwindNorth = crwindNorth(indicesToConserve);
crwindUp = crwindUp(indicesToConserve);
absDensityError = absDensityError(indicesToConserve);
absWindError = absWindError(indicesToConserve);
noiseAffected = logical(round(noiseAffected(indicesToConserve)));
eclipseAffected = logical(round(eclipseAffected(indicesToConserve)));
isMorningPass = logical(round(isMorningPass(indicesToConserve)));
ionThrusterActive = ~logical(round(ionThrusterActive(indicesToConserve)));

[magneticLatitude, magneticLongitude] = convertToMagneticCoordinates(latitude, longitude, altitude);

timestamps10sFixed = timestampsDensityDatenum * secondsInDay;
timestamps1minFixed = timestamps10sFixed(ismember(timestampsDensityDatenum, timestampsAeDatenum));
firstDatenum = datenum(datestr(timestampsDensityDatenum(1), 'yyyy-mm-dd'), 'yyyy-mm-dd');
timestamps10sFixed = timestamps10sFixed - firstDatenum * secondsInDay;
timestamps1minFixed = timestamps1minFixed - firstDatenum * secondsInDay;
timestamps10sFixed = round(timestamps10sFixed);
timestamps1minFixed = round(timestamps1minFixed);
timestamps1min = round((timestampsAeDatenum - firstDatenum) * secondsInDay);

timestampsAbsB = round((timestampsAbsBDatenum - firstDatenum) * secondsInDay);
timestampsEpsilon = round((timestampsEpsilonDatenum - firstDatenum) * secondsInDay);

[magneticLocalTime] = computeMagneticTime(magneticLongitude, doy, timestampsDensityDatenum);

end

function [champData, graceData, de2Data, aeCData, aeEData, aerosData] = computeOtherDensities()

fprintf('%s\n', 'Began reading CHAMP files')

champFiles1 = dir('champ/Density_3deg_0*');
champFiles2 = dir('champ/Density_10_*');

champDensity = [];
champError = [];
champTimestamps = [];
champLatitude = [];
champLongitude = [];
champLst = [];
champAltitude = [];
for i = 1:length(champFiles1)
    load(['champ/', champFiles1(i).name])
    yearVec = repmat([Year.data+2000,1,1,0,0,0], length(Sec.data), 1);
    thisFileTimestamps = datenum(yearVec) + Doy.data + Sec.data/86400 - 1;
    champTimestamps = [champTimestamps; thisFileTimestamps];
    champDensity = [champDensity; Density.data];
    champError = [champError; U_rho.data];
    champLatitude = [champLatitude; Lat.data];
    champLongitude = [champLongitude; Lon.data];
    champAltitude = [champAltitude; Height.data];
    champLst = [champLst; LocTim.data];
end
for i = 1:length(champFiles2)
    load(['champ/',champFiles2(i).name])
    yearVec = [loop.GPSyy'+2000, repmat([1,1,0,0,0], length(loop.GPSsec), 1)];
    thisFileTimestamps = datenum(yearVec) + loop.GPSdoy' + loop.GPSsec'/86400 - 1;
    champTimestamps = [champTimestamps; thisFileTimestamps];
    champDensity = [champDensity; loop.Dstar'];
    champLatitude = [champLatitude; loop.lat'];
    champLongitude = [champLongitude; loop.lon'];
    champAltitude = [champAltitude; loop.height'];
    champLst = [champLst; loop.slt'];
end
champ2010Error = ones(length(champTimestamps)-length(champError),1) * mean(champError);
champError = [champError; champ2010Error]; % !!!!!!!!!!!!!!!!!!!!!
[champTimestamps, order] = unique(champTimestamps);
champDensity = champDensity(order);
champError = champError(order);
champLatitude = champLatitude(order);
champLongitude = champLongitude(order);
champAltitude = champAltitude(order);
champLst = champLst(order);
champData = struct('density', champDensity, 'densityError', champError, 'timestamps', champTimestamps, 'latitude', champLatitude, 'longitude', champLongitude,...
                  'altitude', champAltitude, 'solarTime', champLst);

fprintf('%s\n', 'Began reading GRACE files')
graceFiles1 = dir('grace/Density_graceA_3deg_0*');
graceFiles2 = dir('grace/Density_graceA_10_*');

graceDensity = [];
graceError = [];
graceTimestamps = [];
graceLatitude = [];
graceLongitude = [];
graceLst = [];
graceAltitude = [];
for i = 1:length(graceFiles1)
    load(['grace/',graceFiles1(i).name])
    yearVec = repmat([Year.data+2000,1,1,0,0,0], length(Sec.data), 1);
    thisFileTimestamps = datenum(yearVec) + Doy.data + Sec.data/86400 - 1;
    graceTimestamps = [graceTimestamps; thisFileTimestamps];
    graceDensity = [graceDensity; Density.data];
    graceError = [graceError; U_rho.data];
    graceLatitude = [graceLatitude; Lat.data];
    graceLongitude = [graceLongitude; Lon.data];
    graceAltitude = [graceAltitude; Height.data];
    graceLst = [graceLst; LocTim.data];
end
for i = 1:length(graceFiles2)
    try
        load(['grace/',graceFiles2(i).name])
    catch
        continue
    end
    yearVec = [loop.GPSyy'+2000, repmat([1,1,0,0,0], length(loop.GPSsec'), 1)];
    thisFileTimestamps = datenum(yearVec) + loop.GPSdoy' + loop.GPSsec'/86400 - 1;
    graceTimestamps = [graceTimestamps; thisFileTimestamps];
    graceDensity = [graceDensity; loop.Dstar'];
    graceError = [graceError; Derror.U_rho'];
    graceLatitude = [graceLatitude; loop.lat'];
    graceLongitude = [graceLongitude; loop.lon'];
    graceAltitude = [graceAltitude; loop.height'];
    graceLst = [graceLst; loop.slt'];
end
[graceTimestamps, order] = unique(graceTimestamps);
graceDensity = graceDensity(order);
graceError = graceError(order);
graceLatitude = graceLatitude(order);
graceLongitude = graceLongitude(order);
graceAltitude = graceAltitude(order);
graceLst = graceLst(order);
graceData = struct('density', graceDensity, 'densityError', graceError, 'timestamps', graceTimestamps, 'latitude', graceLatitude, 'longitude', graceLongitude,...
                  'altitude', graceAltitude, 'solarTime', graceLst);

fprintf('%s\n', 'Began reading DE-2 files')       
de2Files = dir('ae_de2_*');
de2N2 = []; N2timestamps = [];
de2O = []; Otimestamps = [];
de2He = []; HEtimestamps = [];
de2Ar = []; ARtimestamps = [];
de2N = []; Ntimestamps = [];
de2Timestamps = [];
de2Latitude = [];
de2Longitude = [];
de2Lst = [];
de2Altitude = [];
for i = 1:length(de2Files)
    de2File = fopen(de2Files(i).name);
    data = textscan(de2File, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f','MultipleDelimsAsOne',1);
    thisFileTimestamps = datenum([data{1}+1900, data{2}, data{3}, data{4}, data{5}, data{6}]);
    de2Timestamps = [de2Timestamps; thisFileTimestamps];
    de2Altitude = [de2Altitude; data{7}];
    de2Latitude = [de2Latitude; data{8}];
    de2Longitude = [de2Longitude; data{9}];
    utHour = data{4} + data{5}/60 + data{6}/3600;
    lst = utHour + data{9} / 15;
    lst(lst >= 24) = lst(lst >= 24) - 24;
    lst(lst < 0) = lst(lst < 0) + 24;
    de2Lst = [de2Lst; lst];
    
    N2 = data{10}; N2timestamps = [N2timestamps; thisFileTimestamps(N2 < 1E30)]; N2(N2 > 1E30) = []; de2N2 = [de2N2; N2];
    O = data{11}; Otimestamps = [Otimestamps; thisFileTimestamps(O < 1E30)]; O(O > 1E30) = []; de2O = [de2O; O];
    He = data{12}; HEtimestamps = [HEtimestamps; thisFileTimestamps(He < 1E30)]; He(He > 1E30) = []; de2He = [de2He; He];
    Ar = data{13}; ARtimestamps = [ARtimestamps; thisFileTimestamps(Ar < 1E30)]; Ar(Ar > 1E30) = []; de2Ar = [de2Ar; Ar];
    N = data{14}; Ntimestamps = [Ntimestamps; thisFileTimestamps(N < 1E30)]; N(N > 1E30) = []; de2N = [de2N; N];
end
[de2Timestamps, order] = unique(de2Timestamps);
de2Latitude = de2Latitude(order);
de2Longitude = de2Longitude(order);
de2Altitude = de2Altitude(order);
de2Lst = de2Lst(order);

densities = struct('N2', de2N2, 'O', de2O, 'He', de2He, 'N', de2N, 'Ar', de2Ar, ...
                 'N2Times', N2timestamps, 'OTimes', Otimestamps, 'HeTimes', HEtimestamps, 'NTimes', Ntimestamps, 'ArTimes', ARtimestamps);

de2Data = struct('dens', densities, ...
                 'timestamps', de2Timestamps, ...
                 'latitude', de2Latitude, 'longitude', de2Longitude,...
                  'altitude', de2Altitude, 'solarTime', de2Lst);

fprintf('%s\n', 'Began reading AE-C/E files')
aeSatellite = 'c';
for i = 1:2
    if strcmp(aeSatellite,'c')
        aeFiles = dir('ae_c_*');
    else
        aeFiles = dir('ae_e_*');
    end
    
    aeTimestamps = [];
    
    naceN2 = []; naceN2timestamps = [];
    naceO = []; naceOtimestamps = [];
    naceHe = []; naceHEtimestamps = [];
    naceAr = []; naceARtimestamps = [];
    
    ossN2 = []; ossN2timestamps = [];
    ossO = []; ossOtimestamps = [];
    ossO2 = []; ossO2timestamps = [];
    ossHe = []; ossHEtimestamps = [];
    ossAr = []; ossARtimestamps = [];
    ossN = []; ossNtimestamps = [];
    
    aeLatitude = [];
    aeLongitude = [];
    aeLst = [];
    aeAltitude = [];
    for j = 1:length(aeFiles)
        aeFile = fopen(aeFiles(j).name);
        if strcmp(aeSatellite,'c')
            data = textscan(aeFile, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','MultipleDelimsAsOne',1);
        else
            data = textscan(aeFile, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','MultipleDelimsAsOne',1);
        end
        thisFileTimestamps = datenum([data{1}+1900, data{2}, data{3}, data{4}, data{5}, data{6}]);
        utHour = data{4} + data{5}/60 + data{6}/3600;
        lst = utHour + data{9} / 15;
        lst(lst >= 24) = lst(lst >= 24) - 24;
        lst(lst < 0) = lst(lst < 0) + 24;

        aeLst = [aeLst; lst];
        aeTimestamps = [aeTimestamps; thisFileTimestamps];
        aeAltitude = [aeAltitude; data{7}];
        aeLatitude = [aeLatitude; data{8}];
        aeLongitude = [aeLongitude; data{9}];
        
        if strcmp(aeSatellite,'c')
            % NACE
            N2 = data{11}; naceN2timestamps = [naceN2timestamps; thisFileTimestamps(N2 < 1E30)]; N2(N2 > 1E30) = []; naceN2 = [naceN2; N2];
            O = data{12}; naceOtimestamps = [naceOtimestamps; thisFileTimestamps(O < 1E30)]; O(O > 1E30) = []; naceO = [naceO; O];
            He = data{13}; naceHEtimestamps = [naceHEtimestamps; thisFileTimestamps(He < 1E30)]; He(He > 1E30) = []; naceHe = [naceHe; He];
            Ar = data{14}; naceARtimestamps = [naceARtimestamps; thisFileTimestamps(Ar < 1E30)]; Ar(Ar > 1E30) = []; naceAr = [naceAr; Ar];

            % OSS
            N2 = data{16}; ossN2timestamps = [ossN2timestamps; thisFileTimestamps(N2 < 1E30)]; N2(N2 > 1E30) = []; ossN2 = [ossN2; N2];
            O2 = data{17}; ossO2timestamps = [ossO2timestamps; thisFileTimestamps(O2 < 1E30)]; O2(O2 > 1E30) = []; ossO2 = [ossO2; O2];
            He = data{18}; ossHEtimestamps = [ossHEtimestamps; thisFileTimestamps(He < 1E30)]; He(He > 1E30) = []; ossHe = [ossHe; He];
            O = data{19}; ossOtimestamps = [ossOtimestamps; thisFileTimestamps(O < 1E30)]; O(O > 1E30) = []; ossO = [ossO; O];
            Ar = data{20}; ossARtimestamps = [ossARtimestamps; thisFileTimestamps(Ar < 1E30)]; Ar(Ar > 1E30) = []; ossAr = [ossAr; Ar];
            N = data{21}; ossNtimestamps = [ossNtimestamps; thisFileTimestamps(N < 1E30)]; N(N > 1E30) = []; ossN = [ossN; N];
        else
            % NACE
            N2 = data{10}; naceN2timestamps = [naceN2timestamps; thisFileTimestamps(N2 < 1E30)]; N2(N2 > 1E30) = []; naceN2 = [naceN2; N2];
            O = data{11}; naceOtimestamps = [naceOtimestamps; thisFileTimestamps(O < 1E30)]; O(O > 1E30) = []; naceO = [naceO; O];
            He = data{12}; naceHEtimestamps = [naceHEtimestamps; thisFileTimestamps(He < 1E30)]; He(He > 1E30) = []; naceHe = [naceHe; He];
            Ar = data{13}; naceARtimestamps = [naceARtimestamps; thisFileTimestamps(Ar < 1E30)]; Ar(Ar > 1E30) = []; naceAr = [naceAr; Ar];

            % OSS
            N2 = data{15}; ossN2timestamps = [ossN2timestamps; thisFileTimestamps(N2 < 1E30)]; N2(N2 > 1E30) = []; ossN2 = [ossN2; N2];
            O2 = data{16}; ossO2timestamps = [ossO2timestamps; thisFileTimestamps(O2 < 1E30)]; O2(O2 > 1E30) = []; ossO2 = [ossO2; O2];
            He = data{17}; ossHEtimestamps = [ossHEtimestamps; thisFileTimestamps(He < 1E30)]; He(He > 1E30) = []; ossHe = [ossHe; He];
            O = data{18}; ossOtimestamps = [ossOtimestamps; thisFileTimestamps(O < 1E30)]; O(O > 1E30) = []; ossO = [ossO; O];
            Ar = data{19}; ossARtimestamps = [ossARtimestamps; thisFileTimestamps(Ar < 1E30)]; Ar(Ar > 1E30) = []; ossAr = [ossAr; Ar];
            N = data{20}; ossNtimestamps = [ossNtimestamps; thisFileTimestamps(N < 1E30)]; N(N > 1E30) = []; ossN = [ossN; N];

        end
    end
    [aeTimestamps, order] = unique(aeTimestamps);
    aeLatitude = aeLatitude(order);
    aeLongitude = aeLongitude(order);
    aeAltitude = aeAltitude(order);
    aeLst = aeLst(order);
    
    naceData = struct('N2', naceN2, 'O', naceO, 'He', naceHe, 'Ar', naceAr, ...
        'N2Times', naceN2timestamps, 'OTimes', naceOtimestamps, 'HeTimes', naceHEtimestamps, 'ArTimes', naceARtimestamps);
    ossData = struct('N2', ossN2, 'O', ossO, 'O2', ossO2, 'He', ossHe, 'Ar', ossAr, 'N', ossN, ...
        'N2Times', ossN2timestamps, 'OTimes', ossOtimestamps, 'O2Times', ossO2timestamps, 'HeTimes', ossHEtimestamps, 'ArTimes', ossARtimestamps, 'NTimes', ossNtimestamps);
    
    if strcmp(aeSatellite,'c')
        aeCData = struct('timestamps', aeTimestamps, 'nace', naceData, 'oss', ossData, 'latitude', aeLatitude, 'longitude', aeLongitude,...
                  'altitude', aeAltitude, 'solarTime', aeLst);
        aeSatellite = 'e';
    else
        aeEData = struct('timestamps', aeTimestamps, 'nace', naceData, 'oss', ossData, 'latitude', aeLatitude, 'longitude', aeLongitude,...
                  'altitude', aeAltitude, 'solarTime', aeLst);
    end
end

fprintf('%s\n', 'Began reading Temperature files')
aeSatellite = 'c';
for i = 1:3
    if strcmp(aeSatellite,'c')
        aeFiles = dir('aeTemp/ae_c_*');
    elseif strcmp(aeSatellite,'e')
        aeFiles = dir('aeTemp/ae_e_*');
    else
        aeFiles = dir('aeTemp/ae_de2_*');
    end
    
    tempTimes = [];
    allTimes = [];
    lat = [];
    alt = [];
    lon = [];
    temp = [];
    lstAllFiles = [];
    
    for j = 1:length(aeFiles)
        aeFile = fopen(['aeTemp/',aeFiles(j).name]);
        data = textscan(aeFile, '%f %f %f %f %f %f %f %f %f %f','MultipleDelimsAsOne',1);
        thisFileTimestamps = datenum([data{1}+1900, data{2}, data{3}, data{4}, data{5}, data{6}]);
        utHour = data{4} + data{5}/60 + data{6}/3600;
        lst = utHour + data{9} / 15;
        lst(lst >= 24) = lst(lst >= 24) - 24;
        lst(lst < 0) = lst(lst < 0) + 24;
        
        lstAllFiles = [lstAllFiles; lst];
        allTimes = [allTimes; thisFileTimestamps];
        alt = [alt; data{7}];
        lat = [lat; data{8}];
        lon = [lon; data{9}];
        
        T = data{10}; tempTimes = [tempTimes; thisFileTimestamps(T < 1E30)]; T(T > 1E30) = []; temp = [temp; T];
    end
    
    if strcmp(aeSatellite,'c')
        dataStruct = aeCData;
    elseif strcmp(aeSatellite,'e')
        dataStruct = aeEData;
    else
        dataStruct = de2Data;
    end
    
    allTimes = [allTimes; dataStruct.timestamps];
    lat = [lat; dataStruct.latitude];
    lon = [lon; dataStruct.longitude];
    alt = [alt; dataStruct.altitude];
    lstAllFiles = [lstAllFiles; dataStruct.solarTime];
    
    [allTimes, order] = unique(allTimes);
    lat = lat(order);
    alt = alt(order);
    lon = lon(order);
    lstAllFiles = lstAllFiles(order);
    
    tempData = struct('T', temp, 'TTimes', tempTimes);
    
    dataStruct.timestamps = allTimes;
    dataStruct.latitude = lat;
    dataStruct.longitude = lon;
    dataStruct.altitude = alt;
    dataStruct.solarTime = lstAllFiles;
    dataStruct.temp = tempData;
    
    if strcmp(aeSatellite,'c')
        aeCData = dataStruct;
        aeSatellite = 'e';
    elseif strcmp(aeSatellite,'e')
        aeEData = dataStruct;
        aeSatellite = 'de2';
    else
        de2Data = dataStruct;
    end
end

fprintf('%s\n', 'Began reading AEROS files')
aerosTimestamps = [];
aerosAlt = [];
aerosLat = [];
aerosLon = [];
aerosLst = [];
for i = 1:4
    if i == 1
        aerosFile = fopen('aerosa_ar.asc');
    elseif i == 2
        aerosFile = fopen('aerosa_he.asc');
    elseif i == 3
        aerosFile = fopen('aerosa_n2.asc');
    elseif i == 4
        aerosFile = fopen('aerosa_o2.asc');
    end
    data = textscan(aerosFile, '%f %f %f %f %f %f %f %f','MultipleDelimsAsOne',1, 'headerlines', 3);
    dataYear = 1900 + floor(data{1} / 1000);
    dataYearDatenums = datenum(horzcat(dataYear, repmat([1,1,0,0,0], length(dataYear), 1)));
    dataDoy = mod(data{1}, 1000);
    timestamps = dataYearDatenums + dataDoy + data{2}/86400 - 1;
    utHour = data{2} / 3600;
    lst = utHour + data{7} / 15;
    lst(lst >= 24) = lst(lst >= 24) - 24;
    lst(lst < 0) = lst(lst < 0) + 24;
    if i == 1
        ArTimes = timestamps; ArDens = data{3}; ArErr = data{4};
    elseif i == 2
        HeTimes = timestamps; HeDens = data{3}; HeErr = data{4};
    elseif i == 3
        N2Times = timestamps; N2Dens = data{3}; N2Err = data{4};
    elseif i == 4
        O2Times = timestamps; O2Dens = data{3}; O2Err = data{4};
    end
    aerosTimestamps = [aerosTimestamps; timestamps];
    aerosAlt = [aerosAlt; data{5}];
    aerosLat = [aerosLat; data{6}];
    aerosLon = [aerosLon; data{7}];
    aerosLst = [aerosLst; lst];
end
[aerosTimestamps, order] = unique(aerosTimestamps);
aerosAlt = aerosAlt(order);
aerosLat = aerosLat(order);
aerosLon = aerosLon(order);
aerosLst = aerosLst(order);
densities = struct('N2', N2Dens, 'Ar', ArDens, 'He', HeDens, 'O2', O2Dens,...
                   'N2Times', N2Times, 'ArTimes', ArTimes, 'HeTimes', HeTimes, 'O2Times', O2Times,...
                   'N2Err', N2Err, 'ArErr', ArErr, 'HeErr', HeErr, 'O2Err', O2Err);
aerosData = struct('timestamps', aerosTimestamps, 'dens', densities, 'latitude', aerosLat, 'longitude', aerosLon, 'solarTime', aerosLst, 'altitude', aerosAlt);

end

function [tiegcmDensityVariableAlt, tiegcmDensity270km] = readTiegcmFile(timestampsDensityDatenum)

tiegcmDensityVariableAlt = ones(size(timestampsDensityDatenum));
tiegcmDensity270km = ones(size(timestampsDensityDatenum));

if exist('tiegcmDens.mat', 'file')
    load tiegcmDens.mat
    ind1 = ismember(timestampsDensityDatenum, tiegcmGoceDatenums);
    ind2 = ismember(tiegcmGoceDatenums, timestampsDensityDatenum);
    tiegcmDensityVariableAlt(ind1) = tiegcmGoceInterp(ind2) * 1e14;
    tiegcmDensity270km(ind1) = tiegcmGoce270km(ind2) * 1e14;
end

end

function indicesToRemove = roundTimestampsToBeginFromNearestThreeHourMark(timestampsDensityDatenum)
%

[y, m, d, hours, minutes, seconds] = datevec(timestampsDensityDatenum(1));
hours = hours + minutes / 60 + seconds / 3600;
hours = 3 * ceil(hours / 3);
nearest3hTime = datenum(y, m, d, hours, 0, 0);

if isempty(find(timestampsDensityDatenum == nearest3hTime, 1))
    while 1
        
        hours = hours + 3;
        nearest3hTime = datenum(y, m, d, hours, 0, 0);
        
        if ~isempty(find(timestampsDensityDatenum == nearest3hTime, 1))
            break;
        end
    end
end

indicesToRemove = timestampsDensityDatenum < nearest3hTime;

end

function [aeIntGoce, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0] = ...
    computeAeIntegrals(ae, timestampsAe, timestampsGoce, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0)
%
lags = [2 4 8 16 21 30 40 50 60];
aeIntGoce = zeros(length(timestampsGoce), length(lags));
aeIntChamp = zeros(length(champData.timestamps), length(lags));
aeIntGrace = zeros(length(graceData.timestamps), length(lags));
aeIntDe2 = zeros(length(de2Data.timestamps), length(lags));
aeIntAEC = zeros(length(aeCData.timestamps), length(lags));
aeIntAEE = zeros(length(aeEData.timestamps), length(lags));
aeIntAeros = zeros(length(aerosData.timestamps), length(lags));
aeIntSaber = zeros(length(saberData.timestamps), length(lags));
aeIntT0 = zeros(length(T0.timestamps), length(lags));
t = (timestampsAe - timestampsAe(1))*24*60; tInterp = t(1):t(end);
aeInterp = interp1(t, ae, tInterp, 'linear', 0); tInterp = tInterp/1440 + timestampsAe(1);
cumulativeAe = cumsum(aeInterp);

oneHour = 60;
for i = 1:length(lags)
    lag = lags(i) * oneHour;
    aeInt = (cumulativeAe(lag + 1 : end) - cumulativeAe(1 : end - lag)) / (lag);
    aeTime = tInterp(lag + 1 : end);
    aeIntGoce(:,i) = interp1(aeTime, aeInt, timestampsGoce, 'linear', 0);
    aeIntChamp(:,i) = interp1(aeTime, aeInt, champData.timestamps, 'linear', 0);
    aeIntGrace(:,i) = interp1(aeTime, aeInt, graceData.timestamps, 'linear', 0);
    aeIntDe2(:,i) = interp1(aeTime, aeInt, de2Data.timestamps, 'linear', 0);
    aeIntAEC(:,i) = interp1(aeTime, aeInt, aeCData.timestamps, 'linear', 0);
    aeIntAEE(:,i) = interp1(aeTime, aeInt, aeEData.timestamps, 'linear', 0);
    aeIntAeros(:,i) = interp1(aeTime, aeInt, aerosData.timestamps, 'linear', 0);
    aeIntSaber(:,i) = interp1(aeTime, aeInt, saberData.timestamps, 'linear', 0);
    aeIntT0(:,i) = interp1(aeTime, aeInt, T0.timestamps, 'linear', 0);
end

goceData.aeInt = aeIntGoce;
champData.aeInt = aeIntChamp;
graceData.aeInt = aeIntGrace;
de2Data.aeInt = aeIntDe2;
aeCData.aeInt = aeIntAEC;
aeEData.aeInt = aeIntAEE;
aerosData.aeInt = aeIntAeros;
saberData.aeInt = aeIntSaber;
T0.aeInt = aeIntT0;

end

function [F10out, F81Aout, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0] = ...
    giveSolarInputForModels(F10, F81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, indexDatenums, F10datenum, timestampsDensityDatenum, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0)
%
fprintf('%s\n', 'Computing solar index values for models')

msisF81A = interp1(F10datenum, F81A, timestampsDensityDatenum, 'linear', 100);
F30A = interp1(F10datenum, F30A, timestampsDensityDatenum, 'linear', 100);

F10out = interp1(F10datenum + 1, F10, timestampsDensityDatenum, 'linear', 100);
F81Aout = interp1(F10datenum + 1, F81A, timestampsDensityDatenum, 'linear', 100);
S10 = interp1(indexDatenums + 1, S10, timestampsDensityDatenum, 'linear', 100);
S81A = interp1(indexDatenums + 1, S81A, timestampsDensityDatenum, 'linear', 100);
F30 = interp1(F10datenum + 1, F30, timestampsDensityDatenum, 'linear', 100);

M10 = interp1(indexDatenums + 2, M10, timestampsDensityDatenum, 'linear', 100);
M81A = interp1(indexDatenums + 2, M81A, timestampsDensityDatenum, 'linear', 100);

Y10 = interp1(indexDatenums + 5, Y10, timestampsDensityDatenum, 'linear', 100);
Y81A = interp1(indexDatenums + 5, Y81A, timestampsDensityDatenum, 'linear', 100);

goceData.F10 = interp1(F10datenum + 1, F10, timestampsDensityDatenum, 'previous', 100);
goceData.F81A = interp1(F10datenum, F81A, timestampsDensityDatenum, 'previous', 100);

champData.F10 = interp1(F10datenum + 1, F10, champData.timestamps, 'previous', 100);
champData.F81A = interp1(F10datenum, F81A, champData.timestamps, 'previous', 100);

graceData.F10 = interp1(F10datenum + 1, F10, graceData.timestamps, 'previous', 100);
graceData.F81A = interp1(F10datenum, F81A, graceData.timestamps, 'previous', 100);

de2Data.F10 = interp1(F10datenum + 1, F10, de2Data.timestamps, 'previous', 100);
de2Data.F81A = interp1(F10datenum, F81A, de2Data.timestamps, 'previous', 100);

aeCData.F10 = interp1(F10datenum + 1, F10, aeCData.timestamps, 'previous', 100);
aeCData.F81A = interp1(F10datenum, F81A, aeCData.timestamps, 'previous', 100);

aeEData.F10 = interp1(F10datenum + 1, F10, aeEData.timestamps, 'previous', 100);
aeEData.F81A = interp1(F10datenum, F81A, aeEData.timestamps, 'previous', 100);

aerosData.F10 = interp1(F10datenum + 1, F10, aerosData.timestamps, 'previous', 100);
aerosData.F81A = interp1(F10datenum, F81A, aerosData.timestamps, 'previous', 100);

saberData.F10 = interp1(F10datenum + 1, F10, saberData.timestamps, 'previous', 100);
saberData.F81A = interp1(F10datenum, F81A, saberData.timestamps, 'previous', 100);

T0.F10 = interp1(F10datenum + 1, F10, T0.timestamps, 'previous', 100);
T0.F81A = interp1(F10datenum, F81A, T0.timestamps, 'previous', 100);

end

function [apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h, ApDaily, am3h, amAver24h, apAllFixed, timestamps3h, timestamps3hFixed, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0] = ...
    giveApValsForMSIS(apGoce, amGoce, ap, timestampsAp, timestamps10sFixed, timestamps1minFixed, timestamps1min, timestampsDensityDatenum, goceData, champData, graceData, de2Data, aeCData, aeEData, aerosData, saberData, T0)
% 
fprintf('%s\n', 'Computing ap values for msis')

threeHinSec = 3 * 60 * 60;
secondsInDay = 24 * 60 * 60;
firstDatenum = datenum(datestr(timestampsDensityDatenum(1), 'yyyy-mm-dd'));
lastDatenum = datenum(datestr(timestampsDensityDatenum(end), 'yyyy-mm-dd'));
timestamps3h = -3 * 8 * threeHinSec + threeHinSec : threeHinSec : ceil(lastDatenum - firstDatenum + 1) * secondsInDay;
firstHour = round(str2double(datestr(timestampsDensityDatenum(1), 'HH')));
timestamps3h(1: firstHour / 3) = [];
timestamps10s = min(timestamps10sFixed) : 10 : max(timestamps10sFixed);
nextThreeHStamp = threeHinSec * ceil(timestamps10s / threeHinSec);

ap10s = reshape(repmat(apGoce', 180 * 6, 1),[],1);
apTime10s = reshape(repmat(timestamps3h, 180 * 6, 1),[],1);
ap10s(1) = [];
apTime10s(1) = [];
apNow = ap10s(ismember(apTime10s, nextThreeHStamp));
ap3h = ap10s(ismember(apTime10s, nextThreeHStamp - threeHinSec));
ap6h = ap10s(ismember(apTime10s, nextThreeHStamp - 2 * threeHinSec));
ap9h = ap10s(ismember(apTime10s, nextThreeHStamp - 3 * threeHinSec));
ap10sSmoothed = reshape(repmat(smooth(apGoce, 9)', 180 * 6, 1),[],1);
ap10sSmoothed(1) = [];
apAver12To33h = ap10sSmoothed(ismember(apTime10s, nextThreeHStamp - 8 * threeHinSec));
apAver36To57h = ap10sSmoothed(ismember(apTime10s, nextThreeHStamp - 16 * threeHinSec));
ApDaily = ap10sSmoothed(ismember(apTime10s, nextThreeHStamp - 4 * threeHinSec));

finalIndices = ismember(timestamps10s, timestamps10sFixed);
apNow = apNow(finalIndices);
ap3h = ap3h(finalIndices);
ap6h = ap6h(finalIndices);
ap9h = ap9h(finalIndices);
apAver12To33h = apAver12To33h(finalIndices);
apAver36To57h = apAver36To57h(finalIndices);
ApDaily = ApDaily(finalIndices);

am3h = interp1(timestamps3h + threeHinSec, amGoce, timestamps10sFixed, 'linear', 'extrap');

averagingMatrix = nan(8, length(amGoce) - 7);

for i = 1:8
    lastIndex = length(amGoce) - 8 + i;
    averagingMatrix(i,:) = amGoce(i:lastIndex)';
end
averagingMatrix = [repmat(averagingMatrix(:,1),1,7) averagingMatrix];
amAver24h = trapz(averagingMatrix)' / 8;

amAver24h = interp1(timestamps3h, amAver24h, timestamps10sFixed, 'linear', 'extrap');

% apNow = interp1(timestamps3h, apAll, timestamps10sFixed, 'linear', 'extrap');
% ap3h = interp1(timestamps3h + threeHinSec, apAll, timestamps10sFixed, 'linear', 'extrap');
% ap6h = interp1(timestamps3h + 2 * threeHinSec, apAll, timestamps10sFixed, 'linear', 'extrap');
% ap9h = interp1(timestamps3h + 3 * threeHinSec, apAll, timestamps10sFixed, 'linear', 'extrap');
% 
% ap24hSmooth = smooth(apAll,8);
% ApDaily = interp1(timestamps3h, ap24hSmooth, timestamps10sFixed, 'linear', 'extrap');
% apAver12To33h = interp1(timestamps3h + 7.5 * threeHinSec, ap24hSmooth, timestamps10sFixed, 'linear', 'extrap');
% apAver36To57h = interp1(timestamps3h + 15.5 * threeHinSec, ap24hSmooth, timestamps10sFixed, 'linear', 'extrap');

timestamps3h(timestamps3h <= firstHour * threeHinSec) = nan(1);
apAllFixed = apGoce(ismember(timestamps3h, timestamps1min));
timestamps3hFixed = timestamps3h(ismember(timestamps3h, timestamps1minFixed))';
timestamps3h = timestamps3h(ismember(timestamps3h, timestamps1min))';

dt = 3 / 24;
ap24hSmooth = smooth(ap, 8);
champData.apNow = interp1(timestampsAp, ap, champData.timestamps, 'previous', 'extrap');
champData.ap3h = interp1(timestampsAp + dt, ap, champData.timestamps, 'previous', 'extrap');
champData.ap6h = interp1(timestampsAp + 2*dt, ap, champData.timestamps, 'previous', 'extrap');
champData.ap9h = interp1(timestampsAp + 3*dt, ap, champData.timestamps, 'previous', 'extrap');
champData.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, champData.timestamps, 'previous', 'extrap');
champData.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, champData.timestamps, 'previous', 'extrap');
champData.ApDaily = interp1(timestampsAp, ap24hSmooth, champData.timestamps, 'previous', 'extrap');

graceData.apNow = interp1(timestampsAp, ap, graceData.timestamps, 'previous', 'extrap');
graceData.ap3h = interp1(timestampsAp + dt, ap, graceData.timestamps, 'previous', 'extrap');
graceData.ap6h = interp1(timestampsAp + 2*dt, ap, graceData.timestamps, 'previous', 'extrap');
graceData.ap9h = interp1(timestampsAp + 3*dt, ap, graceData.timestamps, 'previous', 'extrap');
graceData.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, graceData.timestamps, 'previous', 'extrap');
graceData.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, graceData.timestamps, 'previous', 'extrap');
graceData.ApDaily = interp1(timestampsAp, ap24hSmooth, graceData.timestamps, 'previous', 'extrap');

goceData.apNow = interp1(timestampsAp, ap, goceData.timestamps, 'previous', 'extrap');
goceData.ap3h = interp1(timestampsAp + dt, ap, goceData.timestamps, 'previous', 'extrap');
goceData.ap6h = interp1(timestampsAp + 2*dt, ap, goceData.timestamps, 'previous', 'extrap');
goceData.ap9h = interp1(timestampsAp + 3*dt, ap, goceData.timestamps, 'previous', 'extrap');
goceData.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, goceData.timestamps, 'previous', 'extrap');
goceData.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, goceData.timestamps, 'previous', 'extrap');
goceData.ApDaily = interp1(timestampsAp, ap24hSmooth, goceData.timestamps, 'previous', 'extrap');

de2Data.apNow = interp1(timestampsAp, ap, de2Data.timestamps, 'previous', 'extrap');
de2Data.ap3h = interp1(timestampsAp + dt, ap, de2Data.timestamps, 'previous', 'extrap');
de2Data.ap6h = interp1(timestampsAp + 2*dt, ap, de2Data.timestamps, 'previous', 'extrap');
de2Data.ap9h = interp1(timestampsAp + 3*dt, ap, de2Data.timestamps, 'previous', 'extrap');
de2Data.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, de2Data.timestamps, 'previous', 'extrap');
de2Data.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, de2Data.timestamps, 'previous', 'extrap');
de2Data.ApDaily = interp1(timestampsAp, ap24hSmooth, de2Data.timestamps, 'previous', 'extrap');

aeCData.apNow = interp1(timestampsAp, ap, aeCData.timestamps, 'previous', 'extrap');
aeCData.ap3h = interp1(timestampsAp + dt, ap, aeCData.timestamps, 'previous', 'extrap');
aeCData.ap6h = interp1(timestampsAp + 2*dt, ap, aeCData.timestamps, 'previous', 'extrap');
aeCData.ap9h = interp1(timestampsAp + 3*dt, ap, aeCData.timestamps, 'previous', 'extrap');
aeCData.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, aeCData.timestamps, 'previous', 'extrap');
aeCData.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, aeCData.timestamps, 'previous', 'extrap');
aeCData.ApDaily = interp1(timestampsAp, ap24hSmooth, aeCData.timestamps, 'previous', 'extrap');

aeEData.apNow = interp1(timestampsAp, ap, aeEData.timestamps, 'previous', 'extrap');
aeEData.ap3h = interp1(timestampsAp + dt, ap, aeEData.timestamps, 'previous', 'extrap');
aeEData.ap6h = interp1(timestampsAp + 2*dt, ap, aeEData.timestamps, 'previous', 'extrap');
aeEData.ap9h = interp1(timestampsAp + 3*dt, ap, aeEData.timestamps, 'previous', 'extrap');
aeEData.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, aeEData.timestamps, 'previous', 'extrap');
aeEData.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, aeEData.timestamps, 'previous', 'extrap');
aeEData.ApDaily = interp1(timestampsAp, ap24hSmooth, aeEData.timestamps, 'previous', 'extrap');

aerosData.apNow = interp1(timestampsAp, ap, aerosData.timestamps, 'previous', 'extrap');
aerosData.ap3h = interp1(timestampsAp + dt, ap, aerosData.timestamps, 'previous', 'extrap');
aerosData.ap6h = interp1(timestampsAp + 2*dt, ap, aerosData.timestamps, 'previous', 'extrap');
aerosData.ap9h = interp1(timestampsAp + 3*dt, ap, aerosData.timestamps, 'previous', 'extrap');
aerosData.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, aerosData.timestamps, 'previous', 'extrap');
aerosData.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, aerosData.timestamps, 'previous', 'extrap');
aerosData.ApDaily = interp1(timestampsAp, ap24hSmooth, aerosData.timestamps, 'previous', 'extrap');

saberData.apNow = interp1(timestampsAp, ap, saberData.timestamps, 'previous', 'extrap');
saberData.ap3h = interp1(timestampsAp + dt, ap, saberData.timestamps, 'previous', 'extrap');
saberData.ap6h = interp1(timestampsAp + 2*dt, ap, saberData.timestamps, 'previous', 'extrap');
saberData.ap9h = interp1(timestampsAp + 3*dt, ap, saberData.timestamps, 'previous', 'extrap');
saberData.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, saberData.timestamps, 'previous', 'extrap');
saberData.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, saberData.timestamps, 'previous', 'extrap');
saberData.ApDaily = interp1(timestampsAp, ap24hSmooth, saberData.timestamps, 'previous', 'extrap');

T0.apNow = interp1(timestampsAp, ap, T0.timestamps, 'previous', 'extrap');
T0.ap3h = interp1(timestampsAp + dt, ap, T0.timestamps, 'previous', 'extrap');
T0.ap6h = interp1(timestampsAp + 2*dt, ap, T0.timestamps, 'previous', 'extrap');
T0.ap9h = interp1(timestampsAp + 3*dt, ap, T0.timestamps, 'previous', 'extrap');
T0.apAver12To33h = interp1(timestampsAp + 7.5*dt, ap24hSmooth, T0.timestamps, 'previous', 'extrap');
T0.apAver36To57h = interp1(timestampsAp + 15.5*dt, ap24hSmooth, T0.timestamps, 'previous', 'extrap');
T0.ApDaily = interp1(timestampsAp, ap24hSmooth, T0.timestamps, 'previous', 'extrap');

end

function [correctedDensity, msisDensityVariableAlt, msisDensity270km, msisDensity270kmNoAp, jb2008DensityVariableAlt, jb2008Density270km, jb2008Density270kmNoDtc, dtm2013Density270km, dtm2013DensityVariableAlt, dtm2013Density270kmNoAm, hwmU, hwmV, densityIndex, densityIndex1min, densityIndexNoBg, averagedDensity, averagedDensityNoBg, density3h]...
    = relateMsisToDensity(density, altitude, datenumToJulian, timestampsDatenum, doy, timestampsAeDatenum, timestamps1minFixed, timestamps3hFixed, solarTime, latitude, longitude, F10, F81A, msisF81A, S10, S81A, M10, M81A, Y10, Y81A, F30, F30A, am3h, amAver24h, dtc, ApDaily, apNow, ap3h, ap6h, ap9h, apAver12To33h, apAver36To57h)
      
fprintf('%s\n', 'Computing normalized densities with models. This may take even half an hour.')

msisDensity270km = nan(size(density));
msisDensityVariableAlt = nan(size(density));
msisDensity270kmNoAp = nan(size(density));

modelingIndices = 1:length(density);

altitudeInKm = altitude / 1000;

doyDecimal = doy;
doy = ceil(doy);
doy(doy == 0) = 1;

% second of the day
secondsInDay = 24 * 60 * 60;
seconds = (timestampsDatenum - floor(timestampsDatenum)) * secondsInDay;

julianDay = timestampsDatenum + datenumToJulian;

obsYear = datevec(timestampsDatenum);
obsYear = obsYear(:,1);

targetCount = round(length(modelingIndices) / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running MSIS, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

parfor i = modelingIndices
    [~,~,~,~,~,msisDensityVariableAlt(i),~,~,~,~,~]...
      =nrlmsise_mex(doy(i),seconds(i),altitudeInKm(i),latitude(i),longitude(i),solarTime(i),msisF81A(i),F10(i),...
      ApDaily(i),apNow(i),ap3h(i),ap6h(i),ap9h(i),apAver12To33h(i),apAver36To57h(i));

    [~,~,~,~,~,msisDensity270km(i),~,~,~,~,~]...
      =nrlmsise_mex(doy(i),seconds(i),270,latitude(i),longitude(i),solarTime(i),msisF81A(i),F10(i),...
      ApDaily(i),apNow(i),ap3h(i),ap6h(i),ap9h(i),apAver12To33h(i),apAver36To57h(i));

    [~,~,~,~,~,msisDensity270kmNoAp(i),~,~,~,~,~]...
      =nrlmsise_mex(doy(i),seconds(i),270,latitude(i),longitude(i),solarTime(i),msisF81A(i),F10(i), 3);

    if mod(i, 10000) == 0
     p.progress;
    end
end
p.stop;

jb2008Density270km = nan(size(density));
jb2008DensityVariableAlt = nan(size(density));
jb2008Density270kmNoDtc = nan(size(density));

% p = TimedProgressBar( targetCount, barWidth, ...
%                     'Running JB2008, ETA ', ...
%                     '. Now at ', ...
%                     'Completed in ' );
%                 
% parfor i = modelingIndices
%     [~,~,jb2008DensityVariableAlt(i)] = jb2008_mex(julianDay(i), altitudeInKm(i), latitude(i), longitude(i), F10(i), F81A(i), S10(i),...
%         S81A(i), M10(i), M81A(i), Y10(i), Y81A(i), dtc(i));
%     
%     [~,~,jb2008Density270km(i)] = jb2008_mex(julianDay(i), 270, latitude(i), longitude(i), F10(i), F81A(i), S10(i),...
%         S81A(i), M10(i), M81A(i), Y10(i), Y81A(i), dtc(i));
%     
%     [~,~,jb2008Density270kmNoDtc(i)] = jb2008_mex(julianDay(i), 270, latitude(i), longitude(i), F10(i), F81A(i), S10(i),...
%         S81A(i), M10(i), M81A(i), Y10(i), Y81A(i), 0);
%     
%     if mod(i, 10000) == 0
%      p.progress;
%     end
% end
% p.stop;

dtm2013Density270km = nan(size(density));
dtm2013DensityVariableAlt = nan(size(density));
dtm2013Density270kmNoAm = nan(size(density));
dtm2013_mex();

p = TimedProgressBar( targetCount, barWidth, ...
                    'Running DTM-2013, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
                
for i = modelingIndices
    [~, ~, dtm2013DensityVariableAlt(i),~,~,~] = dtm2013_mex(doyDecimal(i), altitudeInKm(i), latitude(i), longitude(i), ...
        solarTime(i), F30(i), F30A(i), am3h(i), amAver24h(i));
    
    [~, ~, dtm2013Density270km(i),~,~,~] = dtm2013_mex(doyDecimal(i), 270, latitude(i), longitude(i), ...
        solarTime(i), F30(i), F30A(i), am3h(i), amAver24h(i));
    
    [~, ~, dtm2013Density270kmNoAm(i),~,~,~] = dtm2013_mex(doyDecimal(i), 270, latitude(i), longitude(i), ...
        solarTime(i), F30(i), F30A(i), 0.0, 0.0);
    
    if mod(i, 10000) == 0
     p.progress;
    end
end
p.stop;

hwmU = ones(size(density));
hwmV = ones(size(density));

% p = TimedProgressBar( targetCount, barWidth, ...
%                     'Running HWM07, ETA ', ...
%                     '. Now at ', ...
%                     'Completed in ' );
%                 
% parfor i = modelingIndices
%     [hwmU(i), hwmV(i)] = hwm07_mex(obsYear(i), doyDecimal(i), altitudeInKm(i), latitude(i), longitude(i), apNow(i));
%     
%     if mod(i, 10000) == 0
%      p.progress;
%     end
% end
% p.stop;

dtm2013DensityVariableAlt = dtm2013DensityVariableAlt * power(10, 14);
dtm2013Density270km = dtm2013Density270km * power(10, 14);
dtm2013Density270kmNoAm = dtm2013Density270kmNoAm * power(10, 14);

jb2008DensityVariableAlt = jb2008DensityVariableAlt * power(10, 11);
jb2008Density270km = jb2008Density270km * power(10, 11);
jb2008Density270kmNoDtc = jb2008Density270kmNoDtc * power(10, 11);

correctedDensity = density .* msisDensity270km ./ msisDensityVariableAlt;

msisDensityVariableAlt = msisDensityVariableAlt * power(10, 14);
msisDensity270km = msisDensity270km * power(10, 14);
msisDensity270kmNoAp = msisDensity270kmNoAp * power(10, 14);

correctionFactor = mean(msisDensityVariableAlt ./ density);
densityIndex = (correctionFactor .* correctedDensity) - msisDensity270kmNoAp;
densityIndex1min = smooth(densityIndex, 7);
densityIndex1min = densityIndex1min(ismember(timestampsDatenum, timestampsAeDatenum));
densityIndexNoBg = removePeriodicBackground(densityIndex1min, 30, 1, 0);
densityIndexNoBg = normalize(densityIndexNoBg, densityIndex1min);

averagedDensity = smooth(correctedDensity, 7);
averagedDensity = averagedDensity(ismember(timestampsDatenum, timestampsAeDatenum));
averagedDensityNoBg = removePeriodicBackground(averagedDensity, 125, 1, 0);
averagedDensityNoBg = normalize(averagedDensityNoBg, averagedDensity);

density3h = smooth(averagedDensityNoBg, 179);
density3h = density3h(find(ismember(timestamps1minFixed, timestamps3hFixed)) - 90);

end

function [crwindEast, crwindNorth] = computeCrossTrackWind(fullU, fullV, goceCrwindEast, goceCrwindNorth)

goceMag = sqrt(goceCrwindEast.^2 + goceCrwindNorth.^2);
crossTrackX = goceCrwindEast ./ goceMag;
crossTrackY = goceCrwindNorth ./ goceMag;

projectionMag = fullU .* crossTrackX + fullV .* crossTrackY;

crwindEast = projectionMag .* crossTrackX;
crwindNorth = projectionMag .* crossTrackY;

end

function [morningTimestamps10s, eveningTimestamps10s, morningLatitude, eveningLatitude, morningDoy, eveningDoy] = ...
    splitBySolarTime(timestamps10s, latitude, doy, solarTime)
%

morningIndices = find(solarTime <= 12);
eveningIndices = find(solarTime > 12);

morningTimestamps10s = timestamps10s(morningIndices);
morningLatitude = latitude(morningIndices);
morningDoy = doy(morningIndices);

eveningTimestamps10s = timestamps10s(eveningIndices);
eveningLatitude = latitude(eveningIndices);
eveningDoy = doy(eveningIndices);

end

function [timeByLatitude, latitudeBins] = computeTimeCells(timestamps10s, firstDatenum, latitude)
%

analysisLimit = datenum('2013-06-16', 'yyyy-mm-dd');

timeInDays = timestamps10s / 86400 + firstDatenum;
indicesToConserve = timeInDays < analysisLimit;
allTimestamps10s = timestamps10s;
timestamps10s = timestamps10s(indicesToConserve);
latitude = latitude(ismember(allTimestamps10s, timestamps10s));

[minLatitude, maxLatitude] = findInterpolationLimits(latitude);


latitudeBins = minLatitude + 1 : 3 : maxLatitude - 1;
if maxLatitude - 1 > max(latitudeBins)
    latitudeBins(end + 1) = max(latitudeBins) + 3;
end

timeByLatitude = cell(length(latitudeBins), 1);
indicesToRemove = false(length(latitudeBins),1);
parfor i = 1:length(latitudeBins)
    indices = (latitudeBins(i) - 1.5 < latitude & latitude <= latitudeBins(i) + 1.5);
    if isempty(find(indices, 1))
        indicesToRemove(i) = 1;
        continue;
    end
    timeByLatitude{i} = timestamps10s(indices);
end

latitudeBins = latitudeBins(~indicesToRemove);
timeByLatitude = timeByLatitude(~indicesToRemove);

end

function [morningFourierGrid, eveningFourierGrid, morningLatDoy, eveningLatDoy] = computeFourierGrids(morningGrid, morningBins, eveningGrid, eveningBins, doy, morningLatitude, eveningLatitude, morningDoy, eveningDoy, timestamps10sFixed, densityIndex, aeIntegral, F107median)
%

fprintf('%s\n', 'Computing fourier fits.')

targetCount = 2*length(morningBins);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running Parametrizations, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
%aeLimit = aeIntegral < quantile(aeIntegral, 0.99);

binWidth = 10;
doy1day = (1:365)';

parfor i = 1:length(morningBins)
    timeThisLat = morningGrid{i};
    thisLatIndices = ismember(timestamps10sFixed, timeThisLat);
    densityThisLat = densityIndex(thisLatIndices);
    
    doyThisLat = ceil(doy(thisLatIndices));
    
    thisLatDens = zeros(length(doy1day),1);
    for k = 1:length(doy1day)
        indices = doyThisLat == doy1day(k);
        thisLatDens(k,:) = mean(densityThisLat(indices));
    end
   
    morningDensities(:,i) = thisLatDens; 
    
    p.progress;
end

parfor i = 1:length(eveningBins)
    timeThisLat = eveningGrid{i};
    thisLatIndices = ismember(timestamps10sFixed, timeThisLat);
    densityThisLat = densityIndex(thisLatIndices);

    doyThisLat = ceil(doy(thisLatIndices));
    
    thisLatDens = zeros(length(doy1day),1);
    for k = 1:length(doy1day)
        indices = doyThisLat == doy1day(k);
        thisLatDens(k,:) = mean(densityThisLat(indices));
    end

    eveningDensities(:,i) = thisLatDens;
    
    p.progress;
end

% latitude = morningBins;
% morningDensities = zeros(length(doy1day), length(latitude));
% eveningDensities = zeros(length(doy1day), length(latitude));
% morningSolarTime = 7.0;
% eveningSolarTime = 19.0;
% seconds = 0 : 3 * 60 * 60 : 24 * 60 * 60 - 1;
% morningMsis270km = zeros(length(seconds), 1);
% morningMsisNoAp = zeros(length(seconds), 1);
% eveningMsis270km = zeros(length(seconds), 1);
% eveningMsisNoAp = zeros(length(seconds), 1);
% 
% F107A = F107median; % = 107
% F107 = F107median;
% ApDaily = 27;
% apNow = 80;
% ap3h = 207;
% ap6h = 48;
% ap9h = 15;
% apAver12To33h = 7;
% apAver36To57h = 7;
% for i = 1:length(latitude)
%     for j = 1:length(doy1day)
%         for k = 1:length(seconds)
%             morningLongitude = 180 * (morningSolarTime - (seconds(k)/3600)) / 12;
%             eveningLongitude = 180 * (eveningSolarTime - (seconds(k)/3600)) / 12;
%             morningLongitude(morningLongitude < -180) = morningLongitude + 360;
%             eveningLongitude(eveningLongitude > 180) = eveningLongitude - 360;
%             
%             [~,~,~,~,~,morningMsis270km(k),~,~,~,~,~]...
%             =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),morningLongitude,morningSolarTime,F107A,F107,...
%             ApDaily,apNow,ap3h,ap6h,ap9h,apAver12To33h,apAver36To57h);
% 
%             [~,~,~,~,~,morningMsisNoAp(k),~,~,~,~,~]...
%             =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),morningLongitude,morningSolarTime,F107A,F107,3);
%         
%             [~,~,~,~,~,eveningMsis270km(k),~,~,~,~,~]...
%             =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),eveningLongitude,eveningSolarTime,F107A,F107,...
%             ApDaily,apNow,ap3h,ap6h,ap9h,apAver12To33h,apAver36To57h);
% 
%             [~,~,~,~,~,eveningMsisNoAp(k),~,~,~,~,~]...
%             =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),eveningLongitude,eveningSolarTime,F107A,F107,3);            
%         end
%         
%         morningResidue = mean(morningMsis270km - morningMsisNoAp);
%         eveningResidue = mean(eveningMsis270km - eveningMsisNoAp);
%         morningDensities(j,i) = morningResidue * 1e14;        
%         eveningDensities(j,i) = eveningResidue * 1e14;
%     end
%     p.progress;    
% end
% 
% p.stop;

morningDensities = morningDensities';
eveningDensities = eveningDensities';
morningDensities = bsxfun(@minus, morningDensities, mean(morningDensities));
eveningDensities = bsxfun(@minus, eveningDensities, mean(eveningDensities));
morningDensities = bsxfun(@rdivide, morningDensities, std(morningDensities));
eveningDensities = bsxfun(@rdivide, eveningDensities, std(eveningDensities));
morningDensities = bsxfun(@minus, morningDensities, mean(morningDensities,2));
eveningDensities = bsxfun(@minus, eveningDensities, mean(eveningDensities,2));
% h = 1/9*ones(3);
% morningDensities = filter2(h,morningDensities);
% eveningDensities = filter2(h,eveningDensities);
% morningDensities = morningDensities - min(morningDensities(:));
% morningDensities = morningDensities / max(morningDensities(:));
% eveningDensities = eveningDensities - min(eveningDensities(:));
% eveningDensities = eveningDensities / max(eveningDensities(:));

x = [doy1day - 365; doy1day; doy1day + 365];
fitCoeffs(10) = 2 * pi / 365;
parfor i = 1:length(morningBins)
    y = repmat(morningDensities(i,:)', 3, 1);
    fourierFit = fit(x, y, 'fourier4', 'StartPoint', fitCoeffs);
    morningFourierGrid{i} = fourierFit;
    morningPlot(i,:) = fourierFit(doy1day');
end

parfor i = 1:length(eveningBins)
    y = repmat(eveningDensities(i,:)', 3, 1);
    fourierFit = fit(x, y, 'fourier4', 'StartPoint', fitCoeffs);
    eveningFourierGrid{i} = fourierFit;
    eveningPlot(i,:) = fourierFit(doy1day');
end

[morningDoyGrid, morningLatGrid] = meshgrid(1:365, morningBins);
[eveningDoyGrid, eveningLatGrid] = meshgrid(1:365, eveningBins);
%figure('units','normalized','outerposition',[0 0 1 1]);
figure;
subplot(2,1,1)
surf(morningDoyGrid, morningLatGrid, morningPlot)
title('Morning Fourier Grid')
view(2);
shading flat
colorbar;
subplot(2,1,2)
surf(eveningDoyGrid, eveningLatGrid, eveningPlot)
title('Evening Fourier Grid')
view(2);
shading flat
colorbar;

figure;
subplot(2,1,1)
surf(morningDoyGrid, morningLatGrid, morningDensities)
title('Morning Residue Grid')
view(2);
shading flat
colorbar;
subplot(2,1,2)
surf(eveningDoyGrid, eveningLatGrid, eveningDensities)
title('Evening Residue Grid')
view(2);
shading flat
colorbar;

morningLatDoy = interp2(morningDoyGrid, morningLatGrid, morningDensities, morningDoy, morningLatitude, 'spline');
eveningLatDoy = interp2(eveningDoyGrid, eveningLatGrid, eveningDensities, eveningDoy, eveningLatitude, 'spline');

end

function [magneticLocalTime] = computeMagneticTime(magneticLongitude, doy, timestampsDatenum)
%

hours = (timestampsDatenum - floor(timestampsDatenum)) * 24;

subSolarLat = 23.5 * sin(0.0172 * doy - 1.405);
subSolarLon = 15 * (12 - hours);

altitude = ones(size(doy)) * 270e3;
[~, magneticSubSolarLon] = convertToMagneticCoordinates(subSolarLat, subSolarLon, altitude);

magneticLocalTime = 12 + (magneticLongitude - magneticSubSolarLon) / 15;
tooBigTimes = magneticLocalTime >= 24;
tooSmallTimes = magneticLocalTime < 0;
magneticLocalTime(tooBigTimes) = magneticLocalTime(tooBigTimes) - 24;
magneticLocalTime(tooSmallTimes) = magneticLocalTime(tooSmallTimes) + 24;

end
