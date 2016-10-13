function [ballisticOutput] = computeBiRhoAndIntTerms(ensemble, modelFunction, previousTLEs, recentTLEs, dt, maxBdeviation)
% [dt] minutes, [maxBdeviation] = relative (e.g. 0.1*Btrue <= Bi <= 10*Btrue => maxBdeviation = 10;

mu = 3.986004418E14;
omega_E = 7.2921159E-5;

ballisticOutput = struct('Bi', [], 'Bratio', [], 'objectIDs',[],'rhoObs',[],'rhoModel',[],...
    'sig_rho', [],'intProperties',containers.Map('KeyType', 'double', 'ValueType', 'any'));

recentObjects = keys(recentTLEs);
Bi = zeros(length(recentObjects), size(ensemble,2));
Bratio = zeros(length(recentObjects), size(ensemble,2));
objectIDs = zeros(length(recentObjects));
rhoObs = zeros(length(recentObjects));
sig_rho = zeros(length(recentObjects));
rhoModel = zeros(length(recentObjects));

for i = 1:length(recentObjects)
    
    object = recentObjects{i};
    if ~isKey(previousTLEs, object)
        warning(['Object ', object, ' could not be found in the previousTLEs!'])
        continue;
    end
    
    % Mean motions
    n1 = previousTLEs(object).sgp4info.no / 60;
    n2 = recentTLEs(object).sgp4info.no / 60;
    nMean = 0.5*(n1 + n2);
    dn = n2 - n1;
    if dn <= 0
        warning(['Object ', object, ' ignored due to decreased mean motion!'])
        continue;
    end
    
    % Inclination
    incl = previousTLEs(object).sgp4info.inclo; % Radians
    
    % Integral multiplier
    d = (2.0/3.0) * mu^(2.0/3.0) * nMean^(-1.0/3.0) * dn;
    
    % Epochs
    epoch1 = previousTLEs(object).sgp4info.jdsatepoch;
    epoch2 = recentTLEs(object).sgp4info.jdsatepoch;
    M = floor((epoch2-epoch1)*1440/dt) + 2;
    tJulDay = zeros(M,1);
    rVec = zeros(M,3);
    vVec = zeros(M,3);
    
    satrec = previousTLEs(object).sgp4info;
    t = epoch1;
    k = 1;
    while t < epoch2
        [satrec, r, v] = sgp4(satrec, (t - epoch1)*1440);
        rVec(k,:) = [r(1), r(2), r(3)];
        vVec(k,:) = [v(1), v(2), v(3)];
        tJulDay(k) = t;
        t = t + dt/1440;
        k = k + 1;
    end
    [satrec, r, v] = sgp4(satrec, (epoch2 - epoch1)*1440);
    rVec(k,:) = [r(1), r(2), r(3)];
    vVec(k,:) = [v(1), v(2), v(3)];
    tJulDay(k) = epoch2;
    
    ind = tJulDay > 0;
    tJulDay = tJulDay(ind);
    rVec = rVec(ind,:);
    vVec = vVec(ind,:);
    
    [lat, lon, alt] = eciToGeodetic(rVec(:,1), rVec(:,2), rVec(:,3), tJulDay);
    UT = 12 + mod(tJulDay, 1.0);
    lst = lon/15 + UT;
    lst(lst>=24) = lst(lst>=24) - 24;
    lst(lst<0) = lst(lst<0) + 24;
    
    rMag = sqrt(sum(rVec.^2, 2));
    vMag = sqrt(sum(vVec.^2, 2));
    
    windFac = (1 - omega_E*rMag.*cos(incl)./vMag).^2;
    
    tDatenum = zeros(size(tJulDay));
    for j = 1:length(tJulDay)
        [yyyy,mo,dd,hh,mins,ss] = invjday(tJulDay(j));
        tDatenum(j) = datenum([yyyy,mo,dd,hh,mins,ss]);
    end
    
    observationStruct = struct('latitude', lat, 'longitude', lon, 'solarTime', lst,...
                               'altitude', alt, 'timestamps', tDatenum);
    observationStruct = computeVariablesForFit(observationStruct);
    
    densityMatrix = zeros(m, size(ensemble,2));
    
    for j = 1:size(ensemble,2)
        densityMatrix(:,j) = modelFunction(ensemble(:,j), observationStruct);
    end
    
    integralMatrix = bsxfun(@times, ((1E3*vMag).^3).*windFac, densityMatrix);
    
    integrationTimes = (tJulDay - tJulDay(1))*86400;
    integrals = trapz(integrationTimes, integralMatrix);
    
    BiThis = d./integrals;
    BtrueThis = recentTLEs(object).Btrue;
    BratioThis = BiThis./recentTLEs(object).Btrue;
    if ~(BtrueThis/maxBdeviation <= mean(BiThis) && mean(BiThis) <= BtrueThis*maxBdeviation)
        warning(['Object ',object,' disqualified, because Bi/Btrue = ', BratioThis])
    end
    Bi(i,:) = BiThis;
    Bratio(i,:) = BratioThis;
    objectIDs(i) = object;
    intV3Fdt = trapz(integrationTimes, ((1E3*vMag).^3).*windFac);
    rhoModel(i) = mean(integrals) / intV3Fdt;
    rhoObs(i) = d / (recentTLEs(object).Btrue * intV3Fdt);
    sig_rho(i) = (recentTLEs(object).sig_Btrue * d/intV3Fdt) / BtrueThis.^2;
    ballisticOutput.intProperties(object) = struct('V3F',((1E3*vMag).^3).*windFac, 'intTimes',integrationTimes,...
                                                    'observationStruct', observationStruct, 'r', rVec,'v',vVec,'F',windFac,...
                                                    'intV3Fdt', intV3Fdt);
end

ind = objectIDs > 0;
Bi = Bi(ind,:);
Bratio = Bratio(ind,:);
objectIDs = objectIDs(ind);
rhoModel = rhoModel(ind);
rhoObs = rhoObs(ind);
sig_rho = sig_rho(ind);

ballisticOutput.Bi = Bi;
ballisticOutput.Bratio = Bratio;
ballisticOutput.objectIDs = objectIDs;
ballisticOutput.rhoModel = rhoModel;
ballisticOutput.rhoObs = rhoObs;
ballisticOutput.sig_rho = sig_rho;

end