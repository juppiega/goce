function [addStruct] = computeVariablesForFit(addStruct, computeMagDiurnal)

x = cosd(90 - addStruct.latitude);
[magLat, magLon] = convertToMagneticCoordinates(addStruct.latitude, addStruct.longitude,...
                                                addStruct.altitude);
x_mag = cosd(90 - magLat);
addStruct.lv = magLon * pi / 180;

% First degree functions.
P = legendre(1, x);
addStruct.P10 = P(1,:)';
addStruct.P11 = P(2,:)';

% Second degree.
P = legendre(2, x);
addStruct.P20 = P(1,:)';
addStruct.P21 = P(2,:)';
addStruct.P22 = P(3,:)';

% Third degree.
P = legendre(3, x);
addStruct.P30 = P(1,:)';
addStruct.P31 = P(2,:)';
addStruct.P32 = P(3,:)';
addStruct.P33 = P(4,:)';

% Fourth degree.
P = legendre(4, x);
addStruct.P40 = P(1,:)';
addStruct.P43 = P(4,:)';
addStruct.P44 = P(5,:)';

% Fifth degree.
P = legendre(5, x);
addStruct.P50 = P(1,:)';
addStruct.P51 = P(2,:)';
addStruct.P52 = P(3,:)';
addStruct.P53 = P(4,:)';
addStruct.P54 = P(5,:)';

% Sixth degree.
P = legendre(6, x);
addStruct.P60 = P(1,:)';
addStruct.P62 = P(3,:)';
addStruct.P63 = P(4,:)';
addStruct.P64 = P(5,:)';

% Seventh degree.
P = legendre(7, x);
addStruct.P70 = P(1,:)';
addStruct.P71 = P(2,:)';
addStruct.P74 = P(5,:)';

% First degree functions.
P = legendre(1, x_mag);
addStruct.mP10 = P(1,:)';
addStruct.mP11 = P(2,:)';

% Second degree.
P = legendre(2, x_mag);
addStruct.mP20 = P(1,:)';
addStruct.mP21 = P(2,:)';

% Third degree.
P = legendre(3, x_mag);
addStruct.mP30 = P(1,:)';
addStruct.mP31 = P(2,:)';
addStruct.mP32 = P(3,:)';

% Fourth degree.
P = legendre(4, x_mag);
addStruct.mP40 = P(1,:)';
addStruct.mP41 = P(2,:)';

% Fifth degree.
P = legendre(5, x_mag);
addStruct.mP50 = P(1,:)';
addStruct.mP51 = P(2,:)';
addStruct.mP52 = P(3,:)';

% Sixth degree.
P = legendre(6, x_mag);
addStruct.mP60 = P(1,:)';

% Seventh degree.
P = legendre(7, x_mag);
addStruct.mP70 = P(1,:)';

addStruct.F2 = addStruct.F.^2;
addStruct.FA2 = addStruct.FA.^2;
addStruct.FtimesFA = addStruct.F .* addStruct.FA;

% Annual parameter.
if ~isfield(addStruct, 'doy') || length(addStruct.doy) ~= length(x_mag) 
    [yr,~,~,~,~,~] = datevec(addStruct.timestamps);
    yearVec = [yr, repmat([1,1,0,0,0], length(yr), 1)];
    addStruct.doy = addStruct.timestamps - datenum(yearVec) + 1;
end

% Diurnal parameter
addStruct.dv = 2*pi* (addStruct.solarTime) / 24;

if nargin == 1 || computeMagDiurnal
        lstMag = computeMagneticTime(magLon, addStruct.doy, addStruct.timestamps);
    addStruct.dv_mag = 2*pi*lstMag / 24;
else
    addStruct.dv_mag = addStruct.dv;
end

addStruct.yv = 2*pi*(addStruct.doy-1)/365;
addStruct.latitudeTerm = zeros(length(x),1);
addStruct.solarTerm = zeros(length(x),1);
addStruct.annual = zeros(length(x),1);
addStruct.diurnal = zeros(length(x),1);
addStruct.semidiurnal = zeros(length(x),1);
addStruct.terdiurnal = zeros(length(x),1);
addStruct.quaterdiurnal = zeros(length(x),1);
addStruct.geomagnetic = zeros(length(x),1);

if isfield(addStruct,'Z') && length(addStruct.altitude) ~= length(addStruct.Z)
    addStruct = computeGeopotentialHeight(addStruct);
end

end