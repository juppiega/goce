function [addStruct] = computeVariablesForFit(addStruct, computeMagDiurnal)

x = cosd(90 - addStruct.latitude);
if isempty(addStruct.longitude)
    addStruct.longitude = [];
    addStruct.latitude = [];
    addStruct.altitude = [];
end
[magLat, magLon] = convertToMagneticCoordinates(addStruct.latitude, addStruct.longitude,...
                                                addStruct.altitude);
x_mag = cosd(90 - magLat);
addStruct.lv = addStruct.longitude * pi / 180;

% First degree functions.
P = legendre_fixed(1, x);
addStruct.P10 = f(1,P); 
addStruct.P11 = f(2,P);

% Second degree.
P = legendre_fixed(2, x);
addStruct.P20 = f(1,P);
addStruct.P21 = f(2,P);
addStruct.P22 = f(3,P);

% Third degree.
P = legendre_fixed(3, x);
addStruct.P30 = f(1,P);
addStruct.P31 = f(2,P);
addStruct.P32 = f(3,P);
addStruct.P33 = f(4,P);

% Fourth degree.
P = legendre_fixed(4, x);
addStruct.P40 = f(1,P);
addStruct.P41 = f(2,P);
addStruct.P42 = f(3,P);
addStruct.P43 = f(4,P);
addStruct.P44 = f(5,P);

% Fifth degree.
P = legendre_fixed(5, x);
addStruct.P50 = f(1,P);
addStruct.P51 = f(2,P);
addStruct.P52 = f(3,P);
addStruct.P53 = f(4,P);
addStruct.P54 = f(5,P);

% Sixth degree.
P = legendre_fixed(6, x);
addStruct.P60 = f(1,P);
addStruct.P61 = f(2,P);
addStruct.P62 = f(3,P);
addStruct.P63 = f(4,P);
addStruct.P64 = f(5,P);

% Seventh degree.
P = legendre_fixed(7, x);
addStruct.P70 = f(1,P);
addStruct.P71 = f(2,P);
addStruct.P74 = f(5,P);

% First degree functions.
P = legendre_fixed(1, x_mag);
addStruct.mP10 = f(1,P);
addStruct.mP11 = f(2,P);

% Second degree.
P = legendre_fixed(2, x_mag);
addStruct.mP20 = f(1,P);
addStruct.mP21 = f(2,P);

% Third degree.
P = legendre_fixed(3, x_mag);
addStruct.mP30 = f(1,P);
addStruct.mP31 = f(2,P);
addStruct.mP32 = f(3,P);

% Fourth degree.
P = legendre_fixed(4, x_mag);
addStruct.mP40 = f(1,P);
addStruct.mP41 = f(2,P);

% Fifth degree.
P = legendre_fixed(5, x_mag);
addStruct.mP50 = f(1,P);
addStruct.mP51 = f(2,P);
addStruct.mP52 = f(3,P);

% Sixth degree.
P = legendre_fixed(6, x_mag);
addStruct.mP60 = f(1,P);

% Seventh degree.
P = legendre_fixed(7, x_mag);
addStruct.mP70 = f(1,P);

% Annual parameter.
if ~isfield(addStruct, 'doy') || length(addStruct.doy) ~= length(x_mag) 
    [yr,~,~,~,~,~] = datevec(addStruct.timestamps);
    yearVec = [yr, repmat([1,1,0,0,0], length(yr), 1)];
    if ~isempty(addStruct.timestamps)
        addStruct.doy = addStruct.timestamps - datenum(yearVec) + 1;
    else
        addStruct.doy = [];
    end
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

%if ~isfield(addStruct,'Z') || (isfield(addStruct,'Z') && length(addStruct.altitude) ~= length(addStruct.Z))
    addStruct = computeGeopotentialHeight(addStruct);
%end

end

function output = f(i,P)

if isempty(P)
    output = [];
else
    output = P(i,:)';
end

end

function P = legendre_fixed(i,x)

if isempty(x)
    P = [];
else
    P = legendre(i,x);
end

end
