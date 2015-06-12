function [  ] = initializeTIEGCM( datestring, exampleFile)

% Check inputs
if ~(ischar(datestring) && ischar(exampleFile))
    error('datestring and exampleFile must all be strings')
end

% Output file name.
filename = ['tiegcm_init.nc'];
delete(filename);

% Copy file structure from exampleFile
fileSchema = ncinfo(exampleFile);
ncwriteschema(filename, fileSchema);

% Current datenum.
thisDatenum = datenum(datestring);
[~,~,~,hours,minutes,~] = datevec(thisDatenum);
if ~(hours == 0 && minutes == 0)
    error('Not allowed to specify hours, only date, e.g. 2010-12-31')
end
%ncwrite(filename, );

copyVars(filename, exampleFile)

initializeEmpirically(thisDatenum, filename, exampleFile)


end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function [] = initializeEmpirically(thisDatenum, filename, exampleFile)

iriVals = iri2012(thisDatenum, 60, 100, 100);
Kp = iriVals(:,44);
Ap = iriVals(:,43);
ap = iriVals(:,42);
F10A = iriVals(:,41);
F10 = iriVals(:,40);
ncwrite(filename, 'Kp', Kp);
ncwrite(filename, 'f107a', F10A);
ncwrite(filename, 'f107d', F10);

[geopotHeight, geometHeight] = computeMsis(thisDatenum, Ap, F10, F10A, filename, exampleFile);


end

function [geopotHeight, geometHeight] = computeMsis(thisDatenum, Ap, F10, F10A, filename, exampleFile)

lat = ncread(exampleFile, 'lat');
lon = ncread(exampleFile, 'lon');
lev = ncread(exampleFile, 'lev');
ilev = ncread(exampleFile, 'ilev');
grav = ncread(exampleFile, 'grav');
avogadroConst = 6.022E23;
gasConst = 8.314;

pressureLevels = ncread(exampleFile, 'p0') * 100 * exp(-lev);
lbPressure = ncread(exampleFile, 'p0') * 100 * exp(-ilev(1));

geometHeight = zeros(length(lon), length(lat), length(lev));
TN = zeros(length(lon), length(lat), length(lev));
O1 = zeros(length(lon), length(lat), length(lev));
O2 = zeros(length(lon), length(lat), length(lev));
TLBC = zeros(length(lon), length(lat));

doy = floor(thisDatenum) - datenum(datestr(thisDatenum, 'yyyy'), 'yyyy') + 1;
secOfDay = (thisDatenum - floor(thisDatenum)) * 86400;
altInKm = (80:5:1000)';
O1Interp = zeros(length(altInKm),1);
O2Interp = zeros(length(altInKm),1);
TnInterp = zeros(length(altInKm),1);
PressureInterp = zeros(length(altInKm),1);

for i = 1:length(lon)
    lst = secOfDay/(60*60) + 12*lon(i)/180;
    for j = 1:length(lat)
        for k = 1:length(altInKm)
            [He, O1Interp(k), N2, O2Interp(k), Ar, rho, H,N,~,~,TnInterp(k)] = ...
                nrlmsise_mex(doy, secOfDay, altInKm(k), lat(j), lon(i), lst, F10A, F10, Ap);
            molarMass = (He*4 + O1Interp(k)*16 + N2*28 + O2Interp(k)*32 + Ar*40 + H*1 + N*14) / 135;
            PressureInterp(k) = rho * gasConst * TnInterp(k) / (molarMass * 1E-3);
            O1Interp(k) = 16 * O1Interp(k) / (rho * avogadroConst);
            O2Interp(k) = 32 * O2Interp(k) / (rho * avogadroConst);
        end
        geometHeight(i,j,:) = interp1(PressureInterp, altInKm, pressureLevels);
        TN(i,j,:) = interp1(PressureInterp, TnInterp, pressureLevels);
        O1(i,j,:) = interp1(PressureInterp, O1Interp, pressureLevels);
        O2(i,j,:) = interp1(PressureInterp, O2Interp, pressureLevels);
        TLBC(i,j) = interp1(PressureInterp, TnInterp, lbPressure);
    end
end

for i = 1:length(lev)
    geopotHeight(:,:,i) = geomet2geopot(geometHeight(:,:,i)*1E3, grav) * 100;
end

ncwrite(filename, 'TN', TN)
ncwrite(filename, 'TN_NM', TN)
ncwrite(filename, 'TLBC', TLBC)
ncwrite(filename, 'TLBC_NM', TLBC)
ncwrite(filename, 'O1', O1)
ncwrite(filename, 'O1_NM', O1)
ncwrite(filename, 'O2', O2)
ncwrite(filename, 'O2_NM', O2)
ncwrite(filename, 'Z', geopotHeight)

end

function [geopotHeight] = geomet2geopot(geometricHeight, grav)

earthRad = repmat(6371E3, size(geometricHeight, 1), size(geometricHeight, 2));
gravPot = 6.673E-11 * 5.972E24 * ((1./earthRad) - (1./(earthRad + geometricHeight)));
geopotHeight = gravPot / (grav / 100);

end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function [] = copyVars(filename, exampleFile)

copyNames = {'write_date', 'tidi_ncfile', 'swvel', 'swden', 'see_ncfile',...
             'sdtide', 'saber_ncfile', 'p0_model', 'p0', 'ntask_mpi', ...
             'ncep_ncfile', 'mlon', 'mlev', 'mlat', 'mag', 'lon', 'lev', ...
             'lat', 'joulefac', 'ilev', 'hpower', 'h2', 'h1', 'grav' ...
             'ed', 'ec', 'e1', 'e2', 'dtide', 'ctpoten', 'crit2', 'crit1', ...
             'coupled_cmit', 'colfac', 'bzimf', 'byimf', 'bximf', ...
             'alfac', 'alfad', 'al', 'LBC', 'N2D', 'NO', 'NO_NM', 'N4S',...
             'N4S_NM', 'OMEGA', 'POTEN'};

for i = 1:length(copyNames)
    ncwrite(filename, copyNames{i}, ncread(exampleFile, copyNames{i}))
end

end