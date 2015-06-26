function [  ] = initializeTIEGCM( datestring, exampleFile)

% Check inputs
if ~(ischar(datestring) && ischar(exampleFile))
    error('datestring and exampleFile must all be strings')
end

if(matlabpool('size')==0)
    matlabpool(4);
end

% Output file name.
filename = ['tiegcm_init.nc'];
%delete(filename);

% Copy file structure from exampleFile
%fileSchema = ncinfo(exampleFile);
%ncwriteschema(filename, fileSchema);

% Current datenum.
thisDatenum = datenum(datestring);
[~,~,~,hours,minutes,~] = datevec(thisDatenum);
if ~(hours == 0 && minutes == 0)
    error('Not allowed to specify hours, only date, e.g. 2010-12-31')
end

%copyVars(filename, exampleFile)

%initializeEmpirically(thisDatenum, filename, exampleFile)

writeTimeVariables(thisDatenum, filename);

end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function [] = writeTimeVariables(thisDatenum, filename)

year = str2double(datestr(thisDatenum, 'yyyy'));
doy = floor(thisDatenum) - datenum(datestr(thisDatenum, 'yyyy'), 'yyyy') + 1;
[~,~,~,hours,minutes,~] = datevec(thisDatenum);
timestep = 30;

ncwrite(filename, 'year', year);
ncwrite(filename, 'timestep', timestep);
ncwrite(filename, 'ut', 0);
ncwrite(filename, 'time', 0);
modelTime = [doy-1; hours; minutes];
ncwrite(filename, 'mtime', modelTime);
iter = (modelTime(1) * 86400 + modelTime(2) * 3600 + modelTime(3) * 60) / timestep;
ncwrite(filename, 'iter', iter);
ncwrite(filename, 'day', modelTime(1));
ncwrite(filename, 'calendar_advance', 0);

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
computeHwm(thisDatenum, geometHeight, Ap, filename, exampleFile);
computeIri(thisDatenum, geometHeight, filename, exampleFile);

end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function [geopotHeight, geometHeight] = computeMsis(thisDatenum, Ap, F10, F10A, filename, exampleFile)

fprintf('%s\n', 'Computing neutrals with MSIS');

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
altInKm = (0:5:1000)';
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
            molarMass = (He*4 + O1Interp(k)*16 + N2*28 + O2Interp(k)*32 + Ar*40 + H*1 + N*14) / ...
                        (He + O1Interp(k) + N2 + O2Interp(k) + Ar + H + N);
            PressureInterp(k) = 1000.0 * rho * gasConst * TnInterp(k) / (molarMass * 1E-3);
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

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function [] = computeHwm(thisDatenum, geometHeight, Ap, filename, exampleFile)

fprintf('%s\n', 'Computing winds with HWM');

lat = ncread(exampleFile, 'lat');
lon = ncread(exampleFile, 'lon');
lev = ncread(exampleFile, 'lev');

year = str2double(datestr(thisDatenum, 'yyyy'));
doy = floor(thisDatenum) - datenum(datestr(thisDatenum, 'yyyy'), 'yyyy') + 1;

UN = zeros(size(geometHeight));
VN = zeros(size(geometHeight));

for i = 1:length(lon)
    for j = 1:length(lat)
        for k = 1:length(lev)
            [UN(i,j,k), VN(i,j,k)] = hwm07_mex(year, doy, geometHeight(i,j,k), lat(j), lon(i), Ap);
        end
    end
end

UN = UN * 100;
VN = VN * 100;

ncwrite(filename, 'ULBC', UN(:,:,1))
ncwrite(filename, 'ULBC_NM', UN(:,:,1))
ncwrite(filename, 'UN', UN)
ncwrite(filename, 'UN_NM', UN)
ncwrite(filename, 'VLBC', VN(:,:,1))
ncwrite(filename, 'VLBC_NM', VN(:,:,1))
ncwrite(filename, 'VN', VN)
ncwrite(filename, 'VN_NM', VN)

end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function [] = computeIri(thisDatenum, geometHeight, filename, exampleFile)

fprintf('%s\n', 'Computing ionosphere with IRI');

lat = ncread(exampleFile, 'lat');
lon = ncread(exampleFile, 'lon');
lev = ncread(exampleFile, 'lev');

NE = zeros(size(geometHeight));
TE = zeros(size(geometHeight));
TI = zeros(size(geometHeight));

barWidth = 50;
targetCount = length(lev)*length(lat);

p = TimedProgressBar( targetCount, barWidth, ...
                    'Running IRI, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

for j = 1:length(lev)
    parfor i = 1:length(lat)
        passed = 0;
        while passed == 0
            try
                passed = 1;
                iriVals = iri2012(thisDatenum, lat(i), lon, geometHeight(:,i,j));
            catch
                passed = 0;
            end
        end
        NE_temp(:,i) = iriVals(:,1) / 1E6;
        TE_temp(:,i) = iriVals(:,5);
        TI_temp(:,i) = iriVals(:,4);
        
        p.progress;
    end
    NE(:,:,j) = NE_temp;
    TE(:,:,j) = TE_temp;
    TI(:,:,j) = TI_temp;
end
p.stop;

ncwrite(filename, 'NE', NE)
ncwrite(filename, 'TI', TI)
ncwrite(filename, 'TE', TE)

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
             'N4S_NM', 'OMEGA', 'POTEN', 'O2P', 'OP', 'OP_NM'};

for i = 1:length(copyNames)
    ncwrite(filename, copyNames{i}, ncread(exampleFile, copyNames{i}))
end

end