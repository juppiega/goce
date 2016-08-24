function [] = readGuvi()

guviFiles = dir('GUVI/t_GUVI*');

Nfiles = length(guviFiles);
lat = zeros(100*Nfiles,1);
lon = zeros(100*Nfiles,1);
sza = zeros(100*Nfiles,1);
lst = zeros(100*Nfiles,1);
timestamps = zeros(100*Nfiles,1);
alt = zeros(100*Nfiles,1);

Ofinal = zeros(100*Nfiles,1);
N2final = zeros(100*Nfiles,1);
O2final = zeros(100*Nfiles,1);
sigO = zeros(100*Nfiles,1);
sigN2 = zeros(100*Nfiles,1);
sigO2 = zeros(100*Nfiles,1);
O_N2final = zeros(100*Nfiles,1);
opt = optimoptions('lsqnonlin', 'Jacobian', 'off', 'Algorithm', 'Levenberg-Marquardt', 'TolFun', 1E-8, ...
                 'TolX', 1E-8, 'Display', 'off', 'FinDiffType', 'central');

targetCount = round(length(guviFiles) / 10);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Reading GUVI Files, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
             
k = 0;
for i = 1:length(guviFiles)
    load(['GUVI/',guviFiles(i).name]);
    if max(IYD)-min(IYD) > 1
        continue;
    end
    ind = k + (1:length(GLAT));
    lat(ind) = GLAT;
    lon(ind) = GLONG;
    sza(ind) = SZA;
    lst(ind) = STLP;
    fname = guviFiles(i).name;
    nameEnd = strfind(fname, '.mat');
    year = str2num(fname(nameEnd-16:nameEnd-13));
    t0 = datenum([year, 1, 1, 0, 0, 0]);
    timestamps(ind) = t0 - 1 + double(IYD) + SEC/86400;
    
    for j = 1:length(GLAT)
        if any(FLAG(:,j))
            continue;
        end
        ind = ZM(:,j) >= 150 & ZM(:,j) <= 300;
        if sum(ind) <= 3
            continue;
        end
        z0 = min(ZM(ind,j))*1000;
        Z = computeGeopotentialHeight(ZM(ind,j), z0);
        Tdata = T(ind,j);
        f = @(X) tempFitFun(Tdata, Z, X);
        T0ind = find(ind,1);
        dT0 = (T(T0ind+1,j) - T(T0ind-1,j)) / (ZM(T0ind+1,j) - ZM(T0ind-1,j));
        x0 = [max(T(:,j)), min(Tdata), dT0];
        [optCoeff,~,~,flag] = lsqnonlin(f,double(x0),[],[],opt);
        if flag <= 0
            continue 
        end
        
        Onorm = computeNormalizedDensity(OX(ind,j),Z,'O',optCoeff(1),optCoeff(3),optCoeff(2));
        O2norm = computeNormalizedDensity(O2(ind,j),Z,'O2',optCoeff(1),optCoeff(3),optCoeff(2));
        N2norm = computeNormalizedDensity(N2(ind,j),Z,'N2',optCoeff(1),optCoeff(3),optCoeff(2));
        
        v = k + j;
        
        if all(~isnan(SIGOX(ind,j)))    
            sigOnorm = computeNormalizedDensity(SIGOX(ind,j),Z,'O',optCoeff(1),optCoeff(3),optCoeff(2));
            sigO2norm = computeNormalizedDensity(SIGO2(ind,j),Z,'O2',optCoeff(1),optCoeff(3),optCoeff(2));
            sigN2norm = computeNormalizedDensity(SIGN2(ind,j),Z,'N2',optCoeff(1),optCoeff(3),optCoeff(2));

            OinvVar = 1./sigOnorm.^2;
            O2invVar = 1./sigO2norm.^2;
            N2invVar = 1./sigN2norm.^2;
          
            sigO(v) = sqrt(1/sum(OinvVar));
            sigO2(v) = sqrt(1/sum(O2invVar));
            sigN2(v) = sqrt(1/sum(N2invVar));

            Ofinal(v) = sum(OinvVar.*Onorm)/(sum(OinvVar));
            O2final(v) = sum(O2invVar.*O2norm)/(sum(O2invVar));
            N2final(v) = sum(N2invVar.*N2norm)/(sum(N2invVar));
        else
            sigO(v) = 0;
            sigO2(v) = 0;
            sigN2(v) = 0;
            
            Ofinal(v) = mean(Onorm);
            O2final(v) = mean(O2norm);
            N2final(v) = mean(N2norm);
        end
        
        O_N2final(v) = O_N2(j);
        alt(v) = z0;
    end
    if mod(i,10) == 0
        p.progress;
    end
    k = k + length(GLAT);
end
p.stop;

rmInd = alt == 0;
lat(rmInd) = [];
lon(rmInd) = [];
sza(rmInd) = [];
lst(rmInd) = [];
timestamps(rmInd) = [];
alt(rmInd) = [];
Ofinal(rmInd) = [];
N2final(rmInd) = [];
O2final(rmInd) = [];
sigO(rmInd) = [];
sigN2(rmInd) = [];
sigO2(rmInd) = [];
O_N2final(rmInd) = [];

densities = struct('O', Ofinal, 'N2', N2final, 'O2', O2final, 'O_N2', O_N2final, ...
                    'sigO', sigO, 'sigN2', sigN2, 'sigO2', sigO2);

guviData = struct('dens', densities, ...
                 'timestamps', timestamps, ...
                 'latitude', lat, 'longitude', lon,...
                  'altitude', alt, 'solarTime', lst, 'zenithAngle', sza);
              
save('guviData.mat', 'guviData');

end

function resid = tempFitFun(T, Z, coeff)

Tex = coeff(1);
T0 = coeff(2);
dT = coeff(3);
sig = dT/(Tex-T0);

Tpred = Tex - (Tex - T0) .* exp(-sig*Z);

resid = double(T - Tpred);

end

function normalizedData = computeNormalizedDensity(data, Z, name, Tex, dT0, T0)
%global modelLbHeight

u2kg = 1.660538921E-27;
k = 1.38064852E-23;
g = 9.342;  % at 156.4 km
if strcmpi(name, 'O')
    molecMass = 16 * u2kg;
    alpha = 0;
elseif strcmpi(name, 'N2')
    molecMass = 28 * u2kg;
    alpha = 0;
elseif strcmpi(name, 'O2')
    molecMass = 32 * u2kg;
    alpha = 0;
elseif strcmpi(name, 'Ar')
    molecMass = 40 * u2kg;
    alpha = 0;
elseif strcmpi(name, 'He')
    molecMass = 4 * u2kg;
    alpha = -0.38;
else
    error('Incorrect name for gas species!')
end

sigma = dT0 ./ (Tex - T0);
T = Tex - (Tex - T0) .* exp(-sigma .* (Z));
gamma = molecMass * g ./ (k * sigma * 1E-3 .* Tex);
altTerm = (1 + gamma + alpha) .* log(T0 ./ T) - gamma .* sigma .* (Z);
normalizedData = exp(log(data) - altTerm);

end
