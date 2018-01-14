function [] = computeStormResponse (date1, date2,lst,alt)

t1 = datenum(date1);
t2 = datenum(date2);
dt = 0.5/24;

load('ilData.mat','rhoStruct');
ind = t1 <= rhoStruct.timestamps & rhoStruct.timestamps < t2;
aeInt = rhoStruct.aeInt(ind,:);
[timestamps, si] = unique(rhoStruct.timestamps(ind));
aeInt = aeInt(si,:);
plotTimes = timestamps(1):dt:timestamps(end);
aeInt = interp1(timestamps, aeInt, plotTimes);

[X,Y] = meshgrid(plotTimes-plotTimes(1), -90:5:90);
S.latitude = Y(:);
S.timestamps = X(:);
S.longitude = 0;
S.solarTime = lst;
S.altitude = alt;
S.F = mean(rhoStruct.F(ind));
S.FA = mean(rhoStruct.FA(ind));
S.aeInt = interp1(plotTimes-plotTimes(1), aeInt, S.timestamps);
S = computeVariablesForFit(S);
S = computeGeopotentialHeight(S);

load optCoeff
T0 = evalT0(S,T0Coeffs);
dT = evalDT(S,dTCoeffs);
Tex = evalTex(S,optCoeff(TexInd)); 

OlbDens = evalMajorSpecies(S, optCoeff(OInd), 5);
N2lbDens = evalMajorSpecies(S, optCoeff(N2Ind), 6);
HelbDens = evalMajorSpecies(S, optCoeff(HeInd), 5);
%ArlbDens = evalMajorSpecies(rhoStruct, coeffStruct.ArCoeff, numBiasesStruct.Ar);
O2lbDens = exp(optCoeff(O2Ind));

ilRho = computeRho(T0, dT, Tex, S.Z, OlbDens, N2lbDens, HelbDens, 0, O2lbDens);

Z = reshape(ilRho, size(X));

figure;
plot(X,Y,Z,'edgecolor','none');
view(2);
colorbar;
axis tight;


end
