function [normalizedData] = computeNormalizedStructs(name)

if strcmpi(name,'O')
    load('ilData.mat','OStruct');
    S = OStruct;
elseif strcmpi(name,'N2')
    load('ilData.mat','N2Struct');
    S = N2Struct;
elseif strcmpi(name,'He')
    load('ilData.mat','HeStruct');
    S = HeStruct;
end
S = computeVariablesForFit(S);
[T0, dT] = computeMsisDtmLb(S);
Tex = computeMsis(S);

normalizedData = computeNormalizedDensity(S.data, S.Z, name, Tex, dT, T0);
S.data = normalizedData;

save([name,'Normalized.mat'],'normalizedData');

%plot(S.timestamps, log(normalizedData),'.', S.timestamps, S.aeInt(:,7)/400 + 23,'.');

efold = 1:0.5:24;
lat = 5:15:90; dlat = lat(2)-lat(1);
corrs = zeros(length(lat), length(efold));
for i = 1:length(lat)
    Slat = removeDataPoints(S, lat(i)-dlat/2 <= abs(S.latitude) & abs(S.latitude) < lat(i)+dlat/2 ,true,true,true,true);
    corrs(i,:) = computeBestEfold(Slat, efold);
end

[X,Y] = meshgrid(efold, lat);
figure;
surf(X,Y,abs(corrs));
view(2);
colorbar;
title(name)

end

function corrs_mean = computeBestEfold(S, efold)

[sb,se,stormInd] = findStormsForSat(S,'ae',500,0,2,true);

corrs = zeros(length(efold), length(sb));
tau = (1:size(S.aeInt,2))';
for i = 1:length(efold)
    for k = 1:length(sb)
        ind = sb(k):se(k);
        aeInt = interp1(tau, S.aeInt(ind,:)', efold(i));
        corrs(i,k) = corr((S.data(ind)), aeInt','type','spearman');
    end
end

corrs_mean = mean(corrs,2);

end