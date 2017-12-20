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
[~,~,T0, dT] = computeMsisDtmLb(S);
Tex = computeDtm(S);

normalizedData = computeNormalizedDensity(S.data, S.Z, name, Tex, dT, T0);

save([name,'Normalized.mat'],'normalizedData');

%plot(S.timestamps, log(normalizedData),'.', S.timestamps, S.aeInt(:,7)/400 + 23,'.');

efold = 1:0.5:24;
Spolar = removeDataPoints(S, abs(S.latitude) < 50,true,true,true,true);
corrs = computeBestEfold(Spolar, efold);
figure; plot(efold, corrs); title([name,' polar'])

Seq = removeDataPoints(S, abs(S.latitude) > 30,true,true,true,true);
corrs = computeBestEfold(Seq, efold);
figure; plot(efold, corrs); title([name,' equatorial'])

corrs = computeBestEfold(S, efold);
figure; plot(efold, corrs); title([name,' all'])

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