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

save([name,'Normalized.mat'],'normalizedData');

%plot(S.timestamps, log(normalizedData),'.', S.timestamps, S.aeInt(:,7)/400 + 23,'.');

efold = 1:0.5:24;
corrs = computeBestEfold(S, efold);
figure; plot(efold, corrs);

end

function corrs = computeBestEfold(S, efold)

corrs = zeros(size(efold));
tau = (1:size(S.aeInt,2))';
for i = 1:length(efold)
    aeInt = interp1(tau, S.aeInt', efold(i));
    corrs(i) = corr(log(S.data), aeInt');
end

end