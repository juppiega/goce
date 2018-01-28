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
elseif strcmpi(name,'Tex')
    load('ilData.mat','TexStruct');
    S = TexStruct;
end

load('coeffsAll.s2.mat')

S = computeVariablesForFit(S);
%[T0, dT] = computeMsisDtmLb(S);
%Tex = computeMsis(S);
if ~strcmpi(name,'Tex')
T0 = evalT0(S,T0Coeffs);
dT = evalDT(S,dTCoeffs);
Tex = evalTex(S,optCoeff(TexInd));
end

if ~strcmpi(name,'Tex')
normalizedData = computeNormalizedDensity(S.data, S.Z, name, Tex, dT, T0);
S.data = normalizedData;
save([name,'Normalized.mat'],'normalizedData');
end



%plot(S.timestamps, log(normalizedData),'.', S.timestamps, S.aeInt(:,7)/400 + 23,'.');

magLat = convertToMagneticCoordinates(S.latitude, S.longitude,...
                                                S.altitude);

efold = 1:0.5:24;
Spolar = removeDataPoints(S, abs(magLat) < 50,true,true,true,true);
corrs = computeBestEfold(Spolar, efold);
%figure; plot(1./efold, corrs.^2); title([name,' polar'])
corrs_polar = corrs';
polar_mean = sum(corrs.^2 .* 1./efold') / sum(corrs.^2)

Seq = removeDataPoints(S, abs(magLat) > 30,true,true,true,true);
corrs = computeBestEfold(Seq, efold);
%figure; plot(1./efold, corrs.^2); title([name,' equatorial'])
corrs_eq = corrs';
eq_mean = sum(corrs.^2 .* 1./efold') / sum(corrs.^2)

b0 = 0.5 * (eq_mean + polar_mean)

%corrs = computeBestEfold(S, efold);
%figure; plot(1./efold, corrs.^2); title([name,' all'])

x_eq = 0; x_polar = 1.1509;
x = [ones(size(efold))*x_eq, ones(size(efold))*x_polar];
w = [corrs_eq.^2, corrs_polar.^2];
y = [1./efold, 1./efold];
p = fitlm(x',y','weights',w');

end

function corrs_mean = computeBestEfold(S, efold)

[sb,se,stormInd] = findStormsForSat(S,'ae',800,0,2,true);

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