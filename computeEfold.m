function [] = computeEfold()

load('ilData.mat')
model = 'il';

computeBestInt(HeStruct,model);

end

function bestInt = computeBestInt(S,model)

[sb, se, comb] = findStormsForSat(S, 'ae', 400, 0, 2, true);
load coeffsAll
if strcmpi(model,'msis')
    [Tex] = computeMsis(S);
    [T0, dT0] = computeMsisDtmLb(S);
else
    T0 = evalT0(S,T0Coeffs);
    dT0 = evalDT(S,dTCoeffs);
    Tex = evalTex(S,optCoeff(TexInd));
end
S = computeDensityRHS(S, Tex, dT0, T0);

efolds = 1:0.5:24;
tauVec = (1:24)';
corrMat = zeros(length(sb),length(efolds));
corrMat_aver = zeros(length(sb),length(efolds));
rmInd = [];
numDataPoints = [];
numDataPoints_aver = [];
minLat = 0; maxLat = 20;
magLatAbs = abs(convertToMagneticCoordinates(S.latitude, S.longitude,...
                                                S.altitude));
for j = 1:length(sb)
    ind = sb(j):se(j);
    ind(magLatAbs(ind) < minLat | magLatAbs(ind) > maxLat) = [];
    numDataPoints = [numDataPoints, length(ind)];
    rho = S.rhs(ind);
    rho_aver = computeOrbitAverage(rho, S.latitude(ind),S.timestamps(ind));
    if length(rho_aver) < 4
        rmInd = [rmInd,j];
    else
        numDataPoints_aver = [numDataPoints_aver,length(rho_aver)];
    end
    for i = 1:length(efolds)        
        aeInt = interp1(tauVec,S.aeInt(ind,:)',clamp(1.0,efolds(i),24.0))';
        ae_aver = computeOrbitAverage(aeInt, S.latitude(ind),S.timestamps(ind));
        c = corrcoef(rho,aeInt);
        corrMat(j,i) = abs(c(2,1));
        if length(rho_aver) >= 4
            c = corrcoef(rho_aver,ae_aver);
            corrMat_aver(j,i) = abs(c(2,1));            
        end
    end
end
%corrMat(rmInd,:) = [];
corrMat_aver(rmInd,:) = [];

[x,y] = meshgrid(efolds,1:size(corrMat,1));
figure;
surf(x,y,corrMat,'edgecolor','none');
view(2);
colorbar;
axis tight;

[~,maxInd]= max(corrMat,[],2);
bestLag_raw = efolds(maxInd);
i = bestLag_raw ~= 1 & bestLag_raw ~= 24;
bestInt_raw = sum(numDataPoints(i).*bestLag_raw(i)) / sum(numDataPoints(i))

[x,y] = meshgrid(efolds,1:size(corrMat_aver,1));
figure;
surf(x,y,corrMat_aver,'edgecolor','none');
view(2);
colorbar;
axis tight;

[~,maxInd]= max(corrMat_aver,[],2);
bestLag_aver = efolds(maxInd);
i = bestLag_aver ~= 1 & bestLag_aver ~= 24;
bestInt_aver = sum(numDataPoints_aver(i).*bestLag_aver(i)) / sum(numDataPoints_aver(i))

bestInt = 0;

end

function S = computeDensityRHS(S, Tex, dT0, T0)
%global modelLbHeight

u2kg = 1.660538921E-27;
k = 1.38064852E-23;
g = 9.418; 
if strcmpi(S.name, 'O')
    molecMass = 16 * u2kg;
    alpha = 0;
elseif strcmpi(S.name, 'N2')
    molecMass = 28 * u2kg;
    alpha = 0;
elseif strcmpi(S.name, 'O2')
    molecMass = 32 * u2kg;
    alpha = 0;
elseif strcmpi(S.name, 'Ar')
    molecMass = 40 * u2kg;
    alpha = 0;
elseif strcmpi(S.name, 'He')
    molecMass = 4 * u2kg;
    alpha = -0.38;
else
    error('Incorrect name for gas species!')
end

sigma = dT0 ./ (Tex - T0);
T = Tex - (Tex - T0) .* exp(-sigma .* (S.Z));
gamma = molecMass * g ./ (k * sigma * 1E-3 .* Tex);
altTerm = (1 + gamma + alpha) .* log(T0 ./ T) - gamma .* sigma .* (S.Z);
S.rhs = log(S.data) - altTerm;

end