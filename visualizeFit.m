function [] = visualizeFit()

if exist('optCoeff.mat', 'file')
    load optCoeff.mat
else
    error('File optCoeff.mat not found!')
end

% TODO: Ei huomioi numBiases:a!
TexCoeff = optCoeff(TexInd);
OCoeff = optCoeff(OInd);
N2Coeff = optCoeff(N2Ind);
HeCoeff = optCoeff(HeInd);

plotSurfs(270, 100, 100, 1, 'Rho', TexCoeff, OCoeff, N2Coeff, HeCoeff);

end

function [] = plotSurfs(heights, F, FA, doy, paramName, TexCoeff, OCoeff, N2Coeff, HeCoeff)

lat = -90:1:90;
lst = 0:0.1:23.9;
N = length(lat)*length(lst);

S.aeInt = 20 * ones(N, 9);
S.latitude = zeros(N, 1);
S.solarTime = zeros(N,1);
S.F = F * ones(N,1);
S.FA = FA * ones(N,1);
S.doy = doy * ones(N,1);
[lstGrid, latGrid] = meshgrid(lst, lat);
S = computeLatLstGrid(S, lat, lst);
S.numBiases = 0;
S = computeVariablesForFit(S);

if ~strcmpi(paramName, 'Tex')
    [Tex, dT0, T0] = findTemps(S, TexCoeff);
    OlbDens = evalSpecies(S, OCoeff);
    N2lbDens = evalSpecies(S, N2Coeff);
    HelbDens = evalSpecies(S, HeCoeff);
end

for a = 1:length(heights)
    S.Z = ones(N,1) * heights(a);

    if strcmpi(paramName, 'Tex')
        param = evalTex(S, TexCoeff);
        msisParam = computeMsis(S);
    elseif strcmpi(paramName, 'rho')
        param = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens);
        [~,msisParam] = computeMsis(S);
    elseif strcmpi(paramName, 'O')
        [~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens);
        [~,~,msisParam] = computeMsis(S);
    elseif strcmpi(paramName, 'N2')
        [~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens);
        [~,~,~,msisParam] = computeMsis(S);
    elseif strcmpi(paramName, 'He')
        [~,~,~,param] = computeRho(T0, dT0, Tex, S.Z, OlbDens, N2lbDens, HelbDens);   
        [~,~,~,~,msisParam] = computeMsis(S);
    else
        error(['Unknown variable: ', paramName])
    end
    
    param = reshape(param, length(lat), length(lst));
    msisParam = reshape(msisParam, length(lat), length(lst));
    
    figure;
    
    subplot(2,1,1);
    clims = plotSurfSubplot(lstGrid, latGrid, param, paramName, FA, doy, heights(a), 16);
    
    subplot(2,1,2);
    plotSurfSubplot(lstGrid, latGrid, msisParam, ['MSIS ', paramName], FA, doy, heights(a), 16, clims);
    
end

end

function clims = plotSurfSubplot(lstGrid, latGrid, param, paramName, FA, doy, height, fs, clims)

surf(lstGrid, latGrid, param, 'edgecolor', 'none');
view(2);
colorbar;
%colormap jet
axis tight;
xlabel('LST', 'fontsize', fs)
ylabel('lat', 'fontsize', fs)
if (nargin == 9)
    caxis(clims)
end
clims = caxis;
if strcmpi(paramName, 'Tex')
    title(['Tex: F10.7=', num2str(FA), ', doy=', num2str(doy)], 'fontsize', fs)
else
    title([paramName, ': F10.7=', num2str(FA), ', doy=', num2str(doy), ', Z=', num2str(height)], 'fontsize', fs)
end
set(gca, 'fontsize', fs);

end

function [Tex, rho, O, N2, He] = computeMsis(S)

N = length(S.latitude);
Tex = zeros(N,1);
rho = zeros(N,1);
O = zeros(N,1);
N2 = zeros(N,1);
He = zeros(N,1);

R = 6317E3;
S.altitude = zeros(N,1);
S.altitude = (R*S.Z) ./ (R - S.Z*1000);

for i = 1:N
    [He(i), O(i), N2(i),~,~,rho(i),~,~,~,Tex(i)] = nrlmsise_mex(S.doy(i),43200,S.altitude(i),S.latitude(i),180,S.solarTime(i),S.FA(i),S.F(i), 3);
end

rho = rho * 1E3;

end

function S = computeLatLstGrid(S, lat, lst)

for j = 1:length(lst)
    for i = 1:length(lat)
        ind = i + (j-1)*length(lat);
        S.latitude(ind) = lat(i);
        S.solarTime(ind) = lst(j);
    end
end

end

function [Tex, dT0, T0] = findTemps(varStruct, TexCoeff)

origNumBiases = varStruct.numBiases; varStruct.numBiases = 0;
Tex_est = evalTex(varStruct, TexCoeff);
varStruct.numBiases = origNumBiases;

T0 = 507;
dT0 = 12.6;
Tex = clamp(T0+1, Tex_est, 5000);

end