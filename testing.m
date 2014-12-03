%%

% Selects predictors using an F-test. Terms with signifigance at the 95%
% level for over 95% of the latitude bins are approved for the model.
load ftestsGoceFourierRobust
format compact

ftests = vertcat(morningFtest, eveningFtest);

% Minimum value of F-tests if 95% of the lat. bins are taken into account
fnums = quantile(ftests,0.1);
% Sort for easy manual lookup using the variable explorer
[fnums, order] = sort(fnums);
predictorIndices = 1:length(fnums);
predictorIndices = predictorIndices(order);

% Critical value for 95 % confidence is 3.84 for F(1,~10^5);
goodPredictors = predictorIndices(fnums > 2.71); 
if exist('morningFits', 'var')
    coefnames = morningFits{1}.CoefficientNames;
else
    coefnames = oneMorningFit.CoefficientNames;
end
% Print out significant predictors
coefnames(goodPredictors)

%%
save('ftestsGoceFourierRobust.mat', 'morningFtest')
save('ftestsGoceFourierRobust.mat', 'eveningFtest', '-append')
oneMorningFit = morningFits{1};
save('ftestsGoceFourierRobust.mat', 'oneMorningFit', '-append')

%%
doy1day = (1:365)';
latitude = morningBins;
morningDensities = zeros(length(doy1day), length(latitude));
eveningDensities = zeros(length(doy1day), length(latitude));
morningSolarTime = 7.0;
eveningSolarTime = 19.0;
seconds = 0 : 3 * 60 * 60 : 24 * 60 * 60 - 1;
morningMsis270km = zeros(length(seconds), 1);
morningMsisNoAp = zeros(length(seconds), 1);
eveningMsis270km = zeros(length(seconds), 1);
eveningMsisNoAp = zeros(length(seconds), 1);

F107A = 107; % = 107
F107 = 107;
ApDaily = 27;
apNow = 48;
ap3h = 100;
ap6h = 207;
ap9h = 207;
apAver12To33h = 7;
apAver36To57h = 7;
for i = 1:length(latitude)
    for j = 1:length(doy1day)
        for k = 1:length(seconds)
            morningLongitude = 180 * (morningSolarTime - (seconds(k)/3600)) / 12;
            eveningLongitude = 180 * (eveningSolarTime - (seconds(k)/3600)) / 12;
            morningLongitude(morningLongitude < -180) = morningLongitude + 360;
            eveningLongitude(eveningLongitude > 180) = eveningLongitude - 360;
            
            [~,~,~,~,~,morningMsis270km(k),~,~,~,~,~]...
            =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),morningLongitude,morningSolarTime,F107A,F107,...
            ApDaily,apNow,ap3h,ap6h,ap9h,apAver12To33h,apAver36To57h);

            [~,~,~,~,~,morningMsisNoAp(k),~,~,~,~,~]...
            =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),morningLongitude,morningSolarTime,F107A,F107,3);
        
            [~,~,~,~,~,eveningMsis270km(k),~,~,~,~,~]...
            =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),eveningLongitude,eveningSolarTime,F107A,F107,...
            ApDaily,apNow,ap3h,ap6h,ap9h,apAver12To33h,apAver36To57h);

            [~,~,~,~,~,eveningMsisNoAp(k),~,~,~,~,~]...
            =nrlmsise_mex(doy1day(j),seconds(k),270,latitude(i),eveningLongitude,eveningSolarTime,F107A,F107,3);            
        end
        
        morningResidue = mean(morningMsis270km - morningMsisNoAp);
        eveningResidue = mean(eveningMsis270km - eveningMsisNoAp);
        morningDensities(j,i) = morningResidue * 1e14;        
        eveningDensities(j,i) = eveningResidue * 1e14;
    end   
end

morningDensities = morningDensities' - min(morningDensities(:));
morningDensities = morningDensities / max(morningDensities(:));
eveningDensities = eveningDensities' - min(eveningDensities(:));
eveningDensities = eveningDensities / max(eveningDensities(:));

[morningDoyGrid, morningLatGrid] = meshgrid(1:365, morningBins);
[eveningDoyGrid, eveningLatGrid] = meshgrid(1:365, eveningBins);

figure;
subplot(2,1,1)
surf(morningDoyGrid, morningLatGrid, morningDensities)
title('Morning Residue Grid')
view(2);
shading flat
colorbar;
subplot(2,1,2)
surf(eveningDoyGrid, eveningLatGrid, eveningDensities)
title('Evening Residue Grid')
view(2);
shading flat
colorbar;