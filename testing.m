%%

format compact

ftests = vertcat(morningFtest, eveningFtest);
fnums = median(ftests);
[fnums, order] = sort(fnums);
t = 1:length(fnums);
t = t(order);
goodPredictors = t(end-10:end);
coefnames = morningFits{1}.CoefficientNames;
coefnames(goodPredictors)

%%
coefVals = morningFits{1}.Coefficients;
coefIntervals = morningFits{1}.coefCI;
morningFits{1}.NumCoefficients
coefVals = coefVals(:,1);

%%
coefficNames = morningFits{1}.CoefficientNames;
coefVals = morningFits{1}.Coefficients;
coefVals = table2array(coefVals(2,1));