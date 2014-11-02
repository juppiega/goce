%%

% Selects predictors using an F-test. Terms with signifigance at the 95%
% level for over 95% of the latitude bins are approved for the model.

format compact

ftests = vertcat(morningFtest, eveningFtest);

% Minimum value of F-tests if 95% of the lat. bins are taken into account
fnums = quantile(ftests,0.05);
% Sort for easy manual lookup using the variable explorer
[fnums, order] = sort(fnums);
predictorIndices = 1:length(fnums);
predictorIndices = predictorIndices(order);

% Critical value for 95 % confidence is 3.84 for F(1,~10^5);
goodPredictors = predictorIndices(fnums > 3.84); 
coefnames = morningFits{1}.CoefficientNames;
% Print out significant predictors
coefnames(goodPredictors)

%%