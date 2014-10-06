

format compact

ftests = vertcat(morningFtest, eveningFtest);
fnums = median(ftests);
[fnums, order] = sort(fnums);
t = 1:length(fnums);
t = t(order);
goodPredictors = t(end-10:end);
coefnames = morningFits{1}.CoefficientNames;
coefnames(goodPredictors)