% format compact
% 
% ftests = vertcat(morningFtest, eveningFtest);
% fnums = mean(ftests);
% goodPredictors = find(fnums > 1000);
% 
% coefnames(goodPredictors)

aeIntFixed = aeIntegrals(ismember(timestamps1min, timestamps1minFixed),:);
aeIntFixed = aeIntFixed ./ mean(aeIntFixed(:));

fitFormula = 'y ~ x1 + x2 + x3 + x4 + x5 + x9 + x13 + x5:x7 + x5:x8 + x14:x16 + x17:x18 + x3:x3';

proxyMatrix = [rand(length(timestamps1minFixed),1) aeIntFixed];
linearModel = fitlm(proxyMatrix, densityIndex1min, fitFormula);
x5x7 = proxyMatrix(:,5) .* proxyMatrix(:,7);
x5x8 = proxyMatrix(:,5) .* proxyMatrix(:,8);
x14x16 = proxyMatrix(:,14) .* proxyMatrix(:,16);
x17x18 = proxyMatrix(:,17) .* proxyMatrix(:,18);
x3x3 = proxyMatrix(:,3) .^2;
proxyMatrix = [ones(length(timestamps1minFixed),1) proxyMatrix(:,[1 2 3 4 5 9 13]) x5x7 x5x8 x14x16 x17x18 x3x3];
finalProxyDensity = feval(linearModel, proxyMatrix);

coefnames = linearModel.CoefficientNames;