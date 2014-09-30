
% 
% aeIntFixed = aeIntegrals(ismember(timestamps1min, timestamps1minFixed),:);
% aeIntFixed = aeIntFixed ./ mean(aeIntFixed(:));
% 
%fitFormula = 'y ~ x1 + x2 + x3 + x4 + x5 + x8 + x9 + x5:x7 + x5:x8 + x9:x10';

proxyMatrix = [rand(length(timestamps1minFixed),1) aeIntFixed];
x5x7 = proxyMatrix(:,5) .* proxyMatrix(:,7);
x5x8 = proxyMatrix(:,5) .* proxyMatrix(:,8);
x9x10 = proxyMatrix(:,9) .* proxyMatrix(:,10);
proxyMatrix = [proxyMatrix(:,[1 2 3 4 5 8 9]) x5x7 x5x8 x9x10];
linearModel = fitlm(proxyMatrix, densityIndex1min);
finalProxyDensity = feval(linearModel, proxyMatrix);

% format compact
% 
% ftests = vertcat(morningFtest, eveningFtest);
% fnums = mean(ftests);
% [fnums, order] = sort(fnums);
% t = 1:length(fnums);
% t = t(order);
% goodPredictors = t(end-10:end);
% coefnames = linearModel.CoefficientNames;
% coefnames(goodPredictors)