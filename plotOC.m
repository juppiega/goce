load('ilData.mat','rhoStruct')
load('ilComparison.mat')
load('msisDtmComparison.mat')
load('quietIndRhoStruct.mat')

rhoStruct = removeAndFixData(rhoStruct,0);
oc_il = rhoStruct.data(quietInd) ./ ilRho(quietInd);
oc_dtm = rhoStruct.data(quietInd) ./ dtmRho(quietInd);

N = length(oc_il);
plot(1:N,oc_il, 1:N,oc_dtm, 1:N,rhoStruct.aeInt(quietInd,7));
