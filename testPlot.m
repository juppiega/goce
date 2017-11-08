function [] = testPlot(plotInd)

load('ilData.mat','rhoStruct');

[sb,se,ci]=findStormsForSat(rhoStruct,'ae',800,0,2,true); 

for i = 1:length(plotInd)
    ind = sb(plotInd(i)):se(plotInd(i));
    [rho,t_aver] = computeOrbitAverage(rhoStruct.data(ind),rhoStruct.latitude(ind),rhoStruct.timestamps(ind));
    figure;plotyy(rhoStruct.timestamps(ind),rhoStruct.aeInt(ind,5),t_aver,rho)
end

end