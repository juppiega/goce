for i = 1:12

doyBegin = (i-1)*30;
dataMean(i) = mean(data(doyBegin <= doy & doy <= doyBegin+30));
dataErr(i) = std(data(doyBegin <= doy & doy <= doyBegin+30));

end

figure;errorbar(1:12,dataMean,dataErr);