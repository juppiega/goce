meanDens = zeros(365,1);
parfor i = 1:365
    meanDens(i) = mean(densityIndex(ceil(doy)==i & magneticLatitude < 90 & magneticLatitude > 40));
end