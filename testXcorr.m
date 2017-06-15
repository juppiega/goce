function [] = testXcorr

rho = zeros(7,1); rho([3,4,5]) = 1;
ae = zeros(7,1); ae([2,3,4]) = 1;

maxlag = 4;
crossCorr = xcorr(rho,ae,maxlag);

figure;
plot((1:length(crossCorr))-maxlag-1, crossCorr)

end