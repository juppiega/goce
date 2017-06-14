function [result] = computeAEintegralExp(ae1min, timestampsDatenum, lagHours)

dt = diff(timestampsDatenum);
cutInd = find(dt*1440 > 1.5);
cutInd = [0; cutInd; length(ae1min)];
result = zeros(size(ae1min));

sigma = 5;
t = -lagHours*60*sigma : 1 : 0;
crossCorrFunc = exp(t/(lagHours*60));

for i = 1:length(cutInd)-1
    ind = cutInd(i)+1:cutInd(i+1);
    aeThis = ae1min(ind);
    N = length(aeThis);
    crossCorr = xcorr(aeThis, crossCorrFunc);

    numZeroPad = max(N - length(crossCorrFunc),0);
    conserveInd = numZeroPad + (1:N);
    result(ind) = crossCorr(conserveInd);
end

end