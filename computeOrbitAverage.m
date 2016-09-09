function [orbAver, timestampsDatenum] = computeOrbitAverage(x, lat, timestampsDatenum)

eqCrossings = find((lat(1:end-1) .* lat(2:end) < 0) & lat(1:end-1) < 0); % Only ascending
timestamps1min = (timestampsDatenum - timestampsDatenum(1)) * 1440;
eqCrossings = [1; eqCrossings];
orbAver = zeros(length(eqCrossings)-1, 1);

for i = 2:length(eqCrossings)
    ind = eqCrossings(i-1):eqCrossings(i);
    orbAver(i-1) = mean(x(ind));
end

rmInd = [1; find(diff(timestamps1min(eqCrossings)) > 150)];
orbAver(rmInd) = [];
timestampsDatenum = timestampsDatenum(eqCrossings(2:end));
timestampsDatenum(rmInd) = [];

end