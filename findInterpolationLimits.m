function [minAllowedLatitude, maxAllowedLatitude] = findInterpolationLimits(magneticLatitude)
%
maxAllowedLatitude = 90;
numOfSameNumOfCrossings = 0;
previousNumOfCrossings = -1;
while 1
   subtractedLatitude = magneticLatitude - maxAllowedLatitude;
   numOfCrossings = length(find(subtractedLatitude(1:end-1) .* subtractedLatitude(2:end) < 0));
   if numOfCrossings == previousNumOfCrossings
      numOfSameNumOfCrossings = numOfSameNumOfCrossings + 1; 
   else
      numOfSameNumOfCrossings = 0;
   end
   
   if numOfSameNumOfCrossings > 2
       break;
   end
   
   previousNumOfCrossings = numOfCrossings;
   maxAllowedLatitude = maxAllowedLatitude - 1;
end

minAllowedLatitude = -90;
numOfSameNumOfCrossings = 0;
previousNumOfCrossings = -1;
while 1
   subtractedLatitude = magneticLatitude - minAllowedLatitude;
   numOfCrossings = length(find(subtractedLatitude(1:end-1) .* subtractedLatitude(2:end) < 0));
   if numOfCrossings == previousNumOfCrossings
      numOfSameNumOfCrossings = numOfSameNumOfCrossings + 1; 
   else
      numOfSameNumOfCrossings = 0;
   end
   
   if numOfSameNumOfCrossings > 2
       break;
   end
   
   previousNumOfCrossings = numOfCrossings;
   minAllowedLatitude = minAllowedLatitude + 1;
end

end