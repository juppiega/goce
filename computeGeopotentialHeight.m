function [X] = computeGeopotentialHeight(X, z0)
% [X or X.Z] = km

R = 6356770;
if nargin == 1
    z0 = 130E3;
end
if isstruct(X)
    X.Z = (R + z0) * (X.altitude - z0 * 1E-3) ./ (R + X.altitude*1000);
else
    X = (R + z0) * (X - z0 * 1E-3) ./ (R + X*1000);
end


end