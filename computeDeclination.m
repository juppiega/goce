function [declination] = computeDeclination (doy)
% Calculation from:
% U.S. Naval Observatory; U.K. Hydrographic Office, H.M. Nautical Almanac Office (2008),
% The Astronomical Almanac for the Year 2010. U.S. Govt. Printing Office. p. C5.
% Year 2000 is assumed

n = doy - 1.5;
L = mod(280.46 + 0.9856474*n, 360);
g = mod(357.528 + 0.9856003*n, 360);
lambda = L + 1.915*sind(g) + 0.020*sind(2*g);

epsilon = 23.439;
declination = asind(sind(epsilon).*sind(lambda));

end