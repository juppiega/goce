function [lat, lon, alt] = eciToGeodetic(x, y, z, julianDate)
% [x,y,z] = km in ECI, [julianDate] = days, [lat,lon] = degrees in WGS-72, 
% [alt] = km above the WGS-72 ellipsoid.
% SOURCE: Kelso, T.S., Orbital Coordinate Systems, parts I-III, satellite
% times, 2 (1-3), 1996.
    
    UT = mod(julianDate+0.5, 1.0);
    julianDiff = (julianDate - 2451545 - UT) / 36525;
    gmstMidnight = 24110.54841 + 8640184.812866 * julianDiff + 0.093104 * julianDiff.^2 ...
                    - 6.2E-6 * julianDiff.^3;
    gmst = mod(2*pi*gmstMidnight/86400 + 7.292115E-5 * 86400 * UT, 2*pi);
    
    lon = atan2(y,x) - gmst;
    R = sqrt(x.^2 + y.^2);
    
    lat = atan2(z, R);
    latTol = 1E-10;
    f = 1/298.26; % WGS-72 flattening
    a = 6378.135; % WGS-72 equatorial radius
    e2 = 2*f - f^2;
    
    while true
        s = sin(lat);
        C = 1./(sqrt(1 - e2*s.^2));
        latBetter = atan2(z + a*C*e2*s, R); 
        if all(abs(latBetter - lat) < latTol); break; end
        lat = latBetter;
    end

    alt = R./cos(lat) - a*C;
    
    lat = rad2deg(lat);
    lon = rad2deg(lon);
end