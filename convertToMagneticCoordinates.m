function [magneticLatitude, magneticLongitude] = convertToMagneticCoordinates(latitude, longitude, altitude)
% [magneticLatitude] = convertToMagneticCoordinates(latitude, longitude, altitude)

ecefXYZ = geod2ecef(latitude, longitude, altitude)';
ecefToMagTransform = [0.33907 -0.91964 -0.19826; ...
                      0.93826  0.34594  0      ; ...
                      0.06859  0.18602  0.98015];
magXYZ = ecefToMagTransform * ecefXYZ;
r = sqrt(magXYZ(1,:).^2 + magXYZ(2,:).^2 + magXYZ(3,:).^2);
magneticLatitude = pi/2 - acos(magXYZ(3,:) ./ r);
magneticLatitude = magneticLatitude' * 180 / pi;

magneticLongitude = 180 * atan2(magXYZ(2,:), magXYZ(1,:))' / pi;

end
