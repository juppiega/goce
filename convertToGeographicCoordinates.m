function [lat, lon] = convertToGeographicCoordinates(mag_lat, mag_lon)

magXYZ = geod2ecef(mag_lat, mag_lon, zeros(size(mag_lon)))';
ecefToMagTransform = [0.33907 -0.91964 -0.19826; ...
                      0.93826  0.34594  0      ; ...
                      0.06859  0.18602  0.98015];
ecefXYZ = ecefToMagTransform \ magXYZ;
r = sqrt(ecefXYZ(1,:).^2 + ecefXYZ(2,:).^2 + ecefXYZ(3,:).^2);
lat = pi/2 - acos(ecefXYZ(3,:) ./ r);
lat = lat' * 180 / pi;

lon = 180 * atan2(ecefXYZ(2,:), ecefXYZ(1,:))' / pi;


end