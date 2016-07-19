lat = -90:15:90;
lst = 0:1:23;
doy = 1:6:366;

sec = 86400/2;
lon = 0;

alt = 130;
F107 = 150;
F107A = F107;
Ap = 3;

meanT = 0;
meanTgrad = 0;
tempGrad = 0;
divisor = (length(lst)*length(lat)*length(doy));

for i = 1:length(lst)
    for j = 1:length(lat)
        for t = 1:length(doy)
            [~,~,~,~,~,~,~,~,~,~,tempAt] = nrlmsise_mex(doy(t),sec,alt,lat(j),lon,lst(i),F107A,F107,Ap);
            [~,~,~,~,~,~,~,~,~,~,tempBelow] = nrlmsise_mex(doy(t),sec,alt-0.5,lat(j),lon,lst(i),F107A,F107,Ap);
            [~,~,~,~,~,~,~,~,~,~,tempAbove] = nrlmsise_mex(doy(t),sec,alt+0.5,lat(j),lon,lst(i),F107A,F107,Ap);
            meanTgrad = meanTgrad + (tempAbove - tempBelow) / divisor;
            meanT = meanT + tempAt / divisor;
        end
    end
end

disp(meanTgrad)
disp(meanT)

%%

lat = -90:1:90;
lon = -180:1:180;
doy = 365;
U = zeros(length(lat), length(lon));
V = zeros(length(lat), length(lon));

for i = 1:length(lon)
    for j = 1:length(lat)
        [U(j,i), V(j,i)] = hwm07_mex(2015, doy, 0, lat(j), lon(i), 0);
    end
end

[X,Y] = meshgrid(lon, lat);
figure;
windSpeed = sqrt(U.^2 + V.^2);
normU = U ./ windSpeed;
normV = V ./ windSpeed;
W = zeros(size(U));
Z = max(windSpeed(:)) * ones(size(U));

surf(X,Y,windSpeed, 'edgecolor', 'none')
view(2);
axis tight;
colorbar;
hold all;

ind = 1:15:numel(X);
quiver3(X(ind),Y(ind),Z(ind),normU(ind),normV(ind),W(ind),'color', 'k');

%%
clear
geop = @(z) 6371E3*z ./ (6371E3 + z*1000);
zeta = @(z, z0) (z - z0) .* (6371E3 + z0*1000) ./ (6371E3 + z*1000);      %@(z, z0) geop(z) - geop(z0); 
T = @(z, z0, dT, T0, Tex) Tex - (Tex - T0) * exp(-dT * zeta(z, z0) / (Tex - T0));

Tex = 1000;

z = 120:400;
dT1 = 10.0;
T01 = 380;
T1 = T(z, 120, dT1, T01, 1000);
plot(T1, z);
axis tight;
hold all;

dz = 1E-4;
z02 = 200;
dT2 = (T(z02+dz, 120, dT1, T01, 1000) - T(z02-dz, 120, dT1, T01, 1000)) / (2*dz);
T02 = T(z02, 120, dT1, T01, 1000);
T2 = T(z, z02, dT2, T02, 1000);
plot(T2, z,'--');

max(T1-T2)




