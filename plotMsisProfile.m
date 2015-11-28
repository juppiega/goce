function [  ] = plotMsisProfile( lat, lon, hour, doy, F10A, F10, Ap, minAlt, maxAlt )

altitude = minAlt:maxAlt; % km
temperature = zeros(size(altitude));

seconds = hour * 3600;
lst = hour + lon / 15;
if lst >= 24; lst = lst - 24; end
if lst < 0; lst = lst + 24; end

for i = 1:length(altitude)
    [~,~,~,~,~,~,~,~,~,~,temperature(i)] = nrlmsise_mex(doy,seconds,altitude(i),lat,lon,lst,F10A,F10,Ap);
end

plot(temperature, altitude, 'linewidth', 2.0, 'color', 'k')
title('Temperature profile')
xlabel('T / K')
ylabel('Altitude / km')

end

