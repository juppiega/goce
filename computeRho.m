function [rho, OnumDens, N2numDens, HeNumDens] = computeRho(T0, dT0, Tex, Z, OlbDens, N2lbDens, HelbDens)
%global modelLbHeight
sigma = dT0 ./ (Tex - T0);
g = 9.418;
u2kg =  1.660538921E-27;
k = 1.38064852E-23;

T = Tex - (Tex - T0) .* exp(-sigma .* (Z));

gamma_O = 16 * u2kg * g ./ (sigma*1E-3 .* k .* Tex);
f_O = (T0 ./ T).^(1+gamma_O) .* exp(-sigma .* (Z) .* gamma_O);
OnumDens = OlbDens.*f_O; % [1/cm^3]

gamma_N2 = 28 * u2kg * g ./ (sigma*1E-3 .* k .* Tex);
f_N2 = (T0 ./ T).^(1+gamma_N2) .* exp(-sigma .* (Z) .* gamma_N2);
N2numDens = N2lbDens.*f_N2; % [1/cm^3]

gamma_He = 4 * u2kg * g ./ (sigma*1E-3 .* k .* Tex);
f_He = (T0 ./ T).^(1+gamma_He) .* exp(-sigma .* (Z) .* gamma_He);
HeNumDens = HelbDens.*f_He; % [1/cm^3]

rho = (16*OnumDens + 28*N2numDens + 4*HeNumDens) * u2kg * 1E6; % [kg/m^3]

end