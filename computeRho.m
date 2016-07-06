function [rho, OnumDens, N2numDens, HeNumDens, ArNumDens, O2NumDens] = computeRho(T0, dT0, Tex, Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens)
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

gamma_Ar = 40 * u2kg * g ./ (sigma*1E-3 .* k .* Tex);
f_Ar = (T0 ./ T).^(1+gamma_Ar) .* exp(-sigma .* (Z) .* gamma_Ar);
ArNumDens = ArlbDens.*f_Ar; % [1/cm^3]

gamma_O2 = 32 * u2kg * g ./ (sigma*1E-3 .* k .* Tex);
f_O2 = (T0 ./ T).^(1+gamma_O2) .* exp(-sigma .* (Z) .* gamma_O2);
O2NumDens = O2lbDens.*f_O2; % [1/cm^3]

rho = (16*OnumDens + 28*N2numDens + 4*HeNumDens + 40*ArNumDens + 32*O2NumDens) * u2kg * 1E6; % [kg/m^3]

end