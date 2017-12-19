
function normalizedData = computeNormalizedDensity(data, Z, name, Tex, dT0, T0)
%global modelLbHeight

u2kg = 1.660538921E-27;
k = 1.38064852E-23;
g = 9.342;  % at 156.4 km
if strcmpi(name, 'O')
    molecMass = 16 * u2kg;
    alpha = 0;
elseif strcmpi(name, 'N2')
    molecMass = 28 * u2kg;
    alpha = 0;
elseif strcmpi(name, 'O2')
    molecMass = 32 * u2kg;
    alpha = 0;
elseif strcmpi(name, 'Ar')
    molecMass = 40 * u2kg;
    alpha = 0;
elseif strcmpi(name, 'He')
    molecMass = 4 * u2kg;
    alpha = -0.38;
else
    error('Incorrect name for gas species!')
end

sigma = dT0 ./ (Tex - T0);
T = Tex - (Tex - T0) .* exp(-sigma .* (Z));
gamma = molecMass * g ./ (k * sigma * 1E-3 .* Tex);
altTerm = (1 + gamma + alpha) .* log(T0 ./ T) - gamma .* sigma .* (Z);
normalizedData = exp(log(data) - altTerm);

end
