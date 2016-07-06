function [F] = evalMinorSpecies(S, coeff)

F = exp(coeff(1) + G_minor(coeff, S));

end