function [F] = evalMajorSpecies(S, coeff)

F = exp(coeff(1) + G_major(coeff, S));

end