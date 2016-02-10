function [F] = evalSpecies(S, coeff)

F = exp(coeff(1) + G(coeff, S));

end