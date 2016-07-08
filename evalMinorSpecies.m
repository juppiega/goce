function [F] = evalMinorSpecies(S, coeff, numBiases)

F = exp(coeff(1) + G_minor(coeff, S, numBiases));

end