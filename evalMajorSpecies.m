function [F] = evalMajorSpecies(S, coeff, numBiases)

F = exp(coeff(1) + G_major(coeff, S, numBiases));

end