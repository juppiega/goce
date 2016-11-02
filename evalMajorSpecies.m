function [F] = evalMajorSpecies(S, coeff, numBiases)

F = exp(coeff(1) + G_majorTex(coeff, S, numBiases));

end