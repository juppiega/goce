function [F] = evalDT(S, coeff)

F = coeff(1) + G_lbDT(coeff, S);

end