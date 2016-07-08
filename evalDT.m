function [F] = evalDT(S, coeff)

F = coeff(1) * (1 + clamp(-0.9, G_lbDT(coeff, S, 0), 9));

end