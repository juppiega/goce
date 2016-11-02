function [F] = evalDT(S, coeff)

F = coeff(1) * (1 + clamp(-0.9, G_quiet(coeff, S), 9));

end