function [F] = evalT0(S, coeff)

F = coeff(1) * (1 + clamp(-0.9, G_lbT0(coeff, S, 0), 9));

end