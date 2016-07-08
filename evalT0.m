function [F] = evalT0(S, coeff)

F = clamp(300, coeff(1) * (1 + clamp(-0.5, G_lbT0(coeff, S, 0), 2)), 1000);

end