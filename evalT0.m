function [F] = evalT0(S, coeff)

F = coeff(1) * (1 + clamp(-0.9, G_quiet(coeff(2:end), S), 9));

end