function [F] = evalTex(S, coeff)

F = coeff(1) * (1 + clamp(-0.9, G_major(coeff, S, 0), 9));

end