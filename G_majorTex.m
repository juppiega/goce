function [result] = G_majorTex(a, S, numBiases)

k = numBiases + 1; % Counter, which helps adding terms.
result = G_quiet(a(k+1:k+111), S) + G_storm(a(k+112:end), S);

end