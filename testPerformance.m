function [c] = testPerformance()

a = ones(30000000, 20);
b = 2*a;
z = zeros(size(a));

tic;
c = a + b .* (exp(a - 2 * b)) + log(b ./ a) ./ b + a;
toc;

end