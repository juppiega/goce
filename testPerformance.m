function [c] = testPerformance()

a = ones(30000000, 20);
b = 2*a;
z = zeros(size(a));

tic;
c = (a + b .* ((a - 2 * b)) + (b ./ a) ./ b + a .* a ./b) ./ (a + b .* ((a - 2 * b)) + (b ./ a) ./ b + a .* a ./b);
toc;

end