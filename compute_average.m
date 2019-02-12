function [x_p, mean_y] = compute_average(x, y, dx)

x_p = min(x):dx:max(x);
x_p = [x_p,max(x)];

for i = 1:length(x_p)-1
    ind = x_p(i) < x & x < x_p(i+1);
    mean_y(i) = mean(y(ind));
end

x_p = mean([x_p(1:end-1); x_p(2:end)]);

figure;
plot(x_p, mean_y)

end