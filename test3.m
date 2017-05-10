function [] = test3()

a = randi([1,101],6000,1);
b = randi([30,92],1000,1);
c = randi([10,50],700,1);
v = [a;b;c];
v(v==50) = [];

figure;
histogram(v, 20)
set(gca,'fontsize',15);
title('Järjestyslukujen jakauma','fontsize',15)
xlabel('Järjestysluku','fontsize',15)
xlim([0,101])

end