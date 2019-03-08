function plotAE(d1, d2)

load AE_global.mat
load AE_original.mat

t1 = datenum(d1);
t2 = datenum(d2);
i_orig = t1 <= timestampsAeDatenum & timestampsAeDatenum < t2;
i_glob = t1 <= timestamps_AE_global & timestamps_AE_global < t2;

figure;
plot(timestampsAeDatenum(i_orig), ae(i_orig),timestamps_AE_global(i_glob), AE_global(i_glob)); datetick('x')
legend('Original','Global')
corr(ae(i_orig), AE_global(i_glob))

end