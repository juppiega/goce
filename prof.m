alts = 120 : 10 : 800;
for i = 1:length(alts)
    [~,tn(i),~] = jb2008_mex(2457183, alts(i), 70, 25, 100, 100, 100, 100, 100, 100, 100, 100, 0);
end

plot(tn, alts)