function combineAE()

load supermag_2001_2015_south.mat
AU_s = AU; AL_s = AL; t_s = timestamps;

aeFiles = dir('*_north.mat');
AU_n = []; AL_n = []; t_n = [];
for i = 1:length(aeFiles)
    load(aeFiles(i).name)
    AU_n = [AU_n; AU];
    AL_n = [AL_n; AL];
    t_n = [t_n; timestamps];
end

t_comb = intersect(t_n, t_s);
t_full = (0:1:(t_comb(end)-t_comb(1))*24*60)/24/60 + t_comb(1);
AU_s = interp1(t_s, AU_s, t_full);
AL_s = interp1(t_s, AL_s, t_full);
AU_n = interp1(t_n, AU_n, t_full);
AL_n = interp1(t_n, AL_n, t_full);

AE = max(AU_s,AU_n) - min(AL_s,AL_n);

erroneous = AE > 1e4;
t_corr = t_full(~erroneous);
AE_corr = AE(~erroneous);
AE = interp1(t_corr,AE_corr,t_full);

dAE = abs(diff(AE,2));
spike_ind = find(dAE > 1e3);
spikes_exist = true;
k = 1;
while spikes_exist
    i = spike_ind(k);
    if median(AE(i-100:i+100)) < 900 
        %figure; subplot(1,2,1); plot(t_full(i-100:i+100), AE(i-100:i+100))
        ind = ismember(1:length(AE), i-1:i+5);
        t_corr = t_full(~ind);
        AE_corr = AE(~ind);
        AE = interp1(t_corr,AE_corr,t_full);
        %subplot(1,2,2); plot(t_full(i-100:i+100), AE(i-100:i+100))
        dAE = abs(diff(AE,2));
        spike_ind = find(dAE > 1e3);
        k = 0;
    end
    k = k + 1;
    if k > length(spike_ind)
        spikes_exist = false;
    end
end

% for k = 1:length(spike_ind)
%     i = spike_ind(k);
%     figure;
%     plot(t_full(i-100:i+100), AE(i-100:i+100))
% end

timestamps_AE_global = t_full';
AE_global = AE';
save('AE_global.mat','AE_global','timestamps_AE_global')

end