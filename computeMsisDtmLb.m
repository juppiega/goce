function [T0_msis, dT_msis, T0_dtm, dT_dtm] = computeMsisDtmLb(S)

% Annual parameter.
if ~isfield(S, 'doy')
    [yr,~,~,~,~,~] = datevec(S.timestamps);
    yearVec = [yr, repmat([1,1,0,0,0], length(yr), 1)];
    S.doy = S.timestamps - datenum(yearVec) + 1;
end

z0 = 130;
dz = 0.1;

N = length(S.latitude);
T0_msis = zeros(N, 1);
dT_msis = zeros(N, 1);
T0_dtm = zeros(N, 1);
dT_dtm = zeros(N, 1);

targetCount = round(N / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running MSIS, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );
dtm2013_mex();

for i = 1:N
    [~, ~, ~,~,~,~,~,~,~,~,T1] = nrlmsise_mex(S.doy(i),43200,z0-dz,S.latitude(i),S.longitude(i),S.solarTime(i),...
        S.FA(i),S.F(i),S.Ap(i),S.apNow(i),S.ap3h(i),S.ap6h(i),S.ap9h(i),S.ap12To33h(i),S.ap36To57h(i));
    [~, ~, ~,~,~,~,~,~,~,~,T2] = nrlmsise_mex(S.doy(i),43200,z0,S.latitude(i),S.longitude(i),S.solarTime(i),...
        S.FA(i),S.F(i),S.Ap(i),S.apNow(i),S.ap3h(i),S.ap6h(i),S.ap9h(i),S.ap12To33h(i),S.ap36To57h(i));
    [~, ~, ~,~,~,~,~,~,~,~,T3] = nrlmsise_mex(S.doy(i),43200,z0+dz,S.latitude(i),S.longitude(i),S.solarTime(i),...
        S.FA(i),S.F(i),S.Ap(i),S.apNow(i),S.ap3h(i),S.ap6h(i),S.ap9h(i),S.ap12To33h(i),S.ap36To57h(i));
    dT_msis(i) = (T3-T1) / (2*dz);
    T0_msis(i) = T2;
    
    [~,T0,~,~,~,~,dT0] = dtm2013_mex(S.doy(i), z0, S.latitude(i), S.longitude(i), ...
    S.solarTime(i), S.F(i), S.FA(i), S.ap3h(i), S.Ap(i));
    [~,Tupper,~,~,~,~,~] = dtm2013_mex(S.doy(i), z0+dz, S.latitude(i), S.longitude(i), ...
    S.solarTime(i), S.F(i), S.FA(i), S.ap3h(i), S.Ap(i));
    [~,Tlower,~,~,~,~,~] = dtm2013_mex(S.doy(i), z0-dz, S.latitude(i), S.longitude(i), ...
    S.solarTime(i), S.F(i), S.FA(i), S.ap3h(i), S.Ap(i));

    dT_dtm(i) = (Tupper - Tlower)/(2*dz);
    T0_dtm(i) = T0;
    
    if mod(i, 10000) == 0
        p.progress;
    end
end
p.stop;

end