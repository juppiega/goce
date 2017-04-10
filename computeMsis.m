function [Tex, rho, O, N2, He, Ar, O2, T] = computeMsis(S)

% Annual parameter.
if ~isfield(S, 'doy')
    [yr,~,~,~,~,~] = datevec(S.timestamps);
    yearVec = [yr, repmat([1,1,0,0,0], length(yr), 1)];
    S.doy = S.timestamps - datenum(yearVec) + 1;
end

N = length(S.latitude);
Tex = zeros(N,1);
T = zeros(N,1);
rho = zeros(N,1);
O = zeros(N,1);
N2 = zeros(N,1);
He = zeros(N,1);
Ar = zeros(N,1);
O2 = zeros(N,1);

load F30toF10
S.F = pF(1) * S.F + pF(2);
S.FA = pF(1) * S.FA + pF(2);

targetCount = round(N / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running MSIS, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );

for i = 1:N
    [He(i), O(i), N2(i),O2(i),Ar(i),rho(i),~,~,~,Tex(i),T(i)] = nrlmsise_mex(S.doy(i),43200,S.altitude(i),S.latitude(i),S.longitude(i),S.solarTime(i),...
        S.FA(i),S.F(i),S.Ap(i),S.apNow(i),S.ap3h(i),S.ap6h(i),S.ap9h(i),S.ap12To33h(i),S.ap36To57h(i));
    if mod(i, 10000) == 0
        p.progress;
    end
end
p.stop;

rho = rho * 1E3;

end