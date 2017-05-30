function [Tex, rho, O, N2, He] = computeDtm(S)

N = length(S.latitude);
Tex = zeros(N,1);
rho = zeros(N,1);
O = zeros(N,1);
N2 = zeros(N,1);
He = zeros(N,1);

dtm2013_mex();

targetCount = round(N / 10000);
barWidth = 50;
p = TimedProgressBar( targetCount, barWidth, ...
                    'Running DTM, ETA ', ...
                    '. Now at ', ...
                    'Completed in ' );


load F30toF10
S.F = pF(1) * S.F + pF(2);
S.FA = pF(1) * S.FA + pF(2);
for i = 1:N
    [Tex(i),~,rho(i),~,~,d] = dtm2013_mex(S.doy(i), S.altitude(i), S.latitude(i), S.longitude(i), ...
        S.solarTime(i), S.F(i), S.FA(i), S.ap3h(i), S.Ap(i));
    O(i) = d(3); N2(i) = d(4); He(i) = d(2);
    if mod(i, 10000) == 0
        p.progress;
    end
end
p.stop;

rho = rho * 1E3;
u2g = 1.6605402e-24;
O = O / (16 * u2g);
N2 = N2 / (28 * u2g);
He = He / (4 * u2g);


end