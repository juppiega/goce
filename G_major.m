function [result] = G_major(a, S, numBiases)

k = numBiases + 1; % Counter, which helps adding terms.

% Latitude terms.
S.latitudeTerm = a(k+1)*S.P10 + a(k+2)*S.P20 + a(k+3)*S.P30 + a(k+4)*S.P40 + a(k+5)*S.P50 + a(k+6)*S.P60 + ...
                 a(k+7)*S.FA.*S.P10 + a(k+8)*S.FA.*S.P20 + a(k+9)*S.FA.*S.P30 + a(k+10)*S.FA.*S.P40 + ...
                 a(k+11)*S.F.*S.P10 + a(k+12)*S.F.*S.P20 + a(k+13)*S.F.*S.P30 + a(k+14)*S.F.*S.P40;
%a(k+1:k+14) = a(k+1:k+14) - 1;
k = k + 14;

% % Solar activity terms.
S.solarTerm = a(k+1)*S.F + a(k+2)*S.F2 + a(k+3)*S.FA + a(k+4)*S.FA2 + a(k+5)*S.FtimesFA;
%a(k+1:k+5) = a(k+1:k+5) - 1;
k = k + 5;

% Annual terms.
S.annual = (a(k+1) + a(k+2)*S.P10 + a(k+3)*S.P20 + a(k+4)*S.P30 + a(k+5)*S.P40 + a(k+6)*S.FA + a(k+7)*S.F) .* ...
           (a(k+8)*sin(S.yv) + a(k+9)*cos(S.yv) + a(k+10)*sin(2*S.yv) + a(k+11)*cos(2*S.yv) + a(k+12)*sin(3*S.yv) + a(k+13)*cos(3*S.yv) + a(k+14)*sin(4*S.yv) + a(k+15)*cos(4*S.yv));
%a(k+1:k+15) = a(k+1:k+15) - 1;
k = k + 15;

dPy = k + 11;
S.diurnal = ((a(k+1)*S.P11 + a(k+2)*S.P31 + a(k+3)*S.P51 + a(k+4)*S.P71 + a(k+5)*S.FA + a(k+6)*S.FA.^2 + a(k+7)*(S.F - S.FA)) + (a(k+8)*S.P11 + a(k+9)*S.P21 + a(k+10)*S.P31).*(cos(S.yv-pi*a(dPy)))).*cos(S.dv) + ...
            ((a(k+12)*S.P11 + a(k+13)*S.P31 + a(k+14)*S.P51 + a(k+15)*S.P71 + a(k+16)*S.FA + a(k+17)*S.FA.^2 + a(k+18)*(S.F - S.FA)) + (a(k+19)*S.P11 + a(k+20)*S.P21 + a(k+21)*S.P31).*(cos(S.yv-pi*a(dPy)))).*sin(S.dv);
%a(k+1:k+21) = a(k+1:k+21) - 1;
k = k + 21;

S.semidiurnal = (a(k+1)*S.P22 + a(k+2)*S.P32 + a(k+3)*S.P52 + a(k+4)*S.FA + a(k+5)*S.FA.^2 + a(k+6)*(S.F-S.FA) + (a(k+7)*S.P32 + a(k+8)*S.P52).*cos(S.yv-pi*a(dPy))).*cos(2*S.dv) + ...
                (a(k+9)*S.P22 + a(k+10)*S.P32 + a(k+11)*S.P52 + a(k+12)*S.FA + a(k+13)*S.FA.^2 + a(k+14)*(S.F-S.FA) + (a(k+15)*S.P32 + a(k+16)*S.P52).*cos(S.yv-pi*a(dPy))).*sin(2*S.dv);
%a(k+1:k+16) = a(k+1:k+16) - 1;
k = k + 16;

S.terdiurnal = (a(k+1)*S.P33 + a(k+2)*S.P53 + (a(k+3)*S.P43 + a(k+4)*S.P63).*cos(S.yv-pi*a(dPy))).*cos(3*S.dv) + ...
               (a(k+5)*S.P33 + a(k+6)*S.P53 + (a(k+7)*S.P43 + a(k+8)*S.P63).*cos(S.yv-pi*a(dPy))).*sin(3*S.dv);
%a(k+1:k+8) = a(k+1:k+8) - 1;
k = k + 8;

S.quaterdiurnal = a(k+1)*S.P44.*cos(4*S.dv) + a(k+2)*S.P44.*sin(4*S.dv);
%a(k+1:k+2) = a(k+1:k+2) - 1;
k = k + 2;

% ATTEMPT #1
% AE_base = sum(bsxfun(@times, [a(k+1), a(k+2), a(k+3), a(k+4), a(k+5), a(k+6), a(k+7), a(k+8)], S.aeInt(:,1:end-1)),2);
% geom_symmetric = AE_base + (a(k+9)*S.mP20 + a(k+10)*S.mP40 + a(k+11)*S.mP60).*AE_base;
% geom_yearly = (a(k+12)*S.mP10 + a(k+13)*S.mP30 + a(k+14)*S.mP50 + a(k+15)*S.mP70).*AE_base.*cos(S.yv-pi*a(k+16));
% geom_lst = (a(k+17)*S.mP11 + a(k+18)*S.mP31 + a(k+19)*S.mP51).*AE_base.*cos(S.dv_mag-pi*a(k+20));
% S.geomagnetic = geom_symmetric + geom_yearly + geom_lst + a(k+21)*AE_base.^2;

% ATTEMPT #2
% AE_base = sum(bsxfun(@times, [a(k+1), a(k+2), a(k+3), a(k+4), a(k+5), a(k+6), a(k+7), a(k+8)], S.aeInt(:,1:end-1)),2);
% geom_symmetric = AE_base + (a(k+9)*S.mP20 + a(k+10)*S.mP40 + a(k+11)*S.mP60).*AE_base;
% geom_yearly = (a(k+12)*S.mP10 + a(k+13)*S.mP30 + a(k+14)*S.mP50 + a(k+15)*S.mP70).*AE_base.*cos(S.yv-pi*a(k+16));
% geom_lst = (a(k+17)*S.mP11 + a(k+18)*S.mP31 + a(k+19)*S.mP51).*AE_base.*cos(S.dv_mag-pi*a(k+20));
% S.geomagnetic = geom_symmetric + geom_yearly + geom_lst + a(k+21)*exp(a(k+22)*AE_base);


% ATTEMPT #3
AE_base = sum(bsxfun(@times, [a(k+1), a(k+2), a(k+3), a(k+4), a(k+5), a(k+6), a(k+7)], S.aeInt),2);
geom_symmetric = (a(k+8) + a(k+9)*S.mP20 + a(k+10)*S.mP40 + a(k+11)*S.mP60).*AE_base;
dPy = a(k+16);
geom_lon = (a(k+12) + a(k+13)*S.mP21 + a(k+14)*S.mP41).*(1+a(k+15)*cos(S.yv-pi*dPy)).*AE_base.*cos(S.lv-pi*a(k+17)) +...
           (a(k+18) + a(k+19)*S.mP32 + a(k+20)*S.mP52).*(1+a(k+21)*cos(S.yv-pi*dPy)).*AE_base.*cos(2*(S.lv-pi*a(k+22)));
geom_lst = (a(k+23) + a(k+24)*S.mP21 + a(k+25)*S.mP41).*(1+a(k+26)*cos(S.yv-pi*dPy)).*AE_base.*cos(S.dv_mag-pi*a(k+27)) +...
           (a(k+28) + a(k+29)*S.mP32 + a(k+30)*S.mP52).*(1+a(k+31)*cos(S.yv-pi*dPy)).*AE_base.*cos(2*(S.dv_mag-pi*a(k+32)));
geom_solar = (a(k+33) + a(k+34)*S.mP10.*cos(S.yv-pi*dPy) + a(k+35)*S.mP20).*AE_base.*(1./S.FA);
S.geomagnetic = geom_symmetric + geom_lon + geom_lst + geom_solar;

k = k + 35;
            
result = S.latitudeTerm + S.solarTerm + S.annual + S.diurnal + S.semidiurnal + S.terdiurnal + S.quaterdiurnal + S.geomagnetic;

end