function [result] = G(a, S)

k = S.numBiases + 1; % Counter, which helps adding terms.

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

dPh = k + 11;
S.diurnal = ((a(k+1)*S.P11 + a(k+2)*S.P31 + a(k+3)*S.P51 + a(k+4)*S.P71 + a(k+5)*S.FA + a(k+6)*S.FA.^2 + a(k+7)*(S.F - S.FA)) + (a(k+8)*S.P11 + a(k+9)*S.P21 + a(k+10)*S.P31).*(cos(S.yv-pi*a(dPh)))).*cos(S.dv) + ...
            ((a(k+12)*S.P11 + a(k+13)*S.P31 + a(k+14)*S.P51 + a(k+15)*S.P71 + a(k+16)*S.FA + a(k+17)*S.FA.^2 + a(k+18)*(S.F - S.FA)) + (a(k+19)*S.P11 + a(k+20)*S.P21 + a(k+21)*S.P31).*(cos(S.yv-pi*a(dPh)))).*sin(S.dv);
%a(k+1:k+21) = a(k+1:k+21) - 1;
k = k + 21;

S.semidiurnal = (a(k+1)*S.P22 + a(k+2)*S.P32 + a(k+3)*S.P52 + a(k+4)*S.FA + a(k+5)*S.FA.^2 + a(k+6)*(S.F-S.FA) + (a(k+7)*S.P32 + a(k+8)*S.P52).*cos(S.yv-pi*a(dPh))).*cos(2*S.dv) + ...
                (a(k+9)*S.P22 + a(k+10)*S.P32 + a(k+11)*S.P52 + a(k+12)*S.FA + a(k+13)*S.FA.^2 + a(k+14)*(S.F-S.FA) + (a(k+15)*S.P32 + a(k+16)*S.P52).*cos(S.yv-pi*a(dPh))).*sin(2*S.dv);
%a(k+1:k+16) = a(k+1:k+16) - 1;
k = k + 16;

S.terdiurnal = (a(k+1)*S.P33 + a(k+2)*S.P53 + (a(k+3)*S.P43 + a(k+4)*S.P63).*cos(S.yv-pi*a(dPh))).*cos(3*S.dv) + ...
               (a(k+5)*S.P33 + a(k+6)*S.P53 + (a(k+7)*S.P43 + a(k+8)*S.P63).*cos(S.yv-pi*a(dPh))).*sin(3*S.dv);
%a(k+1:k+8) = a(k+1:k+8) - 1;
k = k + 8;

S.quaterdiurnal = a(k+1)*S.P44.*cos(4*S.dv) + a(k+2)*S.P44.*sin(4*S.dv);
%a(k+1:k+2) = a(k+1:k+2) - 1;
k = k + 2;

dPh = k + 14;
S.geomagnetic = sum(bsxfun(@times, [a(k+1), a(k+2), a(k+3), a(k+4), a(k+5), a(k+6), a(k+7), a(k+8)], S.aeInt(:,1:end-1)),2) .* (a(k+9) + a(k+10)*S.P40 + a(k+11)*S.P60 + (a(k+12)*S.P10 + a(k+13)*S.P30).*cos(S.yv-pi*a(dPh))) + ...
                (a(k+15)*S.aeInt(:,1) + a(k+16)*S.aeInt(:,3) + a(k+17)*S.aeInt(:,6)) .* (a(k+18)*S.P11 + (a(k+19)*S.P21 + a(k+20)*S.P31).*cos(S.yv-pi*a(dPh))).*cos(S.dv) + ...
                (a(k+21)*S.aeInt(:,1) + a(k+22)*S.aeInt(:,3) + a(k+23)*S.aeInt(:,6)) .* (a(k+24)*S.P11 + (a(k+25)*S.P21 + a(k+26)*S.P31).*cos(S.yv-pi*a(dPh))).*sin(S.dv);
%a(k+1:k+26) = a(k+1:k+26) - 1;
k = k + 26;
            
result = S.latitudeTerm + S.solarTerm + S.annual + S.diurnal + S.semidiurnal + S.terdiurnal + S.quaterdiurnal + S.geomagnetic;

end