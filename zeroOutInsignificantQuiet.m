function [optCoeff, paramsToFit] = zeroOutInsignificantQuiet(optCoeff, quietInd, paramErrors, significanceTol, S)

optCoeff = optCoeff(ismember(1:length(optCoeff), intersect(quietInd, S.coeffInd)));
paramErrors = paramErrors(ismember(quietInd, S.coeffInd));

k = 0;

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

% Annual symmetric terms.
S.annual = (a(k+1) + a(k+2)*S.P20 + a(k+3)*S.P40).*cos(S.yv-pi*a(k+4)).*(1+a(k+5)*S.FA+a(k+6)*S.FA.^2) + ...
           (a(k+7) + a(k+8)*S.P20).*cos(2*(S.yv-pi*a(k+9))).*(1+a(k+10)*S.FA+a(k+11)*S.FA.^2) + ...
           (a(k+12) + a(k+13)*S.P20).*cos(3*(S.yv-pi*a(k+14))).*(1+a(k+15)*S.FA) + ...
           a(k+16)*cos(4*(S.yv-pi*a(k+17))).*(1+a(k+18)*S.FA);
k = k + 18;
% Annual asymmetric
S.annual = S.annual + (a(k+1)*S.P10 + a(k+2)*S.P30 + a(k+3)*S.P50).*cos(S.yv-pi*a(k+4)).*(1+a(k+5)*S.FA+a(k+6)*S.FA.^2)+...
           (a(k+7)*S.P10 + a(k+8)*S.P30).*cos(2*(S.yv-pi*a(k+9))).*(1+a(k+10)*S.FA+a(k+11)*S.FA.^2) + ...
           a(k+12)*S.P10.*cos(3*(S.yv-pi*a(k+13))).*(1+a(k+14)*S.FA);
           
%a(k+1:k+15) = a(k+1:k+15) - 1;
k = k + 14;

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

dPy = a(k+7);
S.longitudinal = (1.0 + a(k+1)*S.FA).*(a(k+2)*S.P21+a(k+3)*S.P41+a(k+4)*S.P61 + (a(k+5)*S.P11+a(k+6)*S.P31).*cos(S.yv-pi*dPy)).*cos(S.lv)+...
                 (1.0 + a(k+8)*S.FA).*(a(k+9)*S.P21+a(k+10)*S.P41+a(k+11)*S.P61 + (a(k+12)*S.P11+a(k+13)*S.P31).*cos(S.yv-pi*dPy)).*sin(S.lv);
k = k + 13;

end