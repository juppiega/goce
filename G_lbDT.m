function [result] = G_lbDT(a, S, numBiases)

k = numBiases + 1;

S.latitudeTerm = a(k+1)*S.P10 + a(k+2)*S.P20 + a(k+3)*S.P30 + a(k+4)*S.P40 + ...
                 a(k+5)*S.FA.*S.P10 + a(k+6)*S.FA.*S.P20 + a(k+7)*S.FA.*S.P30 + ...
                 a(k+8)*S.F.*S.P10 + a(k+9)*S.F.*S.P20 + a(k+10)*S.F.*S.P30;
k = k + 10;

S.solarTerm = a(k+1)*S.F + a(k+2)*S.F2 + a(k+3)*S.FA + a(k+4)*S.FA2 + a(k+5)*S.FtimesFA;
k = k + 5;

S.annual = (a(k+1) + a(k+2)*S.P10 + a(k+3)*S.P20 + a(k+4)*S.P30 + a(k+5)*S.P40 + a(k+6)*S.FA + a(k+7)*S.F) .* ...
           (a(k+8)*sin(S.yv) + a(k+9)*cos(S.yv) + a(k+10)*sin(2*S.yv) + a(k+11)*cos(2*S.yv));
k = k + 11;

dPy = k + 7;
S.diurnal = ((a(k+1)*S.P11 + a(k+2)*S.P31 + a(k+3)*S.FA + a(k+4)*S.FA.^2) + (a(k+5)*S.P11 + a(k+6)*S.P21).*(cos(S.yv-pi*a(dPy)))).*cos(S.dv) + ...
            ((a(k+8)*S.P11 + a(k+9)*S.P31 + a(k+10)*S.FA + a(k+11)*S.FA.^2) + (a(k+12)*S.P11 + a(k+13)*S.P21).*(cos(S.yv-pi*a(dPy)))).*sin(S.dv);
k = k + 13;

S.semidiurnal = ((a(k+1)*S.P22 + a(k+2)*S.P32 + a(k+3)*S.FA + a(k+4)*S.FA.^2) + (a(k+5)*S.P32).*cos(S.yv-pi*a(dPy))).*cos(2*S.dv) + ...
                ((a(k+6)*S.P22 + a(k+7)*S.P32 + a(k+8)*S.FA + a(k+9)*S.FA.^2) + (a(k+10)*S.P32).*cos(S.yv-pi*a(dPy))).*sin(2*S.dv);
k = k + 10;
            
result = S.latitudeTerm + S.solarTerm + S.annual + S.diurnal + S.semidiurnal;

end
