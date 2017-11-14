function result = G_storm(a, S)

k = 0;

% Activity = a(k+1)*S.aeInt(:,5) + a(k+2)*S.aeInt(:,7) + ...
%            a(k+3)*S.aeInt(:,2).^2 + a(k+4)*S.aeInt(:,3).^2 + ...
%            a(k+5)*S.aeInt(:,4).^2 + a(k+6)*S.aeInt(:,7).^2;
% k = k + 6;

% Activity = a(k+1)*S.aeInt(:,1) + a(k+2)*S.aeInt(:,3) + a(k+3)*S.aeInt(:,4)  + ...
%             a(k+4)*S.aeInt(:,5) + a(k+5)*S.aeInt(:,7) + ...
%             a(k+6)*S.aeInt(:,3).*S.aeInt(:,4) + a(k+7)*S.aeInt(:,3).*S.aeInt(:,5) + ...
%             a(k+8)*S.aeInt(:,5).*S.aeInt(:,7) + a(k+9)*S.aeInt(:,3).^2 + ...
%             a(k+10)*S.aeInt(:,4).^2 + a(k+11)*S.aeInt(:,7).^2;
% k = k + 11;
% 
% mag_lat = (a(k+1)+a(k+2)*S.P20+a(k+3)*S.P40 + (a(k+4)*S.P10+a(k+5)*S.P30).*(1+a(k+6)*S.FA+a(k+7)*S.FA.^2).*cos(S.yv-a(k+8))+...
%            a(k+9)*S.FA+a(k+10)*S.FA.^2) .* Activity;
% k = k + 10;
% 
% mag_lon = (a(k+1)*S.P21 + a(k+2)*S.P41 + a(k+3)*S.P61).*(1+a(k+4)*S.P10.*cos(S.yv+a(k+5))).*cos(S.lv-a(k+6)).*Activity;
% k = k + 6;
% 
% mag_lst = (a(k+1)*S.P11+a(k+2)*S.P31+a(k+3)*S.P51).*cos(S.dv-a(k+4)).*Activity;
% k = k + 4;
% 
% result = mag_lat + mag_lon + mag_lst;

% Activity = a(k+1)*S.aeInt(:,1) + a(k+2)*S.aeInt(:,3) + a(k+3)*S.aeInt(:,4)  + ...
%             a(k+4)*S.aeInt(:,5) + a(k+5)*S.aeInt(:,7) + ...
%             a(k+6)*S.aeInt(:,3).*S.aeInt(:,4) + a(k+7)*S.aeInt(:,3).*S.aeInt(:,5) + ...
%             a(k+8)*S.aeInt(:,5).*S.aeInt(:,7) + a(k+9)*S.aeInt(:,3).^2 + ...
%             a(k+10)*S.aeInt(:,4).^2 + a(k+11)*S.aeInt(:,7).^2;
% k = k + 11;
% 
% mag_lat = (a(k+1)+a(k+2)*S.P20+a(k+3)*S.P40 + (a(k+4)*S.P10+a(k+5)*S.P30).*(1+a(k+6)*S.FA+a(k+7)*S.FA.^2).*cos(S.yv-a(k+8))+...
%            a(k+9)*S.FA+a(k+10)*S.FA.^2) .* Activity;
% k = k + 10;
tauVec = (1:24)';
lagHours = a(k+1);
aeInt = interp1(tauVec,S.aeInt',clamp(1.0,lagHours,24.0))';  
k = k + 1;

%mag_lat = 1E-5*(a(k+1) + a(k+2)*S.P20 + a(k+3)*S.P40).*aeInt.*(1 +1E-3*a(k+4)*aeInt); quadratic
%mag_lat = 1E-5*(a(k+1) + a(k+2)*S.P20 + a(k+3)*S.P40 + (a(k+4)*S.P10 + a(k+5)*S.P30).*cos(S.yv-1E-2*a(k+6)) ).*aeInt;% + 1E-3*a(k+4);
% mag_lat = 1E-5*a(k+1)*aeInt;

an = a(k+8);
san = a(k+10);
order_0 = (a(k+1)+a(k+2)*S.mP20+a(k+3)*S.mP40+a(k+4)*S.mP60 + (a(k+5)*S.mP10+a(k+6)*S.mP30+a(k+7)*S.mP50).*...
    (cos(S.yv-an)+a(k+9)*cos(2*(S.yv-san)))).*(1+a(k+11)*S.F).*aeInt;
k = k + 11;%12

order_1 = (a(k+1)*S.mP11+a(k+2)*S.mP31+a(k+3)*S.mP51 + (a(k+4)*S.mP21+a(k+5)*S.mP41+a(k+6)*S.mP61).*...
    (cos(S.yv-an)+a(k+7)*cos(2*(S.yv-san)))).*(1+a(k+8)*S.F).*aeInt.*cos(S.dv_mag-a(k+9));
k = k + 9;%21

order_2 = (a(k+1)*S.mP22+a(k+2)*S.mP42+a(k+3)*S.mP62 + (a(k+4)*S.mP32+a(k+5)*S.mP52).*...
    (cos(S.yv-an)+a(k+6)*cos(2*(S.yv-san)))).*(1+a(k+7)*S.F).*aeInt.*cos(2*(S.dv_mag-a(k+8)));
k = k + 8;%29

%order_3 = (a(k+1)*S.mP33+a(k+2)*S.mP53 + (a(k+3)*S.mP43+a(k+4)*S.mP63).*...
%    (cos(S.yv-an)+a(k+5)*cos(2*(S.yv-san)))).*(1+a(k+6)*S.F).*aeInt.*cos(3*(S.dv_mag-a(k+7)));
%k = k + 7;%36

result = order_0 + order_1 + order_2;% + order_3;
% result = geomParametrization(S, a(k+1:k+6), S.aeInt(:,1)) +...
%                 geomParametrization(S, a(k+7:k+12), S.aeInt(:,2)) +...
%                 geomParametrization(S, a(k+13:k+18), S.aeInt(:,3)) +...
%                 geomParametrization(S, a(k+19:k+24), S.aeInt(:,5)) +...
%                 geomParametrization(S, a(k+25:k+30), S.aeInt(:,6));
% if isfield(S,'coeffInd') && length(S.coeffInd == 543) ~= 0; 
%     aaaa = 1; 
% end
end