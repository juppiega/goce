function result = G_storm(a, S)

k = 0;

% Activity = a(k+1)*S.aeInt(:,5) + a(k+2)*S.aeInt(:,7) + ...
%            a(k+3)*S.aeInt(:,2).^2 + a(k+4)*S.aeInt(:,3).^2 + ...
%            a(k+5)*S.aeInt(:,4).^2 + a(k+6)*S.aeInt(:,7).^2;
% k = k + 6;

Activity = a(k+1)*S.aeInt(:,1) + a(k+2)*S.aeInt(:,3) + a(k+3)*S.aeInt(:,4)  + ...
            a(k+4)*S.aeInt(:,5) + a(k+5)*S.aeInt(:,7) + ...
            a(k+6)*S.aeInt(:,3).*S.aeInt(:,4) + a(k+7)*S.aeInt(:,3).*S.aeInt(:,5) + ...
            a(k+8)*S.aeInt(:,5).*S.aeInt(:,7) + a(k+9)*S.aeInt(:,3).^2 + ...
            (k+10)*S.aeInt(:,4).^2 + a(k+11)*S.aeInt(:,7).^2;
k = k + 11;

mag_lat = (a(k+1)+a(k+2)*S.P20+a(k+3)*S.P40 + (a(k+4)*S.P10+a(k+5)*S.P30).*(1+a(k+6)*S.FA+a(k+7)*S.FA.^2).*cos(S.dv-a(k+8))+...
           a(k+9)*S.FA+a(k+10)*S.FA.^2) .* Activity;
k = k + 10;

mag_lon = (a(k+1)*S.P21 + a(k+2)*S.P41 + a(k+3)*S.P61).*(1+a(k+4)*S.P10.*cos(S.yv+a(k+5))).*cos(S.lv-a(k+6)).*Activity;
k = k + 6;

mag_lst = (a(k+1)*S.P11+a(k+2)*S.P31+a(k+3)*S.P51).*cos(S.dv-a(k+4)).*Activity;
k = k + 4;

result = mag_lat + mag_lon + mag_lst;

% result = geomParametrization(S, a(k+1:k+6), S.aeInt(:,1)) +...
%                 geomParametrization(S, a(k+7:k+12), S.aeInt(:,2)) +...
%                 geomParametrization(S, a(k+13:k+18), S.aeInt(:,3)) +...
%                 geomParametrization(S, a(k+19:k+24), S.aeInt(:,5)) +...
%                 geomParametrization(S, a(k+25:k+30), S.aeInt(:,6));
% if isfield(S,'coeffInd') && length(S.coeffInd == 543) ~= 0; 
%     aaaa = 1; 
% end
end