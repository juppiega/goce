function result = bulge(a, S)

result = a(1) + a(2)*S.P10 + S.P11.*(a(3)*cos(S.dv) + a(4)*sin(S.dv)) + ...
         a(5)*S.P20 + S.P21.*(a(6)*cos(S.dv) + a(7)*sin(S.dv)) + ...
         S.P22.*(a(8)*cos(2*S.dv) + a(9)*sin(2*S.dv));

end