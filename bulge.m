function result = bulge(a, S)

result = a(1)*S.P10 + S.P11.*(a(2)*cos(S.dv) + a(3)*sin(S.dv)) + ...
         a(4)*S.P20 + S.P21.*(a(5)*cos(S.dv) + a(6)*sin(S.dv)) + ...
         S.P22.*(a(7)*cos(2*S.dv) + a(8)*sin(2*S.dv));

end