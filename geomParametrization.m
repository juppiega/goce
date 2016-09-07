function result = geomParametrization(S, a, aeInt)

result = (a(1) + a(2)*cos(S.yv-a(3)).*S.P10 + a(4)*S.P20 + a(5)*S.P40).*(1 + a(6)*S.FA).*aeInt;

end