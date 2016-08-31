function result = fourier4(S, a)

result = a(1) + a(2)*sin(S.yv) + a(3)*cos(S.yv) + a(4)*sin(2*S.yv) + a(5)*cos(2*S.yv);

end