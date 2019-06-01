function [eqOfTime] = computeEquationOfTime (doy)

n = doy - 1.5;
g = mod(357.528 + 0.9856003*n, 360) * pi / 180;

eqOfTime = (-7.659*sin(g) + 9.863*sin(2*g + 3.5932))/60; % hours

end