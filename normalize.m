function valueNoBg = normalize(valueNoBg, value)
% [aeNoBg, averagedDensityNoBg] = normalize(aeNoBg, averagedDensityNoBg)

valueNoBg = mean(value) / mean(valueNoBg) * valueNoBg;
end