function [] = summariseSignificant (saveFolder)

load optCoeff.mat

Nclass = 9 + 3;
signMat = zeros(4,Nclass);

signMat(1,:) = computeSignificant(optCoeff(TexInd));
signMat(1,:) = computeSignificant(optCoeff(OInd));
signMat(1,:) = computeSignificant(optCoeff(N2Ind));
signMat(1,:) = computeSignificant(optCoeff(HeInd));

outputCell = cell(5,Nclass+1);
outputCell(2:end,2:end) = num2cell(signMat);
outputCell(1,1:end) = {'latitude', 'solar', 'annualSymm', 'annualAsymm', 'diurnal', 'semidiurnal', 'terdiurnal', 'quaterdiurnal',...
                        'longitudinal','S_0','S_1','S_2'};
outputCell(2:end,1) = {'Tex','O','N2','He'};

cell2csv([saveFolder, '/significant.csv'], outputCell);

end

function sign = computeSignificant(param, numBiases)

sign = zeros(1,12);

k = numBiases + 1;
latitude = k+1:k+14; k = k + 14;
solar = k+1:k+5; k = k + 5;
annualSymm = k+1:k+18; k = k + 18;
annualAsymm = k+1:k+14; k = k + 14;
diurnal = k+1:k+21; k = k + 21;
semidiurnal = k+1:k+16; k = k + 16;
terdiurnal = k+1:k+8; k = k + 8;
quaterdiurnal = k+1:k+2; k = k + 2;
longitudinal = k+1:k+13; k = k + 13;
S_0 = k+1:k+9; k = k + 9;
S_1 = k+1:k+8; k = k + 8;
S_2 = k+1:k+7; k = k + 7;

if any(param(latitude) ~= 0)
    sign(1) = 1;
end
if any(param(solar) ~= 0)
    sign(2) = 1;
end
if any(param(annualSymm) ~= 0)
    sign(3) = 1;
end
if any(param(annualAsymm) ~= 0)
    sign(4) = 1;
end
if any(param(diurnal) ~= 0)
    sign(5) = 1;
end
if any(param(semidiurnal) ~= 0)
    sign(6) = 1;
end
if any(param(terdiurnal) ~= 0)
    sign(7) = 1;
end
if any(param(quaterdiurnal) ~= 0)
    sign(8) = 1;
end
if any(param(longitudinal) ~= 0)
    sign(9) = 1;
end
if any(param(S_0) ~= 0)
    sign(10) = 1;
end
if any(param(S_1) ~= 0)
    sign(11) = 1;
end
if any(param(S_2) ~= 0)
    sign(12) = 1;
end

end
