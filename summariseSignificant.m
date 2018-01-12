function [] = summariseSignificant (saveFolder)

load optCoeff.mat

Nclass = 9 + 3;
signMat = cell(4,Nclass);

signMat(1,:) = computeSignificant(optCoeff(TexInd),0);
signMat(2,:) = computeSignificant(optCoeff(OInd),5);
signMat(3,:) = computeSignificant(optCoeff(N2Ind),6);
signMat(4,:) = computeSignificant(optCoeff(HeInd),5);

outputCell = cell(5,Nclass+1);
outputCell(2:end,2:end) = signMat;
outputCell(1,1:end) = {'latitude', 'solar', 'annualSymm', 'annualAsymm', 'diurnal', 'semidiurnal', 'terdiurnal', 'quaterdiurnal',...
                        'longitudinal','S_0','S_1','S_2'};
outputCell(2:end,1) = {'Tex','O','N2','He'};

cell2csv([saveFolder, '/significant.csv'], outputCell);

end

function sign = computeSignificant(param, numBiases)

sign = cell(1,12);

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

sign = fillCell(sign, param, latitude);
sign = fillCell(sign, param, solar);
sign = fillCell(sign, param, annualSymm);
sign = fillCell(sign, param, annualAsymm);
sign = fillCell(sign, param, diurnal);
sign = fillCell(sign, param, semidiurnal);
sign = fillCell(sign, param, terdiurnal);
sign = fillCell(sign, param, quaterdiurnal);
sign = fillCell(sign, param, longitudinal);
sign = fillCell(sign, param, S_0);
sign = fillCell(sign, param, S_1);
sign = fillCell(sign, param, S_2);

end

function sign = fillCell(sign, param, ind)

persistent k;
if isempty(k)
    k = 1;
end

if all(param(ind) ~= 0)
    sign{k} = 'F';
elseif any(param(ind) ~= 0)
    sign{k} = 'P';
end
k = k + 1;

end