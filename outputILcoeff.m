function [] = outputILcoeff()

load optCoeff.mat
numQuiet = 112;
numStorm = max(TexInd) - numQuiet;
N = max(TexInd);
numBiases = struct('O', 5, 'N2', 6,...
     'He', 5, 'Ar', 2, 'O2', 0);

OInd(2:1+numBiases.O) = [];
N2Ind(2:1+numBiases.N2) = [];
HeInd(2:1+numBiases.He) = [];
%ArInd(2:1+numBiases.Ar) = [];

optCoeff = optCoeff([TexInd, OInd, N2Ind, HeInd, O2Ind]);

OInd = TexInd(length(TexInd)) + (1:N);
N2Ind = OInd(end) + (1:N);
HeInd = N2Ind(end) + (1:N);
%ArInd = HeInd(end) + (1:N);
O2Ind = HeInd(end) + 1;

file = fopen('IL2017_coeff_mod.f90', 'w');

fprintf(file,'%s\n','module IL2017_coeff_mod');
fprintf(file,'%s\n','    implicit none');

writeArray(file,'coeff',optCoeff);
writeArray(file,'dTCoeffs',dTCoeffs);
writeArray(file,'T0Coeffs',T0Coeffs);

fprintf(file,'%s\n','    integer :: i');
writeIndices(file, 'TexInd', TexInd)
writeIndices(file, 'OInd', OInd)
writeIndices(file, 'N2Ind', N2Ind)
writeIndices(file, 'HeInd', HeInd)
%writeIndices(file, 'ArInd', ArInd)
writeIndices(file, 'O2Ind', O2Ind)
fprintf(file,'\n');
fprintf(file,'%s\n','end module IL2017_coeff_mod');

fclose(file);


end

function writeArray(file, name, array)

fprintf(file,'\n');
fprintf(file,'    %s%s%s%d%s%16.13e%s\n','real(kind = 8), parameter :: ',name,'(',length(array),') = [',array(1),',&');
for i = 2:length(array)-1
    fprintf(file,'                                 %16.13e%s\n',array(i),',&');
end
fprintf(file,'                                 %16.13e%s\n',array(end),']');
fprintf(file,'\n');

end

function writeIndices(file, name, ind)

fprintf(file,'    %s%s%s%d%s%d%s%d%s\n','integer, parameter :: ',name,'(',length(ind),') = (/(i,i=',min(ind),',',max(ind),')/)');

end