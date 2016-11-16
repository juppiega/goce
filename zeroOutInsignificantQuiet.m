function [optCoeff, ptf] = zeroOutInsignificantQuiet(optCoeff, ptf, quietInd, paramErrors, significanceTol, S)

if iscolumn(optCoeff); optCoeff = optCoeff'; end
if iscolumn(quietInd); quietInd = quietInd'; end
if iscolumn(paramErrors); paramErrors = paramErrors'; end

%optCoeff = optCoeff(ismember(1:length(optCoeff), intersect(quietInd, S.coeffInd)));
%paramErrors = paramErrors(ismember(quietInd, S.coeffInd));
relErr = paramErrors./optCoeff(quietInd);
errTol = 1 - significanceTol;
    function s = significant(ind)
        s = relErr(ind) <= errTol;
    end

% Constant coefficient = first coeff.
ptf = [ptf, S.coeffInd(1)];

k = S.coeffInd(1) + S.numBiases; % Biases newer refitted.

% Latitude terms.
i = k+1:k+14;
qind = quietInd(i);
ptf = [ptf, qind(relErr(i) <= errTol)];
k = k + 14;

% % Solar activity terms.
i = k+1:k+5;
qind = quietInd(i);
ptf = [ptf, qind(relErr(i) <= errTol)];
k = k + 5;

% Annual symmetric terms.
i = k+1:k+32;
qind = quietInd(i);
ptf = [ptf, qind([4,9,14,17,22,27,31])];
% Solar factors:
qsolar = qind([5,6,10,11,15,18,23,24,28,29,32]);
ptf = [ptf, qsolar(significant(qsolar))];
% Latitude factors:
qlat = qind([1:3,7:8,12:13,16,19:21,25:26,30]);
ptf = [ptf, qlat(significant(qlat))];
    function rmAnnual(latInd, allInd)
        if ~any(significant(qind(latInd))); ptf(ismember(ptf,qind(allInd))) = []; end
    end
rmAnnual(1:3, 1:6);
rmAnnual(7:8, 7:11);
rmAnnual(12:13, 12:15);
rmAnnual(16, 16:18);
rmAnnual(19:21, 19:24);
rmAnnual(25:26, 25:29);
rmAnnual(30, 30:32);
k = k + 32;

% Diurnal
i = k+1:k+21;
dPy = k + 11;
qind = quietInd(i);
ptf = [ptf, qind(significant(qind(1:10))), qind(significant(qind(12:21)))];
includeDPy = false;
if any(significant(qind([8:10, 19:21]))); includeDPy = true; end
k = k + 21;

% Semidiurnal
i = k+1:k+16;
qind = quietInd(i);
ptf = [ptf, qind(significant(qind(1:16)))];
if any(significant(qind([7:8, 15:16]))); includeDPy = true; end
k = k + 16;

% Terdiurnal
i = k+1:k+8;
qind = quietInd(i);
ptf = [ptf, qind(significant(qind(1:8)))];
if any(significant(qind([3:4, 7:8]))); includeDPy = true; end
if includeDPy; ptf = [ptf, dPy]; end
k = k + 8;

% Quaterdiurnal
i = k+1:k+2;
qind = quietInd(i);
ptf = [ptf, qind(significant(qind(1:2)))];
k = k + 2;

% Longitudinal
i = k+1:k+13;
qind = quietInd(i);
dPy = k + 7;
ptf = [ptf, qind(significant(qind(2:6))), qind(significant(qind(9:13)))];
if any(significant(qind([5:6, 12:13]))); ptf = [ptf, dPy]; end
rmAnnual(2:6, 1:6);
rmAnnual(9:13, 8:13);
S.longitudinal = (1.0 + a(k+1)*S.FA).*(a(k+2)*S.P21+a(k+3)*S.P41+a(k+4)*S.P61 + (a(k+5)*S.P11+a(k+6)*S.P31).*cos(S.yv-pi*dPy)).*cos(S.lv)+...
                 (1.0 + a(k+8)*S.FA).*(a(k+9)*S.P21+a(k+10)*S.P41+a(k+11)*S.P61 + (a(k+12)*S.P11+a(k+13)*S.P31).*cos(S.yv-pi*dPy)).*sin(S.lv);
k = k + 13;

ptf = sort(ptf);

ocThis = optCoeff(S.coeffInd);
ocThis(setdiff(S.coeffInd, ptf)) = 0;
optCoeff(S.coeffInd) = ocThis;

end