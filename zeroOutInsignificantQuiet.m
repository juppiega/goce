function [optCoeff, ptf] = zeroOutInsignificantQuiet(optCoeff, ptf, quietInd, paramErrors, significanceTol, S, includePhaseSignificantLat)

if iscolumn(optCoeff); optCoeff = optCoeff'; end
if iscolumn(quietInd); quietInd = quietInd'; end
if iscolumn(paramErrors); paramErrors = paramErrors'; end

%optCoeff = optCoeff(ismember(1:length(optCoeff), intersect(quietInd, S.coeffInd)));
%paramErrors = paramErrors(ismember(quietInd, S.coeffInd));
relErr = abs(paramErrors./optCoeff(quietInd));
errTol = 1 - significanceTol;
    function s = anyIsSignificant(ind)
        for j = 1:length(ind)
            if relErr(quietInd == ind(j)) <= errTol
                s = true;
                return;
            end
        end
        s = false;
    end

    function s = signifInd(ind)
        s = zeros(size(ind));
        for j = 1:length(ind)
            if relErr(quietInd == ind(j)) <= errTol
                s(j) = ind(j);
            end
        end
        s(s == 0) = [];
    end

k = find(quietInd == S.coeffInd(1), 1);

% Constant coefficient = first coeff.
ptf = [ptf, S.coeffInd(1)];

k = k + S.numBiases; % Biases newer refitted.

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
ptf = [ptf, signifInd(qsolar)];
% Latitude factors:
qlat = qind([1:3,7:8,12:13,16,19:21,25:26,30]);
ptf = [ptf, signifInd(qlat)];
    function rmAnnual(latInd, allInd)
        if ~anyIsSignificant(qind(latInd)); ptf(ismember(ptf,qind(allInd))) = []; end
    end
rmAnnual(1:3, 1:6);
rmAnnual(7:8, 7:11);
rmAnnual(12:13, 12:15);
rmAnnual(16, 16:18);
rmAnnual(19:21, 19:24);
rmAnnual(25:26, 25:29);
rmAnnual(30, 30:32);
if nargin > 6 && includePhaseSignificantLat
    ptf = [ptf, signifInd(qind([4,9,14,17,22,27,31]))];
    if ~isempty(find(ptf == qind(4))) ptf = [ptf, qind(1:3)]; end
    if ~isempty(find(ptf == qind(9))) ptf = [ptf, qind(7:8)]; end
    if ~isempty(find(ptf == qind(14))) ptf = [ptf, qind(12:13)]; end
    if ~isempty(find(ptf == qind(17))) ptf = [ptf, qind(16)]; end
    if ~isempty(find(ptf == qind(22))) ptf = [ptf, qind(19:21)]; end
    if ~isempty(find(ptf == qind(27))) ptf = [ptf, qind(25:26)]; end
    if ~isempty(find(ptf == qind(31))) ptf = [ptf, qind(30)]; end      
end
k = k + 32;

% Diurnal
i = k+1:k+21;
dPy = quietInd(k + 11);
qind = quietInd(i);
ptf = [ptf, signifInd(qind(1:10)), signifInd(qind(12:21))];
includeDPy = false;
if anyIsSignificant(qind([8:10, 19:21])); includeDPy = true; end
k = k + 21;

% Semidiurnal
i = k+1:k+16;
qind = quietInd(i);
ptf = [ptf, signifInd(qind(1:16))];
if anyIsSignificant(qind([7:8, 15:16])); includeDPy = true; end
k = k + 16;

% Terdiurnal
i = k+1:k+8;
qind = quietInd(i);
ptf = [ptf, signifInd(qind(1:8))];
if anyIsSignificant(qind([3:4, 7:8])); includeDPy = true; end
if includeDPy; ptf = [ptf, dPy]; end
k = k + 8;

% Quaterdiurnal
i = k+1:k+2;
qind = quietInd(i);
ptf = [ptf, signifInd(qind(1:2))];
k = k + 2;

% Longitudinal
i = k+1:k+13;
qind = quietInd(i);
dPy = quietInd(k + 7);
ptf = [ptf, signifInd(qind(2:6)), signifInd(qind(9:13))];
if anyIsSignificant(qind([5:6, 12:13])); ptf = [ptf, dPy]; end
rmAnnual(2:6, 1:6);
rmAnnual(9:13, 8:13);
k = k + 13;

ptf = unique(ptf);

zeroOutInd = setdiff(S.coeffInd, ptf);
zeroOutInd = setdiff(zeroOutInd, S.coeffInd(2:S.numBiases+1));
optCoeff(zeroOutInd) = 0;

end