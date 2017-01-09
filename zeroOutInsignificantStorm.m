function [optCoeff, paramsToFit] = zeroOutInsignificantStorm(optCoeff, paramsToFit, stormInd, paramErrors, significance)

if iscolumn(optCoeff); optCoeff = optCoeff'; end
if iscolumn(stormInd); stormInd = stormInd'; end
if iscolumn(paramErrors); paramErrors = paramErrors'; end

relErr = abs(paramErrors./optCoeff(stormInd));
errTol = 1 - significanceTol;

    function s = anyIsSignificant(ind)
        for j = 1:length(ind)
            if relErr(stormInd == ind(j)) <= errTol
                s = true;
                return;
            end
        end
        s = false;
    end

    function s = signifInd(ind)
        s = zeros(size(ind));
        for j = 1:length(ind)
            if relErr(stormInd == ind(j)) <= errTol
                s(j) = ind(j);
            end
        end
        s(s == 0) = [];
    end

fieldChanges = stormInd(diff(stormInd) > 1);
beginInd = [stormInd(1), fieldChanges+1];
endInd = [fieldChanges, stormInd(end)];



end