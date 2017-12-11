function [paramsToFitShort] = zeroOutInsignificantStorm(optCoeff, relErr, paramsToFitShort, offset, significanceTol)

if iscolumn(optCoeff); optCoeff = optCoeff'; end
numStorm = length(optCoeff) / 4;
indOriginal = (numStorm*offset) + 1 : (offset+1)*numStorm;
optCoeff = optCoeff(indOriginal);
relErr = relErr(indOriginal);
ptf = [];

errTol = 1 - significanceTol;

    function s = anyIsSignificant(ind)
        for j = 1:length(ind)
            if relErr(ind(j)) <= errTol
                s = true;
                return;
            end
        end
        s = false;
    end

    function s = signifInd(ind)
        s = zeros(size(ind));
        for j = 1:length(ind)
            if relErr(ind(j)) <= errTol
                s(j) = ind(j);
            end
        end
        s(s == 0) = [];
    end
    


an = 8;
    
% Order 0
symm = 1:4;
asymm = 5:7;
solar = 9;

if(anyIsSignificant(symm))
    ptf = [ptf, signifInd(symm)];
end
if(anyIsSignificant(asymm))
    ptf = [ptf, signifInd(asymm)];
end
if(anyIsSignificant(solar))
    ptf = [ptf,solar];
end

if(~anyIsSignificant([symm,asymm]))
    ptf = removeFromPtf(solar);
end

% Order 1
offset1 = 9;
symm = (1:3)+offset1;
asymm = (4:6)+offset1;
solar = 7+offset1;
dv_amp = 8+offset1;

if(anyIsSignificant(symm))
    ptf = [ptf, signifInd(symm)];
end
if(anyIsSignificant(asymm))
    ptf = [ptf, signifInd(asymm)];
end
if(anyIsSignificant(solar))
    ptf = [ptf,solar];
end
ptf = [ptf,dv_amp];

if(~anyIsSignificant([symm,asymm]))
    ptf = removeFromPtf([solar,dv_amp]);
end

% Order 2
offset2 = 17;
symm = (1:3)+offset2;
asymm = (4:5)+offset2;
solar = 6+offset2;
dv_amp = 7+offset2;

if(anyIsSignificant(symm))
    ptf = [ptf, signifInd(symm)];
end
if(anyIsSignificant(asymm))
    ptf = [ptf, signifInd(asymm)];
end
if(anyIsSignificant(solar))
    ptf = [ptf,solar];
end
ptf = [ptf,dv_amp];

if(~anyIsSignificant([symm,asymm]))
    ptf = removeFromPtf([solar,dv_amp]);
end

% REST

asymm = [5:7, (4:6)+offset1, (4:5)+offset2];
if(anyIsSignificant(asymm))
    ptf = [ptf, an];
end


    function ptf_fixed = removeFromPtf(ind)
        ptf_fixed = ptf;
        ptf_fixed(ismember(ptf,ind)) = [];
    end

paramsToFitShort = [paramsToFitShort, ptf + offset*numStorm];

paramsToFitShort = unique(paramsToFitShort);
end