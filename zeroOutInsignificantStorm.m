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
san = 10;
    
% Order 0
symm = 1:4;
asymm = 5:7;
san_amp = 9;
solar = 11;

if(anyIsSignificant(symm))
    ptf = [ptf, signifInd(symm)];
end
if(anyIsSignificant(asymm))
    ptf = [ptf, signifInd(asymm)];
end
if(anyIsSignificant(san_amp))
    ptf = [ptf, san_amp];
end
if(anyIsSignificant(solar))
    ptf = [ptf,solar];
end

if(~anyIsSignificant(asymm))
    ptf = removeFromPtf(san_amp);
end
if(~anyIsSignificant([symm,asymm]))
    ptf = removeFromPtf(solar);
end

% Order 1
symm = (1:3)+11;
asymm = (4:6)+11;
san_amp = 7+11;
solar = 8+11;
dv_amp = 9+11;

if(anyIsSignificant(symm))
    ptf = [ptf, signifInd(symm)];
end
if(anyIsSignificant(asymm))
    ptf = [ptf, signifInd(asymm)];
end
if(anyIsSignificant(san_amp))
    ptf = [ptf, san_amp];
end
if(anyIsSignificant(solar))
    ptf = [ptf,solar];
end
if(anyIsSignificant(dv_amp))
    ptf = [ptf,dv_amp];
end

if(~anyIsSignificant(asymm))
    ptf = removeFromPtf(san_amp);
end
if(~anyIsSignificant([symm,asymm]))
    ptf = removeFromPtf([solar,dv_amp]);
end

% Order 2
symm = (1:3)+20;
asymm = (4:5)+20;
san_amp = 6+20;
solar = 7+20;
dv_amp = 8+20;

if(anyIsSignificant(symm))
    ptf = [ptf, signifInd(symm)];
end
if(anyIsSignificant(asymm))
    ptf = [ptf, signifInd(asymm)];
end
if(anyIsSignificant(san_amp))
    ptf = [ptf, san_amp];
end
if(anyIsSignificant(solar))
    ptf = [ptf,solar];
end
if(anyIsSignificant(dv_amp))
    ptf = [ptf,dv_amp];
end

if(~anyIsSignificant(asymm))
    ptf = removeFromPtf(san_amp);
end
if(~anyIsSignificant([symm,asymm]))
    ptf = removeFromPtf([solar,dv_amp]);
end

% For an or san all zero
asymm = [5:7, (4:6)+11, (4:5)+20];
san_amp = [9, 7+11, 6+20];
if(~anyIsSignificant(asymm))
    ptf = removeFromPtf([an,san]);
end
if(~anyIsSignificant(san_amp))
    ptf = removeFromPtf(san);
end


    function ptf_fixed = removeFromPtf(ind)
        ptf_fixed = ptf;
        ptf_fixed(ismember(ptf,ind)) = [];
    end

paramsToFitShort = [paramsToFitShort, ptf + offset*numStorm];

paramsToFitShort = unique(paramsToFitShort);
end