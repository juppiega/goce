function S = removeObviousOutliers(S)

rmTol = 10;
[~,~,O,N2,He,Ar,O2,T] = computeMsis(S); 
if strcmpi(S.name,'O')
    obsToModeled = S.data ./ O;
elseif strcmpi(S.name,'N2')
    obsToModeled = S.data ./ N2;
elseif strcmpi(S.name,'He')
    obsToModeled = S.data ./ He;
elseif strcmpi(S.name,'Ar')
    obsToModeled = S.data ./ Ar;
elseif strcmpi(S.name,'O2')
    obsToModeled = S.data ./ O2;
elseif strcmpi(S.name,'T')
    obsToModeled = S.data ./ T;
end

rmInd = obsToModeled > rmTol | obsToModeled < 1/rmTol;
S = removeDataPoints(S,rmInd,true,true,true,true);

end