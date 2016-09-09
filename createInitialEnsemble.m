function ensemble = createInitialEnsemble(modelString, numMembers, initialGuess, parameter1sigmaError)

load optCoeff.individualAE.mat

if strcmpi(modelString,'dummy')
    rhoGuess = [optCoeff(OInd(1)), optCoeff(N2Ind(1)), optCoeff(HeInd(1))];
    initialGuess =  [507, 12, 1000, rhoGuess];
    stateLength = length(initialGuess);
    lb = [400, 7, 700, rhoGuess*0.99];
    ub = [600, 17, 1400, rhoGuess*1.01];
    ensemble = zeros(stateLength, numMembers);
    for i = 1:stateLength
        d = ub(i)-lb(i);
        m = 0.5*(ub(i)+lb(i));
        ensemble(i,:) = (rand(1, numMembers)-0.5)*d + m;
    end
else
    error('Unidentified model string!')
end

end