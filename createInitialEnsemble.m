function ensemble = createInitialEnsemble(modelString, numMembers, initialGuess, parameter1sigmaError)

load optCoeff.individualAE.mat

if strcmpi(modelString,'dummy') || strcmpi(modelString,'full')
    initialGuess = zeros(11, 1);
    stateLength = length(initialGuess);
    lb = [50, -100, -5, zeros(1,8)-50];
    ub = [150, 100, 5, zeros(1,8)+50];
    %lb = [50, zeros(1,10)];
    %ub = [150, zeros(1,10)];
    ensemble = zeros(stateLength, numMembers);
    for i = 1:stateLength
        d = (ub(i)-lb(i))/2;
        m = 0.5*(ub(i)+lb(i));
%         ensemble(i,:) = (rand(1, numMembers)-0.5)*d + m;
        ensemble(i,:) = randn(1, numMembers)*d + m;
    end
else
    error('Unidentified model string!')
end

end