function [rho, outputStruct] =...
    il_take_time_into_account(states, S, index)

assWindow = 0.5*(S.assTimes(2)-S.assTimes(1));
if any(S.assTimes(1)-assWindow/2 < S.timestamps | S.timestamps > S.assTimes(end)+assWindow/2)
    warning('Some timestamps out of ensemble time range.')
end

rho = zeros(size(S.timestamps));
[~,beginInd] = find(min(abs(S.assTimes - S.timestamps(1))));
[~,endInd] = find(min(abs(S.assTimes - S.timestamps(end))));

for k = beginInd:endInd
    state = states(:,1,k);
    conserveInd = S.assTimes(k)-assWindow/2 <= S.timestamps & S.timestamps < S.assTimes(k)+assWindow/2;
    removeInd = ~conserveInd;
    
    S_this = removeDataPoints(S, removeInd);
    S_this = computeVariablesForFit(S_this);
    
    if state(1) > 0
        S_this.F = state(1);
        S_this.FA = state(1);
    end

    T0 = clamp(200, evalT0(S_this, S_this.T0Coeff), 1000);
    dT0 = clamp(1, evalDT(S_this, S_this.dTCoeff), 30);
    Tex = clamp(T0+1, evalTex(S_this, S_this.TexCoeff), 5000);

    % Add scalar corrections to dT0 and T0
    T0 = clamp(300, T0 + state(2), 700);
    dT0 = clamp(1, dT0 + state(3), 20);

    % Add T2 correction to Tex
    Tex = clamp(T0+1, Tex + bulge(state(4:11), S_this), 5000);

    OlbDens = evalMajorSpecies(S_this, S_this.OCoeff, S_this.O_numBiases);
    N2lbDens = evalMajorSpecies(S_this, S_this.N2Coeff, S_this.N2_numBiases);
    HelbDens = evalMajorSpecies(S_this, S_this.HeCoeff, S_this.He_numBiases);
    ArlbDens = evalMajorSpecies(S_this, S_this.ArCoeff, S_this.Ar_numBiases);
    O2lbDens = exp(S_this.O2Coeff);

    Z = computeGeopotentialHeight(S_this.altitude);

    [rho(conserveInd), O, N2, He, Ar, O2, T] = ...
        computeRho(T0, dT0, Tex, Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens);
    
end

zeroInd = rho==0;
rho(zeroInd) = interp1(1:length(rho), rho, find(zeroInd));
rho = log(rho);

if nargout > 1
    outputStruct = struct('O', O, 'N2', N2, 'He', He, 'Ar', Ar, 'O2', O2, 'Tex', Tex,...
                            'dT', dT0, 'T', T, 'T0', T0);
    error('Not implemented!')
end

end