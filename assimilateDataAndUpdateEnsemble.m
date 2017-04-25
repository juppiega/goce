function [ensemble,innovations,residuals,covariance,P,obsRank] = ...
                    assimilateDataAndUpdateEnsemble(ensemble, observationOperator, ...
                                                    observationsStruct, statistics, computeObsRank)
% INPUT: 
%     ensemble: Ensemble of models. Double matrix size nxN, where n is the number of
%               states and N the number of members (size of the ensemble).
%     
%     observationOperator: function handle of the form f(modelState,
%     observationsStruct), where modelState is one ensemble member and
%     
%     observationsStruct: Struct of the (spatio-temporal) location and value
%                         of observations. Must contain fields "sigma"
%                         describing the 1-sigma uncertainty (not variance) and "data",
%                         which contains the values. "data" and "sigma"
%                         must be in absolute (not relative) units (e.g. Kelvin).
%     statistics:         Compute additional statistics (boolean).
%     IMPORTANT NOTE: ONLY ASSIMILATE ONE TYPE OF DATA AT A TIME! THE DATA
%     IN THE STRUCT MUST HAVE THE SAME UNIT (e.g. kg/m^-3 or K etc.).

% Argument checks.
if length(observationsStruct.data) > 1 &&...
    (isrow(observationsStruct.data) || isrow(observationsStruct.sigma))
    error('Data in the observation struct must be aligned in column vectors!')
end
if length(observationsStruct.data) == 0
    warning('No data assimilated. Need at least one measurement!')
    return;
end

% Constants
numMembers    = size(ensemble, 2);
numVariables  = size(ensemble, 1);
numData       = length(observationsStruct.data);

% Perturb data with N(0,obs.sigma)
dataPerturbations = bsxfun(@times, observationsStruct.sigma, randn(numData, numMembers));
% Make perturbation ensemble mean (rows) unbiased.
dataPerturbations = bsxfun(@minus, dataPerturbations, mean(dataPerturbations,2));

%Perturbed data.
D = repmat(observationsStruct.data, 1, numMembers) + dataPerturbations;

% State perturbations from the ensemble mean.
A = bsxfun(@minus, ensemble, mean(ensemble, 2));

% Simulate observations using the background ensemble members.
Hx = zeros(numData, numMembers);
for i = 1:numMembers
    simObs = observationOperator(ensemble(:,i), observationsStruct, i);
    if any(isnan(simObs) | ~isfinite(simObs) | ~isreal(simObs))
        error('Observation operator returned unphysical value!')
    end
    Hx(:,i) = simObs;
end

if statistics || (nargin>4 && computeObsRank)
    [~,rankedEnsembleAndObs] = sort([Hx, observationsStruct.data], 2);
    [~,obsRank] = find(rankedEnsembleAndObs == numMembers+1);
    innovations = [];
    residuals = [];
    covariance = [];
end

% Simulated data perturbations from ensemble mean.
HA = bsxfun(@minus, Hx, mean(Hx,2));

% Sum of true and simulated observation covariances.
P = (1./(numMembers-1)) * HA*(HA') + diag(observationsStruct.sigma.^2);

% CREATE ANALYSIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ensemble = ensemble + (1./(numMembers-1)) * A * (HA') * (P \ (D - Hx));%!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% ensemble(1,:) = clamp(500, ensemble(1,:), 4000);
% ensemble(2,:) = clamp(log(1E8), ensemble(2,:), log(1E11));
% ensemble(3,:) = clamp(log(1E9), ensemble(2,:), log(1E12));

if statistics
    innovations = mean(bsxfun(@minus, observationsStruct.data, Hx),2);

    residuals = zeros(numData, 1);
    for i = 1:numMembers
        simObs = observationOperator(ensemble(:,i), observationsStruct);
        if any(isnan(simObs) | ~isfinite(simObs) | ~isreal(simObs))
            error('Observation operator returned unphysical value!')
        end
        residuals = residuals + (observationsStruct.data - simObs)/numMembers;
    end

    % Model covariance.
    covariance = (1./(numMembers-1)) * A * (A');
end
                                                
end