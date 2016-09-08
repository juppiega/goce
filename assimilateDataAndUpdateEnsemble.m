function ensemble = assimilateDataAndUpdateEnsemble(ensemble, observationOperator, ...
                                                    observationsStruct)
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
%     IMPORTANT NOTE: ONLY ASSIMILATE ONE TYPE OF DATA AT A TIME! THE DATA
%     IN THE STRUCT MUST HAVE THE SAME UNIT (e.g. kg/m^-3 or K etc.).

% Argument checks.
if isrow(observationsStruct.data) || isrow(observationsStruct.sigma)
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
    Hx(:,i) = observationOperator(ensemble(:,i), observationsStruct);
end

% Simulated data perturbations from ensemble mean.
HA = bsxfun(@minus, Hx, mean(Hx,2));

% Sum of true and simulated observation covariances.
P = (1./(numMembers-1)) * HA*(HA') + diag(observationsStruct.sigma.^2);

% CREATE ANALYSIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ensemble = ensemble + (1./(numMembers-1)) * A * (HA') * (P \ (D - Hx));%!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                
end