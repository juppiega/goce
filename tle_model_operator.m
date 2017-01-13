function [rhoModel] =...
    tle_model_operator(state, obsStruct, index)

rhoModel = obsStruct.rhoModel(:,index);

end