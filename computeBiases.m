function [] = computeBiases (name, ind)

if strcmpi(name,'O')
    load('ilData.mat','OStruct');
    S = OStruct;
elseif strcmpi(name,'N2')
    load('ilData.mat','N2Struct');
    S = N2Struct;
elseif strcmpi(name,'He')
    load('ilData.mat','HeStruct');
    S = HeStruct;
elseif strcmpi(name,'Tex')
    load('ilData.mat','TexStruct');
    S = TexStruct;
end

rmInd = setdiff(1:length(S.data),ind);
S = removeDataPoints(S, rmInd,true,true,true,true);

if strcmpi(name,'O')
    [~,~,obs_dtm] = computeDtm(S);
    [~,~,obs_msis] = computeMsis(S);
elseif strcmpi(name,'N2')
    [~,~,~,obs_dtm] = computeDtm(S);
    [~,~,~,obs_msis] = computeMsis(S);
elseif strcmpi(name,'He')
    [~,~,~,~,obs_dtm] = computeDtm(S);
    [~,~,~,~,obs_msis] = computeMsis(S);
elseif strcmpi(name,'Tex')
    [obs_dtm] = computeDtm(S);
    [obs_msis] = computeMsis(S);
end

OM_dtm = mean(obs_dtm./S.data)
OM_msis = mean(obs_msis./S.data)

end
