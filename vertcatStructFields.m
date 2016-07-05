function S = vertcatStructFields(S)

if size(S,2) ~= 1
    error('Number of columns in the array of structures must be one!')
end

names = unique(fieldnames(S), 'stable');
data = squeeze(struct2cell(S)).';

combinedData = arrayfun(@(i) vertcat(data{:,i}), 1:size(data, 2), 'UniformOutput', false);
S = cell2struct(combinedData, names, 2);

end
