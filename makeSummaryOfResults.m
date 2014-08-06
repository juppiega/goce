function makeSummaryOfResults(results)
%

equinoxesAndSolstices = [-0.02 0.22 0.47 0.73 0.98];
decimalYears = cell2mat(results(2:end, 3));
springStorms = [];
summerStorms = [];
autumnStorms = [];
winterStorms = [];
for i = 1:length(decimalYears)
    [~, season] = min(abs(equinoxesAndSolstices - decimalYears(i)));
    
    switch season
        case 2
            springStorms = [springStorms; i + 1];
        case 3
            summerStorms = [summerStorms; i + 1];
        case 4
            autumnStorms = [autumnStorms; i + 1];
        otherwise
            winterStorms = [winterStorms; i + 1];
    end
end

fprintf('\n%s\n', 'Seasonal distribution:')
fprintf('%s %d\n', 'Spring: ', length(springStorms));
fprintf('%s %d\n', 'Summer: ', length(summerStorms));
fprintf('%s %d\n', 'Autumn: ', length(autumnStorms));
fprintf('%s %d\n\n', 'Winter: ', length(winterStorms));

% nanIndices = cellfun(@isnan, results(2:end, 3:end), 'uniformoutput', 0);
% [row,col] = find(cell2mat(nanIndices));
% means = cellfun(@nanmean, results(2:end, row+2));
% results(nanIndices) = means;
% otherNumericVariables = cellfun(@double, results(2:end, 3:end));
% numericVariablesMean = mean(otherNumericVariables);
% numericVariablesStd = std(otherNumericVariables);
% [rows, ~] = size(results);
% results(rows + 1, 3:end) = num2cell(numericVariablesMean);
% results(rows + 2, 3:end) = num2cell(numericVariablesStd);
% results(rows + 1, 1:2) = {'Mean', 'NaN'};
% results(rows + 2, 1:2) = {'Std', 'NaN'};

smallVariationsColumns = strfind(results(1,:),'SH/NH');
smallVariationsColumns = find(~cellfun(@isempty,smallVariationsColumns));
morningVariations = cell2mat(results(2:end,smallVariationsColumns(1)));
eveningVariations = cell2mat(results(2:end,smallVariationsColumns(2)));

figure;
plot(decimalYears, eveningVariations, 'r.', 'MarkerSize', 10);
xlim([0 1]);
title('South std / North std of small variations')
ylabel('South std / North std');
xlabel('Decimal Year');
legend('Morning', 'Evening')
grid on

cell2csv('goceResults.csv', results);

end