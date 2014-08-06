function results = plotCorrelation(xvals, yvals, xvalName, yvalName, plotFigures, results)
%

r = corr(xvals, yvals);

[rowNum, ~] = size(results);
if rowNum > 0
    emptyCells = cellfun(@isempty,results);
    [~, emptyColPositions] = find(emptyCells);
    colNum = min(emptyColPositions);
    if rowNum == 2
        colNum = length(results(rowNum,:)) + 1;
        results{1, colNum} = ['r: ',xvalName];
    end
    results{rowNum, colNum} = r;
end

if plotFigures ~= 0
    r2 = r * r;
    fprintf('%s %f\n', 'Pearson correlation: ', r)
    fprintf('%s %d %s\n', 'Thus, ', round(r2 * 100), ['% of variation in ', yvalName, ' can be explained by changes in ', xvalName])
    
    figure;
    plot(xvals, yvals, '.')
    p = polyfit(xvals, yvals, 1);
    m = p(1);
    b = p(2);
    ylimits = get(gca, 'ylim');
    xlimits = get(gca, 'xlim');
    regLineYAxisCross = m * xlimits(1) + b;
    regLineXAxisCross = (ylimits(2) - b) / m;
    line([xlimits(1) regLineXAxisCross], [regLineYAxisCross ylimits(2)]);
    title([yvalName, ' vs ', xvalName ,' correlation'])
    xlabel(xvalName)
    ylabel(yvalName)
    
    textYLocation = ylimits(2) - 0.05 * (ylimits(2) - ylimits(1));
    textXLocation = xlimits(2) - 0.05 * (xlimits(2) - xlimits(1));
    str1 = [' r =  ', num2str(r, '%07.4f')];
    str2 = ['r^2 = ', num2str(r2, '%07.4f')];
    textString = [str1 ; str2];
    text(textXLocation, textYLocation, textString, 'FontSize', 12, ...
        'VerticalAlignment','top', 'HorizontalAlignment','right', 'EdgeColor', [0 0 0]);
end

end