function plotPooledSMI_CDF(pooledData, targetArea, compareBy)
    % plotPooledSMI_CDF: Plots CDF of SMI values for a given area.
    % compareBy: 'Total' (one curve) or 'Days' (separate curve per day).
    
    if ~isfield(pooledData, targetArea)
        error('Area %s not found in pooled data.', targetArea);
    end
    
    figure('Color', 'w', 'Position', [200 200 500 450]);
    hold on;
    
    if strcmpi(compareBy, 'Total')
        % Plot single curve for the entire area
        [f, x] = ecdf(pooledData.(targetArea).AllSMI);
        plot(x, f, 'k', 'LineWidth', 2.5);
        legendLabels = {sprintf('%s (n=%d)', targetArea, length(pooledData.(targetArea).AllSMI))};
        
    elseif strcmpi(compareBy, 'Days')
        % Plot separate curve for each Day of Experience
        dayFields = fieldnames(pooledData.(targetArea).Days);
        colors = lines(length(dayFields));
        legendLabels = {};
        
        for d = 1:length(dayFields)
            daySMI = pooledData.(targetArea).Days.(dayFields{d});
            if isempty(daySMI), continue; end
            
            [f, x] = ecdf(daySMI);
            plot(x, f, 'Color', colors(d,:), 'LineWidth', 2);
            legendLabels{end+1} = sprintf('%s: %s (n=%d)', targetArea, dayFields{d}, length(daySMI));
        end
    end
    
    % Aesthetics matching Diamanti et al.
    line([0 0], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); % Zero line
    line([-1 1], [0.5 0.5], 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); % Median guide
    
    grid on;
    xlim([-1 1]); ylim([0 1]);
    xlabel('Spatial Modulation Index (SMI)');
    ylabel('Cumulative Probability');
    title(sprintf('Spatial Modulation in %s', targetArea));
    legend(legendLabels, 'Location', 'southeast', 'Box', 'off');
end