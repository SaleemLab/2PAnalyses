function plotSMI_PercentileAnalysis(sessionMatrix)
    allTuning = [];
    allSMI = [];
    allPrefBins = [];
    
    for s = 1:length(sessionMatrix)
        base = sessionMatrix(s).ConditionData.Baseline;
        if ~isfield(base, 'SMI'), continue; end
        
        data = base.LapActivity;
        % Use Mean Even for the tuning curves (independent of the peak-finding set)
        meanEven = squeeze(mean(data(:, 2:2:end, :), 2, 'omitnan'));
        meanOdd = squeeze(mean(data(:, 1:2:end, :), 2, 'omitnan'));
        
        validMask = ~isnan(base.SMI);
        
        allTuning = [allTuning; meanEven(validMask, :)];
        allSMI = [allSMI; base.SMI(validMask)];
        
        % Calculate prefBins from Odd trials (consistent with SMI calc)
        [~, prefBins] = max(meanOdd(validMask, 16:185), [], 2); 
        allPrefBins = [allPrefBins; prefBins + 15]; 
    end
    
    % 1. Calculate Percentiles
    p_vals = prctile(allSMI, [25, 50, 75]);
    
    % Find ROIs closest to these values
    [~, idx25] = min(abs(allSMI - p_vals(1)));
    [~, idx50] = min(abs(allSMI - p_vals(2)));
    [~, idx75] = min(abs(allSMI - p_vals(3)));
    repIndices = [idx25, idx50, idx75];
    
    % 2. Normalize Heatmap
    % Subtract min and divide by max per ROI to get [0 1] range
    normTuning = allTuning - min(allTuning, [], 2);
    normTuning = normTuning ./ max(normTuning, [], 2);
    
    % Sort by preferred position
    [~, sortIdx] = sort(allPrefBins);
    sortedTuning = normTuning(sortIdx, :);
    
    figure('Color', 'w', 'Position', [100 100 1000 500]);
    
    % --- SUBPLOT 1: Population Heatmap ---
    subplot(1, 2, 1);
    imagesc(sortedTuning);
    colormap(gca, 'parula');
    title(sprintf('V1 Population sorted by Peak (n=%d)', size(sortedTuning, 1)));
    xlabel('Position (cm)'); ylabel('Sorted ROIs');
    
    % Fixed Colorbar logic
    c = colorbar('southoutside');
    c.Label.String = 'Normalized Activity'; 
    
    % --- SUBPLOT 2: Percentile Examples ---
    subplot(1, 2, 2);
    hold on;
    colors = [0.6 0.6 0.6; 0.2 0.6 1; 1 0.2 0.2]; % Gray, Blue, Red
    labels = {sprintf('25th (SMI: %.2f)', p_vals(1)), ...
              sprintf('50th (SMI: %.2f)', p_vals(2)), ...
              sprintf('75th (SMI: %.2f)', p_vals(3))};
    
    for i = 1:3
        targetIdx = repIndices(i);
        % Smooth the example curve slightly for better visualization
        trace = allTuning(targetIdx, :);
        plot(trace, 'Color', colors(i,:), 'LineWidth', 2.5);
    end
    
    legend(labels, 'Location', 'northeast', 'Box', 'off');
    title('Representative Tuning Curves');
    xlabel('Position (cm)'); ylabel('Activity (dF/F)');
    xlim([1 200]); grid on;
end