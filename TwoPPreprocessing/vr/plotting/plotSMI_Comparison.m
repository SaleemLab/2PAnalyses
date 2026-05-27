function plotSMI_Comparison(sessionMatrix)
    % plotSMI_Comparison: Overlays Averaged SMI vs Single-Trial SMI
    
    allMeanSMI = [];
    allSingleSMI = [];
    
    % 1. Pool data from all sessions
    for s = 1:length(sessionMatrix)
        base = sessionMatrix(s).ConditionData.Baseline;
        
        % Pool Mean-Averaged values
        if isfield(base, 'SMI')
            vals = base.SMI(:);
            allMeanSMI = [allMeanSMI; vals(~isnan(vals))];
        end
        
        % Pool Single-Trial values (flattening the matrix)
        if isfield(base, 'SingleTrialSMI')
            vals = base.SingleTrialSMI(:);
            allSingleSMI = [allSingleSMI; vals(~isnan(vals))];
        end
    end
    
    % 2. Create Plot
    figure('Color', 'w', 'Position', [100 100 550 500]);
    hold on;
    
    % Plot Single-Trial (The "Noisy" distribution)
    [f_single, x_single] = ecdf(allSingleSMI);
    plot(x_single, f_single, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'LineStyle', '--');
    
    % Plot Mean-Averaged (The "Clean" distribution)
    [f_mean, x_mean] = ecdf(allMeanSMI);
    plot(x_mean, f_mean, 'k', 'LineWidth', 2.5);
    
    % 3. Reference Lines matching the paper
    line([0 0], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); % Zero line
    line([-1 1], [0.5 0.5], 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); % Median guide
    
    % 4. Aesthetics
    xlabel('Spatial Modulation Index (SMI)');
    ylabel('Cumulative Probability');
    title('V1 Spatial Modulation: Averaging vs. Trials');
    xlim([-1 1]); ylim([0 1]);
    grid on;
    
    legend({sprintf('Single-Trial (n=%d)', length(allSingleSMI)), ...
            sprintf('Mean-Averaged (n=%d)', length(allMeanSMI))}, ...
            'Location', 'southeast', 'Box', 'off');
            
    % Diagnostic Printout for your supervisor
    fprintf('--- Comparison Results ---\n');
    fprintf('Mean-Averaged Median: %.3f\n', median(allMeanSMI));
    fprintf('Single-Trial Median:  %.3f\n', median(allSingleSMI));
    fprintf('%% Values < 0 (Averaged): %.1f%%\n', 100*mean(allMeanSMI < 0));
    fprintf('%% Values < 0 (Single):   %.1f%%\n', 100*mean(allSingleSMI < 0));
end