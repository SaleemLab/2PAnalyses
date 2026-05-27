function plotExampleROIsBySMI(sessionMatrix)
    % 1. Pool SMI and map back to sessions
    allSMI = [];
    roiLookup = []; % [sessionID, roiID]
    
    for s = 1:length(sessionMatrix)
        base = sessionMatrix(s).ConditionData.Baseline;
        if ~isfield(base, 'SMI'), continue; end
        
        valid = find(~isnan(base.SMI));
        for i = 1:length(valid)
            allSMI = [allSMI; base.SMI(valid(i))];
            roiLookup = [roiLookup; s, valid(i)];
        end
    end
    
    % 2. Identify Percentile ROIs
    p = prctile(allSMI, [25, 50, 75]);
    reps = zeros(3,1);
    for k = 1:3
        [~, reps(k)] = min(abs(allSMI - p(k)));
    end
    
    % 3. Plotting
    figure('Color', 'w', 'Position', [50 100 1300 700]);
    titles = {'25th Percentile', '50th Percentile', '75th Percentile'};
    
    for k = 1:3
        sID = roiLookup(reps(k), 1);
        rID = roiLookup(reps(k), 2);
        
        % Get data
        trialData = squeeze(sessionMatrix(sID).ConditionData.Baseline.LapActivity(rID, :, :));
        meanCurve = mean(trialData, 1, 'omitnan');
        
        % --- Subplot TOP: Mean Tuning Curve ---
        subplot(3, 3, k); 
        plot(meanCurve, 'k', 'LineWidth', 2);
        title(sprintf('%s (SMI: %.2f)', titles{k}, p(k)));
        ylabel('Mean dF/F');
        set(gca, 'XTick', [40 80 120 160], 'XTickLabel', {}, 'XGrid', 'on');
        xlim([1 size(trialData,2)]);
        
        % --- Subplot BOTTOM: Trial Heatmap ---
        subplot(3, 3, [k+3, k+6]); % Spans the bottom two rows
        imagesc(trialData);
        colormap(gca, 'parula');
        xlabel('Position (cm)');
        if k == 1, ylabel('Laps'); end
        
        set(gca, 'XTick', [40 80 120 160]);
        grid on; set(gca, 'GridColor', 'w', 'GridAlpha', 0.2);
        
        % Scale heatmap for visibility
        caxis([0 prctile(trialData(:), 98) + 1e-9]); 
        
        % Link the X-axes for zooming
        if k > 1, linkaxes(findobj(gcf,'Type','axes'),'x'); end
    end
end