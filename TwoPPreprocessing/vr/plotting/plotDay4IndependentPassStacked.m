function plotDay4IndependentPassStacked(allData)
    saveDir = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\DistributionsAllCriteria';
    if ~exist(saveDir, 'dir'), mkdir(saveDir); end
    
    % Thresholds
    T_CORR = 0.4;
    T_VAR_RANGE = 1;
    T_VAR_VAR = 20;

    % Filter for Day 4
    if isnumeric(allData(1).Day)
        dayIdx = [allData.Day] == 5;
    else
        dayIdx = str2double({allData.Day}) == 5;
    end
    day4Data = allData(dayIdx);

    mice = unique({day4Data.MouseID});
    numMice = length(mice);
    
    % Categories: Total followed by specific independent PASSING counts
    filterLabels = {'Total ROIs', 'Peak Shuff', 'Range Shuff', ...
                    'Halves Corr', 'Variance/TuningVariance Ratio', 'Variance/TuningRange Ratio', 'Pass ALL'};
    
    % Rows = Mice, Columns = 7 categories
    mouseCounts = zeros(numMice, 7);

    for i = 1:numMice
        mouseID = mice{i};
        mouseIdx = strcmp({day4Data.MouseID}, mouseID);
        mouseSessions = day4Data(mouseIdx);
        
        for s = 1:length(mouseSessions)
            sess = mouseSessions(s);
            
            % 1. Total
            nTotal = sess.NumCells;
            
            % 2. Pass Peak Shuffling
            isPeak = sess.isSignificantByPeakShuffling == 1;
            
            % 3. Pass Range Shuffling 
            isRangeShuff = sess.isSignificantByRange == 1;
            
            % 4. Pass Stability (Correlation >= threshold)
            isStable = sess.lapCorr_HalvesRho >= T_CORR;
            
            % 5. Pass Range Ratio (Ratio <= threshold)
            isRangeRatio = sess.ratioVarToTuningRange <= T_VAR_RANGE;
            
            % 6. Pass Var Ratio (Ratio <= threshold)
            isVarRatio = sess.ratioVarToTuningVar <= T_VAR_VAR;
            
            % 7. Pass ALL (The actual survivors)
            isAll = isPeak & isRangeShuff & isStable & isRangeRatio & isVarRatio;
            
            % Aggregate counts for the mouse
            mouseCounts(i, :) = mouseCounts(i, :) + ...
                [nTotal, sum(isPeak), sum(isRangeShuff), sum(isStable), ...
                 sum(isRangeRatio), sum(isVarRatio), sum(isAll)];
        end
    end

    % --- Plotting ---
    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.7 0.6]);
    
    % Transpose for stacked bar: Mice become segments of each bar
    b = bar(mouseCounts', 'stacked', 'FaceColor', 'flat');
    
    % Apply mouse-specific colors
    colors = lines(numMice);
    for k = 1:numMice
        b(k).CData = colors(k, :);
        b(k).DisplayName = sprintf('Mouse %s', mice{k});
    end

    % Formatting
    set(gca, 'XTickLabel', filterLabels, 'TickLabelInterpreter', 'none', 'FontSize', 10);
    xtickangle(30);
    ylabel('Number of Boutons (Pooled)', 'FontSize', 12);
    title('Day 5 ROI Yield', 'FontSize', 14);
    
    
    set(gca, 'Layer', 'top', 'GridAlpha', 0.3); 
    legend('Location', 'northeastoutside');

    % Save
    saveName = fullfile(saveDir, 'Day4_Independent_Pass_Stacked.png');
    exportgraphics(fig, saveName, 'Resolution', 300);
    fprintf('Pass-rate stacked plot saved to: %s\n', saveName);
end