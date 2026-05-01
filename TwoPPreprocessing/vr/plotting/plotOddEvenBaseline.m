function fig = plotOddEvenBaseline(allData, targetArea, days, savePath)
    % Pooled Odd/Even Plotting - No-Fail Version
    
    pooledOdd = [];
    pooledEven = [];
    
    % Force targetArea to string
    targetAreaStr = string(targetArea);
    
    for i = 1:length(allData)
        % Check Area (Case-insensitive string check)
        if isfield(allData(i), 'TargetArea_ROI')
            currentArea = string(allData(i).TargetArea_ROI);
        else
            currentArea = string(allData(i).TargetArea);
        end
        
        % If Area doesn't match, skip session
        if ~strcmpi(currentArea, targetAreaStr), continue; end
        
        % Check Day (Matches if any provided day matches current day)
        % We use string comparison here to handle "200" vs 200
        if ~any(string(allData(i).Day) == string(days)), continue; end
        
        % Access ConditionData
        if ~isfield(allData(i), 'ConditionData') || isempty(allData(i).ConditionData)
            continue;
        end
        
        cData = allData(i).ConditionData;
        cNames = fieldnames(cData);
        
        % Find Baseline: look for 'base', 'norm', or just take the first field
        idxB = find(contains(lower(cNames), 'base') | contains(lower(cNames), 'norm'), 1);
        if isempty(idxB), idxB = 1; end
        baseName = cNames{idxB};
        
        % Get ROI Indices (Use highlyCorr if it exists, otherwise all)
        if isfield(allData(i), 'highlyCorrBoutons') && ~isempty(allData(i).highlyCorrBoutons)
            roiIdx = allData(i).highlyCorrBoutons;
        else
            roiIdx = 1:allData(i).NumCells;
        end
        
        % Pull activity and calculate means
        try
            baseAct = cData.(baseName).LapActivity(roiIdx, :, :);
            
            % Odd/Even Split
            oMean = squeeze(mean(baseAct(:, 1:2:end, :), 2, 'omitnan'));
            eMean = squeeze(mean(baseAct(:, 2:2:end, :), 2, 'omitnan'));
            
            pooledOdd  = [pooledOdd;  oMean];
            pooledEven = [pooledEven; eMean];
        catch
            continue; % Skip if LapActivity is missing or wrong shape
        end
    end
    
    % Final check
    if isempty(pooledOdd)
        error('STILL NO DATA: Check if allData(%d).ConditionData.%s.LapActivity exists.', i, baseName);
    end
    
    % --- Normalization & Sorting ---
    oMin = min(pooledOdd, [], 2, 'omitnan');
    oMax = max(pooledOdd, [], 2, 'omitnan');
    rangeVal = oMax - oMin;
    rangeVal(rangeVal == 0) = 1;
    
    oddNorm  = (pooledOdd - oMin) ./ rangeVal;
    evenNorm = (pooledEven - oMin) ./ rangeVal;
    
    [~, peakPos] = max(oddNorm, [], 2);
    [~, sortIdx] = sort(peakPos);
    
    % --- Plotting ---
    fig = figure('Color', 'w', 'Position', [100 100 900 450]);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
    
    data = {oddNorm(sortIdx, :), evenNorm(sortIdx, :)};
    lbls = {'Odd Laps', 'Even Laps'};

    for k = 1:2
        ax = nexttile();
        imagesc(data{k});
        colormap(ax, flipud(gray)); clim([0 1]);

        title(lbls{k}, 'FontSize', 14);
        set(ax, 'YDir', 'normal', 'TickDir', 'out', 'FontSize', 12);
        xline([40 80 120 160], 'k--', 'LineWidth', 1.5);
        xticks([1 40 80 120 160 200]);
        xticklabels({'1','40','80','120','160','200'});

        if k == 1
            ylabel('Somas', 'FontWeight', 'bold');
        else
            set(ax, 'YTick', []);
            cb = colorbar;
            ylabel(cb, '\DeltaF/F (norm.)', 'FontSize', 12);
            set(cb, 'Ticks', [0 0.5 1], 'TickLabels', {'0', '0.5', '1'});
        end
    end

    xlabel(t, 'Position (cm)', 'FontSize', 14, 'FontWeight', 'bold');
    
    if nargin > 3 && ~isempty(savePath)
        exportgraphics(fig, savePath, 'Resolution', 300);
    end
end