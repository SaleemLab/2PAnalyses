function fig = plotConditionTuning(allData, targetArea, days, savePath)
    % plotConditionTuning: Individual plot normalization, fixed Baseline-Odd sorting.
    
    % 1. Identify manipulation conditions
    allConds = {};
    for i = 1:length(allData)
        if isfield(allData(i), 'ConditionData') && ~isempty(allData(i).ConditionData)
            names = fieldnames(allData(i).ConditionData);
            manips = names(~contains(lower(names), 'baseline') & ~contains(lower(names), 'norm'));
            allConds = [allConds, manips'];
        end
    end
    manipConds = unique(allConds, 'stable');
    plotLabels = [{'Baseline Odd', 'Baseline Even'}, manipConds];
    numCols = length(plotLabels);
    numDays = length(days);

    % Colors
    omitColor = [0.9 0.2 0.2]; % Red
    swapColor = [0.2 0.4 0.9]; % Blue
    baseColor = [0 0 0];       % Black
    
    fig = figure('Color', 'w', 'Position', [50 100 300*numCols + 120 400*numDays]);
    t = tiledlayout(numDays, numCols, 'TileSpacing', 'none', 'Padding', 'compact');
    
    for d = 1:numDays
        currentDay = days(d);
        targetDayIdx = [allData.Day] == currentDay;
        
        areaIdx = strcmpi(string({allData.TargetArea_ROI}), string(targetArea));
        sessList = allData(targetDayIdx & areaIdx);
        if isempty(sessList), continue; end
        
        dayData = cell(1, numCols);
        for s = 1:length(sessList)
            cData = sessList(s).ConditionData;
            cNames = fieldnames(cData);
            baseIdx = find(contains(lower(cNames), 'baseline') | contains(lower(cNames), 'norm'), 1);
            if isempty(baseIdx), continue; end
            baseName = cNames{baseIdx};
            
            % ROI selection
            if isfield(sessList(s), 'highlyCorrBoutons') && ~isempty(sessList(s).highlyCorrBoutons)
                idx = sessList(s).highlyCorrBoutons;
            else
                idx = 1:sessList(s).NumCells;
            end
            
            % Baseline Odd/Even
            fullBase = cData.(baseName).LapActivity(idx, :, :);
            dayData{1} = [dayData{1}; squeeze(mean(fullBase(:, 1:2:end, :), 2, 'omitnan'))];
            dayData{2} = [dayData{2}; squeeze(mean(fullBase(:, 2:2:end, :), 2, 'omitnan'))];
            
            % Manipulations
            for m = 1:length(manipConds)
                mName = manipConds{m};
                if isfield(cData, mName)
                    dayData{m+2} = [dayData{m+2}; squeeze(mean(cData.(mName).LapActivity(idx, :, :), 2, 'omitnan'))];
                end
            end
        end
        
        % --- SORTING ANCHOR (Baseline Odd) ---
        if isempty(dayData{1}), continue; end
        baseOdd = dayData{1};
        % Self-normalize just the reference for sorting purposes
        refNorm = (baseOdd - min(baseOdd,[],2,'omitnan')) ./ (max(baseOdd,[],2,'omitnan') - min(baseOdd,[],2,'omitnan'));
        [~, peakIdx] = max(refNorm, [], 2);
        [~, sortIdx] = sort(peakIdx);
        
        % --- PLOTTING ---
        for col = 1:numCols
            ax = nexttile();
            currRawData = dayData{col};
            currLabel = lower(plotLabels{col});
            
            if contains(currLabel, 'omit')
                activeColor = omitColor;
            elseif contains(currLabel, 'swap')
                activeColor = swapColor;
            else
                activeColor = baseColor;
            end
            
            if ~isempty(currRawData)
                % SELF-NORMALIZATION: Each plot scaled to its own 0-1 range
                cMin = min(currRawData, [], 2, 'omitnan');
                cMax = max(currRawData, [], 2, 'omitnan');
                cRange = cMax - cMin;
                cRange(cRange == 0) = 1;
                normData = (currRawData - cMin) ./ cRange;
                
                imagesc(normData(sortIdx, :));
                colormap(ax, flipud(gray)); clim([0 1]);
                set(ax, 'YDir', 'normal', 'TickDir', 'out', 'FontSize', 10);
                
                % Title
                if d == 1
                    tH = title(plotLabels{col}, 'Interpreter', 'none');
                    tH.Color = activeColor;
                end
                
                % Landmark Ticks and Lines
                landmarks = [40, 80, 120, 160];
                xticks(landmarks);
                xticklabels({'40', '80', '120', '160'});
                
                for l = 1:4
                    isHighlighted = contains(currLabel, num2str(l));
                    lColor = [0 0 0]; lWidth = 0.5;
                    if isHighlighted
                        lColor = activeColor; lWidth = 1.8;
                    end
                    xline(landmarks(l), '--', 'Color', lColor, 'LineWidth', lWidth);
                end
            end
            
            if col == 1
                ylabel(sprintf('Day %d\n(n=%d Boutons)', currentDay, size(dayData{1},1)), 'FontWeight', 'bold');
            else
                set(ax, 'YTickLabel', []);
            end
            
            if col == numCols
                cb = colorbar;
                ylabel(cb, '\DeltaF/F (norm.)');
            end
        end
    end
    
    xlabel(t, 'Position (cm)', 'FontSize', 12, 'FontWeight', 'bold');
    
    if nargin > 3 && ~isempty(savePath)
        exportgraphics(fig, savePath, 'Resolution', 300);
    end
end