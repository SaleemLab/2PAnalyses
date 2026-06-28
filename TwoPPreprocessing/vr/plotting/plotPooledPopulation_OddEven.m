function figHandle = plotPooledPopulation_OddEven(allData, targetArea, varargin)
% plotPooledPopulation_OddEven: Baseline function plotting side-by-side Odd
% and Even lap profiles. Uses 'smoothLapActivity' for trace smoothing and 
% Includes population means at the bottom for both odd and even 
    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    addRequired(p, 'targetArea', @(x) ischar(x) || isstring(x));
    
    addParameter(p, 'DaysToPlot', [1, 2, 3, 4, 5], @isnumeric);
    addParameter(p, 'TypeToPlot', 'Boutons', @(x) ischar(x) || isstring(x));
    addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x)); 
    addParameter(p, 'ApplySmoothing', true, @islogical);
    addParameter(p, 'FontName', 'Arial', @(x) ischar(x) || isstring(x)); 
    
    parse(p, allData, targetArea, varargin{:});
    
    targetFont = p.Results.FontName;
    
    if all(cellfun(@isempty, {allData.TargetArea})), [allData.TargetArea] = deal(char(targetArea)); end
    if isfield(allData, 'Type') && ~isfield(allData, 'TypeImaged')
        [allData.TypeImaged] = allData.Type;
    elseif all(cellfun(@isempty, {allData.TypeImaged}))
        [allData.TypeImaged] = deal(char(p.Results.TypeToPlot));
    end
    
    daysToPlot = p.Results.DaysToPlot;
    daysToPlot(daysToPlot == 200) = []; 
    nDays = length(daysToPlot);
    
    plotQueue = struct('matrixToPlot', {}, 'numN', {}, 'titleStr', {}, 'thisCol', {}, 'rawMatrix', {});
    titleColors = {'r', [0.4 0.7 0.2], 'b', 'm', 'k'};
    
    % --- DATA PROCESSING LOOP ---
    for d = 1:nDays
        day = daysToPlot(d);
        thisCol = titleColors{mod(d-1,5)+1};
        
        if day == 5
            dayMask = ([allData.Day] == 5 | [allData.Day] == 200);
            titleDayStr = '5 + 200'; 
        else
            dayMask = ([allData.Day] == day);
            titleDayStr = sprintf('%d', day);
        end
        
        daySessions = allData(dayMask & ...
                              strcmpi(string({allData.TargetArea}), string(targetArea)) & ...
                              strcmpi(string({allData.TypeImaged}), string(p.Results.TypeToPlot)));
        
        if isempty(daySessions), continue; end
        
        allOdd = []; allEven = [];
        for s = 1:length(daySessions)
            thisSess = daySessions(s);
            if ~isfield(thisSess, 'ConditionData') || isempty(thisSess.ConditionData), continue; end
            
            condNames = fieldnames(thisSess.ConditionData);
            if ismember('Baseline', condNames), activeCond = 'Baseline';
            elseif ismember('Default', condNames), activeCond = 'Default';
            else, activeCond = condNames{1}; end
            
            if ~isfield(thisSess.ConditionData, activeCond), continue; end
            
            lapActivity = thisSess.ConditionData.(activeCond).LapActivity;
            if isfield(thisSess, 'FilteredROIs') && ~isempty(thisSess.FilteredROIs)
                idx = thisSess.FilteredROIs;
            else
                idx = 1:size(lapActivity, 1);
            end
            if isempty(idx), continue; end
            
            lapActivity = lapActivity(idx, :, :);
            
            if p.Results.ApplySmoothing
                lapActivity = smoothLapActivity(lapActivity);
            end
            
            nTotalLaps = size(lapActivity, 2);
            allOdd  = vertcat(allOdd, squeeze(mean(lapActivity(:, 1:2:nTotalLaps, :), 2, 'omitnan')));
            allEven = vertcat(allEven, squeeze(mean(lapActivity(:, 2:2:nTotalLaps, :), 2, 'omitnan')));
        end
        
        if isempty(allEven), continue; end
        if size(allOdd, 2) == 1, allOdd = allOdd'; allEven = allEven'; end
        
        minOdd = min(allOdd, [], 2); 
        maxOdd = max(allOdd, [], 2); 
        rangeOdd = maxOdd - minOdd;
        rangeOdd(rangeOdd == 0) = 1;
        
        normOdd  = (allOdd  - minOdd) ./ rangeOdd; 
        normEven = (allEven - minOdd) ./ rangeOdd;
        normOdd(isnan(normOdd)) = 0; 
        normEven(isnan(normEven)) = 0;
        
        [~, peaks] = max(normOdd, [], 2); 
        [~, sIdx] = sort(peaks);
        numN = size(allEven, 1);
        
        qIdx = length(plotQueue) + 1;
        plotQueue(qIdx).matrixToPlot = normOdd(sIdx, :);
        plotQueue(qIdx).rawMatrix = normOdd; 
        plotQueue(qIdx).numN = numN;
        plotQueue(qIdx).titleStr = sprintf('Day %s (Odd)', titleDayStr);
        plotQueue(qIdx).thisCol = thisCol;
        
        qIdx = length(plotQueue) + 1;
        plotQueue(qIdx).matrixToPlot = normEven(sIdx, :);
        plotQueue(qIdx).rawMatrix = normEven;
        plotQueue(qIdx).numN = numN;
        plotQueue(qIdx).titleStr = sprintf('Day %s (Even)', titleDayStr);
        plotQueue(qIdx).thisCol = thisCol;
    end
    
    if isempty(plotQueue), figHandle = []; disp('No valid baseline plots generated.'); return; end
    
    nCols = length(plotQueue);
    
    % layout setup
    figWidth = 280 * nCols + 180;
    figHandle = figure('Position', [100 100 figWidth 580], 'Color', 'w');
    
    leftMargin   = 0.12; 
    rightMargin  = 0.12; 
    topMargin    = 0.08;
    bottomMargin = 0.14; 
    gapCols      = 0.05; 
    gapRows      = 0.06; 
    
    availableH = 1 - topMargin - bottomMargin - gapRows;
    heatmapH   = availableH * 0.72; 
    traceH     = availableH * 0.28; 
    
    availableW = 1 - leftMargin - rightMargin - (nCols-1)*gapCols;
    colW       = availableW / nCols;
    
    heatmapAxes = gobjects(1, nCols);
    traceAxes   = gobjects(1, nCols);
    
    for iCol = 1:nCols
        q = plotQueue(iCol);
        colLeft = leftMargin + (iCol-1) * (colW + gapCols);
        
        % heatmaps
        ax_heat = axes('Position', [colLeft, bottomMargin + traceH + gapRows, colW, heatmapH]);
        heatmapAxes(iCol) = ax_heat;
        
        renderHeatmap(ax_heat, q.matrixToPlot, q.numN, p.Results.TypeToPlot, targetFont);
        title(ax_heat, q.titleStr, 'Color', q.thisCol, 'FontSize', 12, 'FontName', targetFont, 'FontWeight', 'normal');
        
        % popoulation traces 
        ax_trace = axes('Position', [colLeft, bottomMargin, colW, traceH]);
        traceAxes(iCol) = ax_trace;
        
        meanProfile = mean(q.rawMatrix, 1, 'omitnan');
        stdProfile  = std(q.rawMatrix, 0, 1, 'omitnan');
        semProfile  = stdProfile ./ sqrt(sum(~isnan(q.rawMatrix), 1));
        x_vector = 1:length(meanProfile);
        
        hold(ax_trace, 'on');
        x_patch = [x_vector, fliplr(x_vector)];
        y_patch = [(meanProfile + semProfile), fliplr(meanProfile - semProfile)];
        nan_mask = isnan(x_patch) | isnan(y_patch);
        x_patch(nan_mask) = []; y_patch(nan_mask) = [];
        
        fill(ax_trace, x_patch, y_patch, 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        plot(ax_trace, x_vector, meanProfile, 'Color', 'k', 'LineWidth', 1.8);
        
        landmarks = [40 80 120 160];
        for lIdx = 1:4
            xline(ax_trace, landmarks(lIdx), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
        end
        
        xlabel(ax_trace, 'Position (cm)', 'FontName', targetFont, 'FontSize', 12);
        xticks(ax_trace, [40 80 120 160]);
        
        if iCol == 1
%             ylabel(ax_trace, 'Mean \Delta F/F', 'FontName', targetFont, 'FontSize', 12);
            ylabel(ax_trace, 'Mean \Delta F/F', 'FontName', targetFont, 'FontSize', 12);
        else
            set(ax_trace, 'YColor', 'none');
        end
        
        ylim(ax_trace, [0.2, max(meanProfile + semProfile) * 1.15]);
        
        % Break state hold before executing styling properties
        hold(ax_trace, 'off');
        drawnow; 
        
        defaultAxesProperties(ax_trace, true);
        
        % Colorbar settings
        if iCol == nCols
            cb = colorbar(ax_heat, 'Position', [colLeft + colW + 0.02, bottomMargin + traceH + gapRows, 0.02, heatmapH]);
            cb.Ticks = [0.25 0.50 0.75]; 
            cb.TickLabels = {'0.25', '0.50', '0.75'};
            cb.TickDirection = 'out'; 
            cb.Box = 'off';
            cb.FontName = targetFont;
            cb.FontSize = 10;
            
            cb.Label.String = 'Activity (norm.)'; 
            cb.Label.FontName = targetFont;
            cb.Label.FontSize = 12; 
            cb.Label.Rotation = 90;
            cb.Label.VerticalAlignment = 'bottom';
        end
    end
    
    % 
    linkaxes([heatmapAxes, traceAxes], 'x');
    xPadding = 10; 
    xlim(heatmapAxes(1), [1 - xPadding, size(plotQueue(1).matrixToPlot, 2) + xPadding]);
    
    if ~isempty(p.Results.SavePath)
        saveFigureFormats(figHandle, p.Results.SavePath);
    end
end

function renderHeatmap(ax, displayMatrix, numN, cellType, targetFont)
    imagesc(ax, displayMatrix);
    colormap(ax, flipud(gray)); 
    
    set(ax, 'CLim', [0 1], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
    set(ax, 'YTick', [], 'FontName', targetFont, 'FontSize', 12, 'YColor', 'none');
    set(ax, 'XTick', [], 'XTickLabel', [], 'XColor', 'none'); 
    
    hold(ax, 'on');
    landmarks = [40 80 120 160];
    for lIdx = 1:4
        xline(ax, landmarks(lIdx), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
    end
    
    text(ax, -12, numN, sprintf('%d %s', numN, lower(cellType)), ...
        'Rotation', 90, ...
        'FontName', targetFont, ...
        'FontSize', 12, ...
        'FontAngle', 'italic', ...
        'Color', [0.3 0.3 0.3], ...
        'HorizontalAlignment', 'right', ... 
        'VerticalAlignment', 'middle', ...
        'Clipping', 'off');
    hold(ax, 'off');
end