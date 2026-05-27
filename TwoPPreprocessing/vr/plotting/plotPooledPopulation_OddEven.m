function figHandle = plotPooledPopulation_OddEven(allData, targetArea, varargin)
% plotPooledPopulation_OddEven: Baseline function plotting side-by-side Odd
% and Even lap profiles. Uses 'smoothLapActivity' for trace smoothing and 
% 'saveFigureFormats' for vector-safe graphics export. Sets uniform font parameters.
    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    addRequired(p, 'targetArea', @(x) ischar(x) || isstring(x));
    
    addParameter(p, 'DaysToPlot', [1, 2, 3, 4, 5], @isnumeric);
    addParameter(p, 'TypeToPlot', 'Boutons', @(x) ischar(x) || isstring(x));
    addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x)); 
    addParameter(p, 'ApplySmoothing', true, @islogical);
    addParameter(p, 'FontName', 'Arial', @(x) ischar(x) || isstring(x)); % Custom font configuration
    
    parse(p, allData, targetArea, varargin{:});
    
    targetFont = p.Results.FontName;
    
    % 
    if all(cellfun(@isempty, {allData.TargetArea})), [allData.TargetArea] = deal(char(targetArea)); end
    if isfield(allData, 'Type') && ~isfield(allData, 'TypeImaged')
        [allData.TypeImaged] = allData.Type;
    elseif all(cellfun(@isempty, {allData.TypeImaged}))
        [allData.TypeImaged] = deal(char(p.Results.TypeToPlot));
    end
    
    daysToPlot = p.Results.DaysToPlot;
    daysToPlot(daysToPlot == 200) = []; 
    nDays = length(daysToPlot);
    
    plotQueue = struct('matrixToPlot', {}, 'numN', {}, 'titleStr', {}, 'thisCol', {});
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
            allOdd  = vertcat(allOdd, squeeze(mean(lapActivity(:, 1:2:nTotalLaps, :), 2, 'omitnan'))); %#ok<AGROW>
            allEven = vertcat(allEven, squeeze(mean(lapActivity(:, 2:2:nTotalLaps, :), 2, 'omitnan'))); %#ok<AGROW>
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
        plotQueue(qIdx).numN = numN;
        % plotQueue(qIdx).titleStr = sprintf('Day %s (Odd, n=%d)', titleDayStr, numN);
        plotQueue(qIdx).thisCol = thisCol;
        
        qIdx = length(plotQueue) + 1;
        plotQueue(qIdx).matrixToPlot = normEven(sIdx, :);
        plotQueue(qIdx).numN = numN;
        % plotQueue(qIdx).titleStr = sprintf('Day %s (Even)', titleDayStr);
        plotQueue(qIdx).thisCol = thisCol;
    end
    
    if isempty(plotQueue), figHandle = []; disp('No valid baseline plots generated.'); return; end
    
    % --- GRID DRAWING ---
    nRows = 1;
    nCols = length(plotQueue);
    
    % Expanded width dimensions to protect the upgraded 12pt fonts from edge compression
    figHandle = figure('Position', [100 100 230*nCols + 150 400], 'Color', 'w');
    t = tiledlayout(figHandle, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'loose');
    
    for iTile = 1:nCols
        ax = nexttile(t);
        q = plotQueue(iTile);
        
        renderHeatmap(ax, q.matrixToPlot, q.numN, p.Results.TypeToPlot, targetFont);
        title(ax, q.titleStr, 'Color', q.thisCol, 'FontSize', 12, 'FontName', targetFont, 'FontWeight', 'normal');
        
        if iTile == nCols
            cb = colorbar(ax);
            cb.Layout.Tile = 'east'; 
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
            cb.Label.Units = 'normalized';
            cb.Label.Position = [4, 0.5, 0]; 
            cb.Label.HorizontalAlignment = 'center';
            cb.Label.VerticalAlignment = 'bottom';
        end
    end
    
    if ~isempty(p.Results.SavePath)
        saveFigureFormats(figHandle, p.Results.SavePath);
    end
end

% --- HELPER PLOTTING SUBFUNCTION ---
function renderHeatmap(ax, displayMatrix, numN, cellType, targetFont)
    imagesc(ax, displayMatrix);
    colormap(ax, flipud(gray)); 
    
    set(ax, 'CLim', [0.25 0.75], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
    set(ax, 'YTick', [], 'FontName', targetFont, 'FontSize', 12, 'YColor', 'none');
    xlabel(ax, 'Position (cm)', 'FontName', targetFont, 'FontSize', 12);
    xticks(ax, [40 80 120 160]);
    
    % --- FIXED: ADDED GRAY DOTTED LANDMARK LINES ---
    hold(ax, 'on');
    landmarks = [40 80 120 160];
    for lIdx = 1:4
        xline(ax, landmarks(lIdx), ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    end
    
    % Vertical text alignment labels along the Y-axis left flank
    text(ax, -12, numN, sprintf('%d %s', numN, lower(cellType)), ...
        'Rotation', 90, ...
        'FontName', targetFont, ...
        'FontSize', 12, ...
        'FontAngle', 'italic', ...
        'Color', [0.3 0.3 0.3], ...
        'HorizontalAlignment', 'right', ... 
        'VerticalAlignment', 'middle', ...
        'Clipping', 'off');
   
    % if numN >= 100
    %     y_start = numN * 0.2; 
    %     y_end   = y_start + 100; 
    %     y_mid   = y_start + 100;  
    % 
    %     x_line = -12; 
    %     x_text = -20; 
    % 
    %     unitStr = lower(cellType);
    %     scaleBarText = sprintf('100 %s', unitStr);
    % 
    %     line(ax, [x_line x_line], [y_start y_end], 'Color', 'k', 'LineWidth', 2, 'Clipping', 'off');
    % 
    %     text(ax, x_text, y_mid+40, scaleBarText, 'Rotation', 90, ...
    %         'FontName', targetFont, ...
    %         'HorizontalAlignment', 'right', ... 
    %         'VerticalAlignment', 'middle', ...   
    %         'FontSize', 12, 'FontWeight', 'normal', 'Clipping', 'off');
    % end
end