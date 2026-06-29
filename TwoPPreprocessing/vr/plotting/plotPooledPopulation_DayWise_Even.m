function figHandle = plotPooledPopulation_DayWise_Even(allData, targetArea, varargin)
% plotPooledPopulation_DayWise_OddEven: Cross-day timeline visualizer.
% Plots only Even trials SORTED by local matching Odd peaks. Every day panel 
% is independently local range normalized, maintaining perfect square aspect ratios,
% customized gray dotted landmark lines, and a vertical colorbar on the final tile.
%
% Reusable Dependencies: smoothLapActivity, saveFigureFormats

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
    cellType = p.Results.TypeToPlot;
    
    % --- FIELD ALIGNMENT ---
    if all(cellfun(@isempty, {allData.TargetArea})), [allData.TargetArea] = deal(char(targetArea)); end
    if isfield(allData, 'Type') && ~isfield(allData, 'TypeImaged')
        [allData.TypeImaged] = allData.Type;
    elseif all(cellfun(@isempty, {allData.TypeImaged}))
        [allData.TypeImaged] = deal(char(cellType));
    end
    
    daysToPlot = p.Results.DaysToPlot;
    daysToPlot(daysToPlot == 200) = []; 
    nDays = length(daysToPlot);
    
    dayQueue = struct('matrixEven', {}, 'sortIdx', {}, 'numN', {}, 'titleStr', {}, 'titleColor', {});
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
                              strcmpi(string({allData.TypeImaged}), string(cellType)));
        
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
        
        % --- LOCAL PEAK ORDER CALCULATION FROM ODD TRIALS ---
        minOdd = min(allOdd, [], 2, 'omitnan'); 
        maxOdd = max(allOdd, [], 2, 'omitnan'); 
        rangeOdd = maxOdd - minOdd; rangeOdd(rangeOdd == 0) = 1;
        
        normOddLocal = (allOdd - minOdd) ./ rangeOdd;
        normOddLocal(isnan(normOddLocal)) = 0;
        [~, peaks] = max(normOddLocal, [], 2); 
        [~, localSortIdx] = sort(peaks);
        
        numN = size(allEven, 1);
        
        qIdx = length(dayQueue) + 1;
        dayQueue(qIdx).matrixEven  = allEven;
        dayQueue(qIdx).sortIdx     = localSortIdx; 
        dayQueue(qIdx).numN        = numN;
        % dayQueue(qIdx).titleStr    = sprintf('Day %s\n(Even)', titleDayStr);
        dayQueue(qIdx).titleColor  = thisCol;
    end
    
    if isempty(dayQueue), figHandle = []; disp('No valid day-wise records compiled.'); return; end
    
    % --- PLOT RENDERING MATRIX GRID ---
    nPlots = length(dayQueue);
    
    figHandle = figure('Position', [100 100 240*nPlots + 150 360], 'Color', 'w');
    t = tiledlayout(figHandle, 1, nPlots, 'TileSpacing', 'compact', 'Padding', 'loose');
    
    for i = 1:nPlots
        ax = nexttile(t);
        q = dayQueue(i);
        
        % --- LOCAL RANGE NORMALIZATION FOR EVEN TRIALS ---
        minEven = min(q.matrixEven, [], 2, 'omitnan');
        maxEven = max(q.matrixEven, [], 2, 'omitnan');
        rangeEven = maxEven - minEven; rangeEven(rangeEven == 0) = 1;
        
        normalizedEven = (q.matrixEven - minEven) ./ rangeEven;
        normalizedDisplay = normalizedEven(q.sortIdx, :);
        normalizedDisplay(isnan(normalizedDisplay)) = 0;
        
        imagesc(ax, normalizedDisplay);
        colormap(ax, flipud(gray));
        
        axis(ax, 'square');
        
        set(ax, 'CLim', [0.25 0.75], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
        set(ax, 'YTick', [], 'FontName', targetFont, 'FontSize', 12, 'YColor', 'none');
        xlabel(ax, 'Position (cm)', 'FontName', targetFont, 'FontSize', 12);
        xticks(ax, [40 80 120 160]);
        
        title(ax, q.titleStr, 'Color', q.titleColor, 'FontSize', 12, 'FontName', targetFont, 'FontWeight', 'normal');
        
        % --- ADDED: CUSTOM BACKGROUND LANDMARK LINES ---
        hold(ax, 'on');
        landmarks = [40 80 120 160];
        for lIdx = 1:4
            xline(ax, landmarks(lIdx), ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        end
        
        % Vertical Labels positioned along the Y-axis left flank
        text(ax, -12, q.numN, sprintf('%d %s', q.numN, lower(cellType)), ...
            'Rotation', 90, ...
            'FontName', targetFont, ...
            'FontSize', 12, ...
            'FontAngle', 'italic', ...
            'Color', [0.3 0.3 0.3], ...
            'HorizontalAlignment', 'right', ... 
            'VerticalAlignment', 'middle', ...
            'Clipping', 'off');
        
        if i == nPlots
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
            cb.Label.Position = [4.0, 0.5, 0]; 
            cb.Label.HorizontalAlignment = 'center';
            cb.Label.VerticalAlignment = 'bottom';
        end
    end
    

    if ~exist(p.Results.SavePath, 'dir'), mkdir(p.Results.SavePath); end
    saveFigureFormats(figHandle, fullfile(p.Results.SavePath, 'peakbin_halves_within_lap_rsp_visp'));
end
