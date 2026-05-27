function figHandle = plotPooledPopulation_AcrossConditions(allData, targetArea, varargin)
% plotPooledPopulation_AcrossConditions: Strict Intersection Version.
% This version guarantees that ALL columns contain the EXACT same ROIs by completely 
% dropping any session that doesn't contain the full battery of required conditions.

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
    
    % warmColors = [0.9 0.2 0.2; 1.0 0.6 0.0; 0.8 0.0 0.4];  
    % coolColors = [0.0 0.45 0.74; 0.0 0.8 0.8; 0.2 0.2 0.6]; 
    

    warmColors = [
    0.541, 0.012, 0.012;  % Blood
    0.824, 0.016, 0.176;  % Cherry
    0.898, 0.357, 0.357   % Faded Red
    ];  
    

   
    coolColors = [
    0.000, 0.400, 1.000;  % Azul
    0.000, 0.500, 1.000   % Azure
    ];
    
    
    if all(cellfun(@isempty, {allData.TargetArea})), [allData.TargetArea] = deal(char(targetArea)); end
    if isfield(allData, 'Type') && ~isfield(allData, 'TypeImaged')
        [allData.TypeImaged] = allData.Type;
    elseif all(cellfun(@isempty, {allData.TypeImaged}))
        [allData.TypeImaged] = deal(char(cellType));
    end
    
    daysToPlot = p.Results.DaysToPlot;
    daysToPlot(daysToPlot == 200) = []; 
    nDays = length(daysToPlot);
    
    condDataStore = struct('rawMatrices', {}, 'conditionNames', {}, 'condTypes', {}, 'titleColors', {});
    
    for d = 1:nDays
        day = daysToPlot(d);
        if day == 5
            dayMask = ([allData.Day] == 5 | [allData.Day] == 200);
        else
            dayMask = ([allData.Day] == day);
        end
        
        rawSessions = allData(dayMask & ...
                              strcmpi(string({allData.TargetArea}), string(targetArea)) & ...
                              strcmpi(string({allData.TypeImaged}), string(cellType)));
        
        if isempty(rawSessions), continue; end
        
        % Dynamic baseline name discovery
        allCondsInDay = {};
        for s = 1:length(rawSessions)
            if isfield(rawSessions(s), 'ConditionData') && ~isempty(rawSessions(s).ConditionData)
                allCondsInDay = unique([allCondsInDay; fieldnames(rawSessions(s).ConditionData)]);
            end
        end
        baseIdx = find(contains(lower(allCondsInDay), 'baseline') | contains(lower(allCondsInDay), 'default'), 1);
        if isempty(baseIdx), baseIdx = 1; end
        baseName = allCondsInDay{baseIdx};
        
        % Set up the required list based on imaging type
        if strcmpi(cellType, 'Boutons')
            targetSwaps = {'Swap_2_3'};
            targetOmits = {'Omit_2', 'Omit_3'};
        else
            targetSwaps = {'Swap_2_3', 'Swap_3_4'};
            targetOmits = {'Omit_2', 'Omit_3', 'Omit_4'};
        end
        requiredConditions = [{baseName}, targetSwaps, targetOmits];
        
        % --- CRITICAL FIX: STAGE-GATE FILTER FOR SESSIONS ---
        % Keep a session ONLY if it contains 100% of the required conditions
        validSessionMask = false(1, length(rawSessions));
        for s = 1:length(rawSessions)
            if isfield(rawSessions(s), 'ConditionData') && ~isempty(rawSessions(s).ConditionData)
                sessConds = fieldnames(rawSessions(s).ConditionData);
                if all(ismember(requiredConditions, sessConds))
                    validSessionMask(s) = true;
                end
            end
        end
        daySessions = rawSessions(validSessionMask);
        
        if isempty(daySessions)
            warning('No sessions contained all required conditions for Day %d. Skipping.', day);
            continue;
        end
        
        % Process the remaining cleanly matched dataset
        orderedConds = requiredConditions;
        omitCount = 1; swapCount = 1;
        
        for c = 1:length(orderedConds)
            cName = orderedConds{c};
            isBaseBlock = (c == 1);
            nameLow = lower(cName);
            
            if isBaseBlock
                cTypeOdd = 'baseline_odd'; cTypeEven = 'baseline_even';
                colColorOdd = [0 0 0]; colColorEven = [0 0 0];
            elseif contains(nameLow, 'swap')
                cTypeOdd = 'swap'; cTypeEven = 'swap';
                colColorOdd = coolColors(mod(swapCount-1, size(coolColors,1))+1, :);
                colColorEven = colColorOdd; swapCount = swapCount + 1;
            elseif contains(nameLow, 'omit')
                cTypeOdd = 'omit'; cTypeEven = 'omit';
                colColorOdd = warmColors(mod(omitCount-1, size(warmColors,1))+1, :);
                colColorEven = colColorOdd; omitCount = omitCount + 1;
            else
                cTypeOdd = 'other'; cTypeEven = 'other';
                colColorOdd = [0.5 0.5 0.5]; colColorEven = colColorOdd;
            end
            
            dayMatrixOdd = []; dayMatrixEven = [];
            
            for s = 1:length(daySessions)
                thisSess = daySessions(s);
                
                if isfield(thisSess, 'FilteredROIs') && ~isempty(thisSess.FilteredROIs)
                    idx = thisSess.FilteredROIs;
                else
                    idx = 1:size(thisSess.ConditionData.(baseName).LapActivity, 1);
                end
                
                lapActivity = thisSess.ConditionData.(cName).LapActivity;
                lapActivity = lapActivity(idx, :, :);
                
                if p.Results.ApplySmoothing
                    lapActivity = smoothLapActivity(lapActivity);
                end
                
                nTotalLaps = size(lapActivity, 2);
                dayMatrixOdd  = vertcat(dayMatrixOdd, squeeze(mean(lapActivity(:, 1:2:nTotalLaps, :), 2, 'omitnan'))); %#ok<AGROW>
                dayMatrixEven = vertcat(dayMatrixEven, squeeze(mean(lapActivity(:, 2:2:nTotalLaps, :), 2, 'omitnan'))); %#ok<AGROW>
            end
            
            if size(dayMatrixOdd, 2) == 1, dayMatrixOdd = dayMatrixOdd'; dayMatrixEven = dayMatrixEven'; end
            
            if isBaseBlock
                qEnd = length(condDataStore) + 1;
                condDataStore(qEnd).rawMatrices = dayMatrixOdd;
                condDataStore(qEnd).conditionNames = 'Baseline Odd';
                condDataStore(qEnd).condTypes = cTypeOdd;
                condDataStore(qEnd).titleColors = colColorOdd;
                
                qEnd = length(condDataStore) + 1;
                condDataStore(qEnd).rawMatrices = dayMatrixEven;
                condDataStore(qEnd).conditionNames = 'Baseline Even';
                condDataStore(qEnd).condTypes = cTypeEven;
                condDataStore(qEnd).titleColors = colColorEven;
            else
                combinedOtherLaps = mean(cat(3, dayMatrixOdd, dayMatrixEven), 3, 'omitnan');
                qEnd = length(condDataStore) + 1;
                condDataStore(qEnd).rawMatrices = combinedOtherLaps;
                condDataStore(qEnd).conditionNames = strrep(cName, '_', ' ');
                condDataStore(qEnd).condTypes = cTypeOdd;
                condDataStore(qEnd).titleColors = colColorOdd;
            end
        end
    end
    
    if isempty(condDataStore), figHandle = []; disp('No sessions matched the filter.'); return; end
    
    % --- MASTER UNIFIED SORTING ---
    referenceMatrix = condDataStore(1).rawMatrices;
    minOdd = min(referenceMatrix, [], 2, 'omitnan'); 
    maxOdd = max(referenceMatrix, [], 2, 'omitnan'); 
    rangeOdd = maxOdd - minOdd; rangeOdd(rangeOdd == 0) = 1;
    
    normOddRef = (referenceMatrix - minOdd) ./ rangeOdd;
    normOddRef(isnan(normOddRef)) = 0;
    [~, peaks] = max(normOddRef, [], 2);
    [~, sortIdx] = sort(peaks);
    numN = size(referenceMatrix, 1);
    
    % Strip out Baseline Odd panel from plotting queues
    condDataStore(1) = [];
    nPlots = length(condDataStore);
    
    figHandle = figure('Position', [50 100 230*nPlots + 150 400], 'Color', 'w');
    t = tiledlayout(figHandle, 1, nPlots, 'TileSpacing', 'compact', 'Padding', 'loose');
    
    for i = 1:nPlots
        ax = nexttile(t);
        currentCond = condDataStore(i);
        localMatrix = currentCond.rawMatrices;
        
        normalizedDisplay = (localMatrix - minOdd) ./ rangeOdd;
        normalizedDisplay(isnan(normalizedDisplay)) = 0;
        
        imagesc(ax, normalizedDisplay(sortIdx, :));
        colormap(ax, flipud(gray));
        
        set(ax, 'CLim', [0.25 0.75], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
        set(ax, 'YTick', [], 'FontName', targetFont, 'FontSize', 12, 'YColor', 'none');
        xlabel(ax, 'Position (cm)', 'FontName', targetFont, 'FontSize', 12);
        xticks(ax, [40 80 120 160]);
        
        title(ax, sprintf('%s', currentCond.conditionNames), ...
            'Color', currentCond.titleColors, 'FontSize', 12, 'FontName', targetFont, 'FontWeight', 'normal');
        
        % Landmarks
        landmarks = [40 80 120 160];
        currentNameLow = lower(currentCond.conditionNames);
        hold(ax, 'on');
        for lIdx = 1:4
            posBin = landmarks(lIdx);
            isTargetLandmark = false; 

            if ~contains(currentCond.condTypes, 'baseline') && contains(currentNameLow, num2str(lIdx))
                isTargetLandmark = true;
            end
            if strcmp(currentCond.condTypes, 'omit') && ~any(regexp(currentNameLow, '[1-4]')) && posBin == 120
                isTargetLandmark = true;
            end

            % 
            if isTargetLandmark
                xline(ax, posBin, '--', 'Color', currentCond.titleColors, 'LineWidth', 2.0);
            else
                xline(ax, posBin, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            end
        end

        if i == 1
            text(ax, -10, numN, sprintf('%d %s', numN, lower(cellType)), ...
                'Rotation', 90, 'FontName', targetFont, 'FontSize', 12, 'FontAngle', 'italic', ...
                'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Clipping', 'off');

            % if numN >= 50
            %     y_start = numN * 0.2; y_end = y_start + 100; y_mid = y_start + 100;  
            %     line(ax, [-12 -12], [y_start y_end], 'Color', 'k', 'LineWidth', 2, 'Clipping', 'off');
            %     text(ax, -20, y_mid+40, sprintf('100 %s', lower(cellType)), 'Rotation', 90, ...
            %         'FontName', targetFont, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...   
            %         'FontSize', 12, 'FontWeight', 'normal', 'Clipping', 'off');
            % end
        end
        
        if i == nPlots
            cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Ticks = [0.25 0.50 0.75]; 
            cb.TickLabels = {'0.25', '0.50', '0.75'}; cb.TickDirection = 'out'; cb.Box = 'off';
            cb.FontName = targetFont; cb.FontSize = 10;
            cb.Label.String = 'Activity (norm.)'; cb.Label.FontName = targetFont; cb.Label.FontSize = 12; 
            cb.Label.Rotation = 90; cb.Label.Units = 'normalized'; cb.Label.Position = [4, 0.5, 0]; 
            cb.Label.HorizontalAlignment = 'center'; cb.Label.VerticalAlignment = 'bottom';
        end
    end
    
    if ~isempty(p.Results.SavePath), saveFigureFormats(figHandle, p.Results.SavePath); end
end