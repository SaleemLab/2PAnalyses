function figHandle = plotThreeConditions_DifferenceIncluded(allData, targetArea, varargin)
% plotThreeConditions_DifferenceIncluded: Heatmap-Only Production Layout (3-Condition Version).
% Row 1: Absolute Normalized Heatmaps (Grayscale, 4 Columns [BaseEven + 3 Conditions])
% Row 2: Relative Delta Change Heatmaps (Red/White/Blue, 3 Columns)
% Row 3: Mean Delta Change Line Plots (3 Columns)
% Layout mirrors plotConditions_DifferenceIncluded, squeezed for a 4/3-column set.

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
    
    warmColors = [
        0.541, 0.012, 0.012;  % Omit 2 (Deep Blood Red)
        0.824, 0.016, 0.176   % Omit 3 (Cherry Red)
    ];  
    coolColors = [
        0.000, 0.400, 1.000   % Swap 2 3 (Azul Blue)
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
        
        allCondsInDay = {};
        for s = 1:length(rawSessions)
            if isfield(rawSessions(s), 'ConditionData') && ~isempty(rawSessions(s).ConditionData)
                allCondsInDay = unique([allCondsInDay; fieldnames(rawSessions(s).ConditionData)]);
            end
        end
        baseIdx = find(contains(lower(allCondsInDay), 'baseline') | contains(lower(allCondsInDay), 'default'), 1);
        if isempty(baseIdx), baseIdx = 1; end
        baseName = allCondsInDay{baseIdx};
        
        % Three-condition requirement set: 1 swap + 2 omits
        targetSwaps = {'Swap_2_3'};
        targetOmits = {'Omit_2', 'Omit_3'}; 
        requiredConditions = [{baseName}, targetSwaps, targetOmits];
        
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
        if isempty(daySessions), continue; end
        
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
            end
            
            dayMatrixOdd = []; dayMatrixEven = [];
            
            for s = 1:length(daySessions)
                thisSess = daySessions(s);
                if isfield(thisSess, 'FilteredROIs') && ~isempty(thisSess.FilteredROIs)
                    idx = thisSess.FilteredROIs;
                else
                    idx = 1:size(thisSess.ConditionData.(baseName).LapActivity, 1);
                end
                
                lapActivity = thisSess.ConditionData.(cName).LapActivity(idx, :, :);
                if p.Results.ApplySmoothing, lapActivity = smoothLapActivity(lapActivity, 25); end
                
                nTotalLaps = size(lapActivity, 2);
                
                meanOddVals = squeeze(mean(lapActivity(:, 1:2:nTotalLaps, :), 2, 'omitnan'));
                meanEvenVals = squeeze(mean(lapActivity(:, 2:2:nTotalLaps, :), 2, 'omitnan'));
                
                if size(meanOddVals, 2) == 1, meanOddVals = meanOddVals'; meanEvenVals = meanEvenVals'; end
                
                dayMatrixOdd  = vertcat(dayMatrixOdd, meanOddVals); 
                dayMatrixEven = vertcat(dayMatrixEven, meanEvenVals); 
            end
            
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
    
    if isempty(condDataStore), figHandle = []; disp('No sessions matched.'); return; end
    
    %% MASTER SORTING & SUBSET ISOLATION 
    fullBaselineOdd = condDataStore(1).rawMatrices;
    minFullOdd = min(fullBaselineOdd, [], 2, 'omitnan'); 
    maxFullOdd = max(fullBaselineOdd, [], 2, 'omitnan'); 
    rangeFullOdd = maxFullOdd - minFullOdd; rangeFullOdd(rangeFullOdd == 0) = 1;
    
    normFullOddRef = (fullBaselineOdd - minFullOdd) ./ rangeFullOdd;
    normFullOddRef(isnan(normFullOddRef)) = 0;
    [~, fullPeaks] = max(normFullOddRef, [], 2);
    [~, fullSortIdx] = sort(fullPeaks);
    numN = size(fullBaselineOdd, 1);

    tuningRange = max(normFullOddRef, [], 2) - min(normFullOddRef, [], 2);
    keepCells = tuningRange > 0.3;
    fullSortIdx_filtered = fullSortIdx(keepCells(fullSortIdx)); %#ok<NASGU>
    numN_filtered = sum(keepCells); %#ok<NASGU>
    fprintf('Cells with clear baseline tuning: %d / %d\n', numN_filtered, numN);
    
    % EXTRACT BASELINE REFERENCE ARRAYS 
    baseEvenIdx = find(contains(lower({condDataStore.conditionNames}), 'baseline even'), 1);
    rawBaseEven = condDataStore(baseEvenIdx).rawMatrices;
    normBaseEven = (rawBaseEven - minFullOdd) ./ rangeFullOdd;
    normBaseEven(isnan(normBaseEven)) = 0;
    
    %% 
    % Subplot positions sized for a 4-column absolute row and a 3-column delta/line set.
    row1_Positions = {[0.10, 0.64, 0.16, 0.28], [0.32, 0.64, 0.16, 0.28], [0.54, 0.64, 0.16, 0.28], [0.76, 0.64, 0.16, 0.28]}; % 4 Absolute Columns
    row2_Positions = {[0.32, 0.28, 0.16, 0.28], [0.54, 0.28, 0.16, 0.28], [0.76, 0.28, 0.16, 0.28]};                         % 3 Delta Columns (same height as Row 1)
    row3_Positions = {[0.32, 0.06, 0.16, 0.14], [0.54, 0.06, 0.16, 0.14], [0.76, 0.06, 0.16, 0.14]};                         % 3 Mean Delta Line Columns
    
    figHandle = figure('Position', [50 50 1150 780], 'Color', 'w');
    
    % 
    cMapLen = 256;
    b = [linspace(1, 1, cMapLen/2), linspace(1, 0, cMapLen/2)];
    g = [linspace(0, 1, cMapLen/2), linspace(1, 0, cMapLen/2)];
    r = [linspace(0, 1, cMapLen/2), linspace(1, 1, cMapLen/2)];
    redWhiteBlueMap = [b', g', r']; 
    
    row1_Indices = [2, 3, 4, 5]; % Base Even + 3 active tracking manipulations
    row2_Tokens  = {'swap.*2.*3', 'omit.*2', 'omit.*3'};
    displayNames = {'Swap 2 3', 'Omit 2', 'Omit 3'};
    colors = {coolColors(1,:), warmColors(1,:), warmColors(2,:)};

    %% --- ROW 1: ABSOLUTE INFERRED SPIKE AMPLITUDE HEATMAPS (GRAYSCALE) ---
  
    for i = 1:4
        axAbs = axes('Position', row1_Positions{i}); 
        currentCond = condDataStore(row1_Indices(i));
        
        localIsolateMatrix = currentCond.rawMatrices;
        normalizedDisplay = (localIsolateMatrix - minFullOdd) ./ rangeFullOdd;
        normalizedDisplay(isnan(normalizedDisplay)) = 0;
        
        imagesc(axAbs, [0, 200], [1, numN], normalizedDisplay(fullSortIdx, :));
        colormap(axAbs, flipud(gray));
        set(axAbs, 'CLim', [0.25 0.75], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
        set(axAbs, 'XTick', [40 80 120 160], 'XTickLabel', {'40', '80', '120', '160'}, 'FontName', targetFont, 'FontSize', 10);
        
        if i == 1
            title(axAbs, 'Base Even', 'Color', [0 0 0], 'FontSize', 11, 'FontName', targetFont, 'FontWeight', 'bold');
            ylabel(axAbs, sprintf('Sorted Background %s\n(n = %d)', cellType, numN), 'FontName', targetFont, 'FontSize', 12);
            text(axAbs, -15, numN, sprintf('%d background %s', numN, lower(cellType)), ...
                'Rotation', 90, 'FontName', targetFont, 'FontSize', 11, 'FontAngle', 'italic', ...
                'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Clipping', 'off');
        else
            set(axAbs, 'YColor', 'none', 'YTick', []);
            cleanTitle = strrep(currentCond.conditionNames, '_', ' ');
            title(axAbs, cleanTitle, 'Color', currentCond.titleColors, 'FontSize', 11, 'FontName', targetFont, 'FontWeight', 'bold');
        end
        
        hold(axAbs, 'on');
        for posBin = [40 80 120 160]
            xline(axAbs, posBin, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        end
    end
    
    % Absolute Grayscale Colorbar Flank
    axAbsCB = axes('Position', [0.94, 0.64, 0.015, 0.28]); 
    set(axAbsCB, 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
    colormap(axAbsCB, flipud(gray));
    cbAbs = colorbar(axAbsCB, 'Location', 'eastoutside');
    cbAbs.Position = [0.94, 0.64, 0.015, 0.28];
    cbAbs.Ticks = [0.25 0.50 0.75]; cbAbs.TickLabels = {'0.25', '0.50', '0.75'};
    cbAbs.TickDirection = 'out'; cbAbs.Box = 'off'; cbAbs.FontName = targetFont; cbAbs.FontSize = 10;
    cbAbs.Label.String = 'Inferred Spike Amplitude (norm.)'; cbAbs.Label.FontName = targetFont; cbAbs.Label.FontSize = 11;

    %%RELATIVE CHANGE DELTA HEATMAP PANELS (RED/WHITE/BLUE) ---
   
    for i = 1:3
        matchIdx = find(~cellfun(@isempty, regexpi({condDataStore.conditionNames}, row2_Tokens{i})), 1);
        if isempty(matchIdx), continue; end
        currentCond = condDataStore(matchIdx);
        
        localMatrix = currentCond.rawMatrices;
        normCond = (localMatrix - minFullOdd) ./ rangeFullOdd;
        normCond(isnan(normCond)) = 0; 
        
        % Subtraction Logic (Condition - Baseline Even)
        diffMat = normCond - normBaseEven;
        
        axHM = axes('Position', row2_Positions{i});
        imagesc(axHM, [0, 200], [1, numN], diffMat(fullSortIdx, :));
        colormap(axHM, flipud(redWhiteBlueMap)); % Red = Increase (+), Blue = Decrease (-)
        set(axHM, 'CLim', [-0.4 0.4], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
        
        set(axHM, 'XTick', [40 80 120 160], 'XTickLabel', {});
        title(axHM, ['\Delta ' displayNames{i}], 'Color', colors{i}, 'FontSize', 11, 'FontName', targetFont, 'FontWeight', 'bold');
        
        landmarks = [40 80 120 160];
        currentNameLow = lower(currentCond.conditionNames);
        hold(axHM, 'on');
        for lIdx = 1:4
            posBin = landmarks(lIdx);
            isTargetLandmark = false; 
            
         
            if contains(currentNameLow, 'swap')
                tokenParts = regexp(row2_Tokens{i}, '\d', 'match');
                if any(strcmpi(num2str(lIdx), tokenParts)), isTargetLandmark = true; end
            else
                if contains(currentNameLow, num2str(lIdx)), isTargetLandmark = true; end
            end
            
            if isTargetLandmark
                xline(axHM, posBin, '--', 'Color', colors{i}, 'LineWidth', 2.0);
            else
                xline(axHM, posBin, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            end
        end
        
        if i == 1
            ylabel(axHM, sprintf('Sorted Background %s\n(n = %d)', cellType, numN), 'FontName', targetFont, 'FontSize', 12);
        else
            set(axHM, 'YColor', 'none', 'YTick', []);
        end

        %% --- ROW 3: MEAN DELTA CHANGE LINE PLOT ---
        meanDiffTrace = mean(diffMat, 1, 'omitnan');
        nBins = numel(meanDiffTrace);
        xPositions = linspace(0, 200, nBins);

        axMean = axes('Position', row3_Positions{i});
        hold(axMean, 'on');
        yline(axMean, 0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        plot(axMean, xPositions, meanDiffTrace, '-', 'Color', colors{i}, 'LineWidth', 1.5);
        xlim(axMean, [1, 200]);
        set(axMean, 'Box', 'off', 'TickDir', 'out');
        % Set LAST, after all other formatting, so nothing resets it back to auto
        set(axMean, 'YLim', [-0.3, 0.3], 'YLimMode', 'manual','XTick', [-0.3,0, 0.3]);
        set(axMean, 'XTick', [40 80 120 160], 'XTickLabel', {'40', '80', '120', '160'}, 'FontName', targetFont, 'FontSize', 10);
        xlabel(axMean, 'Position (cm)', 'FontName', targetFont, 'FontSize', 11);

        if i == 1
            ylabel(axMean, sprintf('Mean activity\n(Cond - Base)'), 'FontName', targetFont, 'FontSize', 11);
        end
    end
    
    %% 
    axCB = axes('Position', [0.94, 0.28, 0.015, 0.28]); 
    set(axCB, 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
    colormap(axCB, flipud(redWhiteBlueMap));
    
    cbDiff = colorbar(axCB, 'Location', 'eastoutside');
    cbDiff.Position = [0.94, 0.28, 0.015, 0.28];
    cbDiff.Ticks = [-0.4 0 0.4]; cbDiff.TickLabels = {'Decrease (-)', 'No Change (0)', 'Increase (+)'};
    cbDiff.TickDirection = 'out'; cbDiff.Box = 'off'; cbDiff.FontName = targetFont; cbDiff.FontSize = 10;
    cbDiff.Label.String = '\Delta Inferred Spike Amplitude (Condition - Baseline Even)';
    cbDiff.Label.FontName = targetFont; cbDiff.Label.FontSize = 11;
    
    if ~isempty(p.Results.SavePath), saveFigureFormats(figHandle, p.Results.SavePath); end
end