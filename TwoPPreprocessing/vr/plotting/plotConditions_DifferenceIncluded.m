function figHandle = plotConditions_DifferenceIncluded(allData, targetArea, varargin)
% plotSixConditions_CombinedHeatmaps: Ultimate Heatmap-Only Production Layout.
% Row 1: Absolute Normalized Heatmaps (Grayscale, 6 Columns [BaseEven + 5 Conditions])
% Row 2: Relative Delta Change Heatmaps (Red/White/Blue, 5 Columns)
% Squeezed layout dimensions to match thin, vertically elongated subplots.

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
        0.824, 0.016, 0.176;  % Omit 3 (Cherry Red)
        0.950, 0.350, 0.000   % Omit 4 (Burnt Orange)
    ];  
    coolColors = [
        0.000, 0.400, 1.000;  % Swap 2 3 (Azul Blue)
        0.000, 0.700, 0.600   % Swap 3 4 (Teal Blue)
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
        
        % 
        targetSwaps = {'Swap_2_3', 'Swap_3_4'};
        targetOmits = {'Omit_2', 'Omit_3', 'Omit_4'}; 
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
                if p.Results.ApplySmoothing, lapActivity = smoothLapActivity(lapActivity); end
                
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
    
    % EXTRACT BASELINE REFERENCE ARRAYS 
    baseEvenIdx = find(contains(lower({condDataStore.conditionNames}), 'baseline even'), 1);
    rawBaseEven = condDataStore(baseEvenIdx).rawMatrices;
    normBaseEven = (rawBaseEven - minFullOdd) ./ rangeFullOdd;
    normBaseEven(isnan(normBaseEven)) = 0;
    
    %% 
    % Subplot widths reduced to 0.11 to cleanly space out a 6-column layout.
    row1_Positions = {[0.08, 0.54, 0.11, 0.38], [0.21, 0.54, 0.11, 0.38], [0.34, 0.54, 0.11, 0.38], [0.47, 0.54, 0.11, 0.38], [0.60, 0.54, 0.11, 0.38], [0.73, 0.54, 0.11, 0.38]}; % 6 Absolute Columns
    row2_Positions = {[0.21, 0.10, 0.11, 0.38], [0.34, 0.10, 0.11, 0.38], [0.47, 0.10, 0.11, 0.38], [0.60, 0.10, 0.11, 0.38], [0.73, 0.10, 0.11, 0.38]};                         % 5 Delta Columns
    
    figHandle = figure('Position', [30 30 1450 780], 'Color', 'w');
    
    % 
    cMapLen = 256;
    b = [linspace(1, 1, cMapLen/2), linspace(1, 0, cMapLen/2)];
    g = [linspace(0, 1, cMapLen/2), linspace(1, 0, cMapLen/2)];
    r = [linspace(0, 1, cMapLen/2), linspace(1, 1, cMapLen/2)];
    redWhiteBlueMap = [b', g', r']; 
    
    row1_Indices = [2, 3, 4, 5, 6, 7]; % Base Even + 5 active tracking manipulations
    row2_Tokens  = {'swap.*2.*3', 'swap.*3.*4', 'omit.*2', 'omit.*3', 'omit.*4'};
    displayNames = {'Swap 2 3', 'Swap 3 4', 'Omit 2', 'Omit 3', 'Omit 4'};
    colors = {coolColors(1,:), coolColors(2,:), warmColors(1,:), warmColors(2,:), warmColors(3,:)};

    %% --- ROW 1: ABSOLUTE INFERRED SPIKE AMPLITUDE HEATMAPS (GRAYSCALE) ---
  
    for i = 1:6
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
    axAbsCB = axes('Position', [0.87, 0.54, 0.012, 0.38]); 
    set(axAbsCB, 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
    colormap(axAbsCB, flipud(gray));
    cbAbs = colorbar(axAbsCB, 'Location', 'eastoutside');
    cbAbs.Position = [0.87, 0.54, 0.012, 0.38];
    cbAbs.Ticks = [0.25 0.50 0.75]; cbAbs.TickLabels = {'0.25', '0.50', '0.75'};
    cbAbs.TickDirection = 'out'; cbAbs.Box = 'off'; cbAbs.FontName = targetFont; cbAbs.FontSize = 10;
    cbAbs.Label.String = 'Inferred Spike Amplitude (norm.)'; cbAbs.Label.FontName = targetFont; cbAbs.Label.FontSize = 11;

    %%RELATIVE CHANGE DELTA HEATMAP PANELS (RED/WHITE/BLUE) ---
   
    for i = 1:5
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
        
        set(axHM, 'XTick', [40 80 120 160], 'XTickLabel', {'40', '80', '120', '160'}, 'FontName', targetFont, 'FontSize', 10);
        title(axHM, ['\Delta ' displayNames{i}], 'Color', colors{i}, 'FontSize', 11, 'FontName', targetFont, 'FontWeight', 'bold');
        xlabel(axHM, 'Position (cm)', 'FontName', targetFont, 'FontSize', 11);
        
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
    end
    
    %% 
    axCB = axes('Position', [0.87, 0.10, 0.012, 0.38]); 
    set(axCB, 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
    colormap(axCB, flipud(redWhiteBlueMap));
    
    cbDiff = colorbar(axCB, 'Location', 'eastoutside');
    cbDiff.Position = [0.87, 0.10, 0.012, 0.38];
    cbDiff.Ticks = [-0.4 0 0.4]; cbDiff.TickLabels = {'Decrease (-)', 'No Change (0)', 'Increase (+)'};
    cbDiff.TickDirection = 'out'; cbDiff.Box = 'off'; cbDiff.FontName = targetFont; cbDiff.FontSize = 10;
    cbDiff.Label.String = '\Delta Inferred Spike Amplitude (Condition - Baseline Even)';
    cbDiff.Label.FontName = targetFont; cbDiff.Label.FontSize = 11;
    
    if ~isempty(p.Results.SavePath), saveFigureFormats(figHandle, p.Results.SavePath); end
end