function figHandle = plotBackgroundROIs_DifferenceOnly(allData, targetArea, varargin)

    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    addRequired(p, 'targetArea', @(x) ischar(x) || isstring(x));
    addParameter(p, 'DaysToPlot', [1, 2, 3, 4, 5], @isnumeric);
    addParameter(p, 'TypeToPlot', 'Boutons', @(x) ischar(x) || isstring(x));
    addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ApplySmoothing', true, @islogical);
    addParameter(p, 'FontName', 'Arial', @(x) ischar(x) || isstring(x));
    addParameter(p, 'RequiredConditions', {}, @iscell);
    parse(p, allData, targetArea, varargin{:});

    targetFont = p.Results.FontName;
    cellType   = p.Results.TypeToPlot;

    warmColors = [0.541, 0.012, 0.012; 0.824, 0.016, 0.176; 0.600, 0.000, 0.400];
    coolColors = [0.000, 0.400, 1.000; 0.000, 0.650, 0.650];

    if all(cellfun(@isempty, {allData.TargetArea})), [allData.TargetArea] = deal(char(targetArea)); end
    if isfield(allData, 'Type') && ~isfield(allData, 'TypeImaged')
        [allData.TypeImaged] = allData.Type;
    elseif all(cellfun(@isempty, {allData.TypeImaged}))
        [allData.TypeImaged] = deal(char(cellType));
    end

    daysToPlot = p.Results.DaysToPlot;
    daysToPlot(daysToPlot == 200) = [];
    nDays = length(daysToPlot);

    condDataStore = struct('rawMatrices', {}, 'conditionNames', {}, 'condTypes', {}, 'titleColors', {}, 'hasData', {});

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

        targetSwaps    = {'Swap_2_3', 'Swap_3_4'};
        targetOmits    = {'Omit_2', 'Omit_3', 'Omit_4'};
        allTargetConds = [targetSwaps, targetOmits];

        % Build required conditions list
        requiredForSelection = p.Results.RequiredConditions;
        if isempty(requiredForSelection)
            requiredForSelection = {baseName};
        else
            requiredForSelection = unique([{baseName}, requiredForSelection(:)']);
        end

        fprintf('  Requiring conditions: %s\n', strjoin(requiredForSelection, ', '));

        % Valid sessions must have all required conditions
        validSessionMask = false(1, length(rawSessions));
        for s = 1:length(rawSessions)
            if isfield(rawSessions(s), 'ConditionData') && ~isempty(rawSessions(s).ConditionData)
                sessConds = fieldnames(rawSessions(s).ConditionData);
                if all(ismember(requiredForSelection, sessConds))
                    validSessionMask(s) = true;
                end
            end
        end

        nValid = sum(validSessionMask);
        fprintf('  Day %d: %d/%d sessions pass required conditions\n', day, nValid, length(rawSessions));

        daySessions = rawSessions(validSessionMask);
        if isempty(daySessions), continue; end

        orderedConds = [{baseName}, allTargetConds];
        omitCount = 1; swapCount = 1;

        for c = 1:length(orderedConds)
            cName   = orderedConds{c};
            isBase  = (c == 1);
            nameLow = lower(cName);

            if isBase
                colColorOdd = [0 0 0]; colColorEven = [0 0 0];
                cTypeOdd = 'baseline_odd'; cTypeEven = 'baseline_even';
            elseif contains(nameLow, 'swap')
                colColorOdd  = coolColors(mod(swapCount-1, size(coolColors,1))+1, :);
                colColorEven = colColorOdd; cTypeOdd = 'swap'; cTypeEven = 'swap';
                swapCount    = swapCount + 1;
            elseif contains(nameLow, 'omit')
                colColorOdd  = warmColors(mod(omitCount-1, size(warmColors,1))+1, :);
                colColorEven = colColorOdd; cTypeOdd = 'omit'; cTypeEven = 'omit';
                omitCount    = omitCount + 1;
            end

            dayMatrixOdd = []; dayMatrixEven = [];
            hasData      = false;

            for s = 1:length(daySessions)
                thisSess  = daySessions(s);
                sessConds = fieldnames(thisSess.ConditionData);

                if isfield(thisSess, 'FilteredROIs') && ~isempty(thisSess.FilteredROIs)
                    idx = thisSess.FilteredROIs;
                else
                    idx = 1:size(thisSess.ConditionData.(baseName).LapActivity, 1);
                end

                nBins = size(thisSess.ConditionData.(baseName).LapActivity, 3);
                nLaps = size(thisSess.ConditionData.(baseName).LapActivity, 2);

                if ismember(cName, sessConds)
                    lapActivity = thisSess.ConditionData.(cName).LapActivity(idx, :, :);
                    if p.Results.ApplySmoothing, lapActivity = smoothLapActivity(lapActivity); end
                    hasData = true;
                else
                    lapActivity = nan(numel(idx), nLaps, nBins);
                end

                nTotalLaps   = size(lapActivity, 2);
                meanOddVals  = squeeze(mean(lapActivity(:, 1:2:nTotalLaps, :), 2, 'omitnan'));
                meanEvenVals = squeeze(mean(lapActivity(:, 2:2:nTotalLaps, :), 2, 'omitnan'));
                if size(meanOddVals,2) == 1, meanOddVals = meanOddVals'; meanEvenVals = meanEvenVals'; end

                dayMatrixOdd  = vertcat(dayMatrixOdd,  meanOddVals); 
                dayMatrixEven = vertcat(dayMatrixEven, meanEvenVals); 
            end

            if isBase
                qEnd = length(condDataStore) + 1;
                condDataStore(qEnd).rawMatrices    = dayMatrixOdd;
                condDataStore(qEnd).conditionNames = 'Baseline Odd';
                condDataStore(qEnd).condTypes      = cTypeOdd;
                condDataStore(qEnd).titleColors    = colColorOdd;
                condDataStore(qEnd).hasData        = true;

                qEnd = length(condDataStore) + 1;
                condDataStore(qEnd).rawMatrices    = dayMatrixEven;
                condDataStore(qEnd).conditionNames = 'Baseline Even';
                condDataStore(qEnd).condTypes      = cTypeEven;
                condDataStore(qEnd).titleColors    = colColorEven;
                condDataStore(qEnd).hasData        = true;
            else
                combinedLaps = mean(cat(3, dayMatrixOdd, dayMatrixEven), 3, 'omitnan');
                qEnd = length(condDataStore) + 1;
                condDataStore(qEnd).rawMatrices    = combinedLaps;
                condDataStore(qEnd).conditionNames = strrep(cName, '_', ' ');
                condDataStore(qEnd).condTypes      = cTypeOdd;
                condDataStore(qEnd).titleColors    = colColorOdd;
                condDataStore(qEnd).hasData        = hasData;
            end
        end
    end

    if isempty(condDataStore), figHandle = []; disp('No sessions matched.'); return; end

    %% Sorting and normalisation
    fullBaselineOdd  = condDataStore(1).rawMatrices;
    minFullOdd       = min(fullBaselineOdd, [], 2, 'omitnan');
    maxFullOdd       = max(fullBaselineOdd, [], 2, 'omitnan');
    rangeFullOdd     = maxFullOdd - minFullOdd; rangeFullOdd(rangeFullOdd == 0) = 1;

    normFullOddRef   = (fullBaselineOdd - minFullOdd) ./ rangeFullOdd;
    normFullOddRef(isnan(normFullOddRef)) = 0;
    [~, fullPeaks]   = max(normFullOddRef, [], 2);
    [~, fullSortIdx] = sort(fullPeaks);
    totalSomas       = size(fullBaselineOdd, 1);

    topRows           = fullSortIdx(max(4300,1):totalSomas);
    candidateIndices  = topRows;
    meanVal           = mean(fullBaselineOdd(candidateIndices,:), 2, 'omitnan');
    maxVal            = max(fullBaselineOdd(candidateIndices,:),  [], 2, 'omitnan');
    modulationDepth   = meanVal ./ maxVal;
    isTonicBackground = (modulationDepth > 0.40);
    isolatedGlobalIndices = candidateIndices(isTonicBackground);

    if isempty(isolatedGlobalIndices)
        error('Sparsity thresholds too aggressive — no background cells survived.');
    end

    minIsolateOdd   = minFullOdd(isolatedGlobalIndices);
    maxIsolateOdd   = maxFullOdd(isolatedGlobalIndices);
    rangeIsolateOdd = maxIsolateOdd - minIsolateOdd; rangeIsolateOdd(rangeIsolateOdd == 0) = 1;

    isolateBaselineOdd = fullBaselineOdd(isolatedGlobalIndices, :);
    smoothedForPeak    = movmean(isolateBaselineOdd, 3, 2);
    [~, peakBins]      = max(smoothedForPeak, [], 2);
    [~, localSortIdx]  = sort(peakBins, 'ascend');

    numN        = length(isolatedGlobalIndices);
    numBins     = size(isolateBaselineOdd, 2);
    cmPositions = linspace(0, 200, numBins);

    baseEvenIdx  = find(contains(lower({condDataStore.conditionNames}), 'baseline even'), 1);
    rawBaseEven  = condDataStore(baseEvenIdx).rawMatrices(isolatedGlobalIndices, :);
    normBaseEven = (rawBaseEven - minIsolateOdd) ./ rangeIsolateOdd;
    normBaseEven(isnan(normBaseEven)) = 0;

    %% Layout
    nCondCols   = 5;
    nRow1Cols   = 6;
    leftMargin  = 0.07;
    rightMargin = 0.10;
    colGap      = 0.01;
    cbW         = 0.015;
    cbGap       = 0.01;
    availW      = 1 - leftMargin - rightMargin - (nRow1Cols-1)*colGap;
    colW        = availW / nRow1Cols;
    row1Bot     = 0.68; row1H = 0.25;
    row2Bot     = 0.38; row2H = 0.25;
    row3Bot     = 0.07; row3H = 0.20;
    colLefts    = leftMargin + (0:nRow1Cols-1) * (colW + colGap);
    cbLeft      = 1 - rightMargin + cbGap;

    cMapLen = 256;
    b = [linspace(1,1,cMapLen/2), linspace(1,0,cMapLen/2)];
    g = [linspace(0,1,cMapLen/2), linspace(1,0,cMapLen/2)];
    r = [linspace(0,1,cMapLen/2), linspace(1,1,cMapLen/2)];
    redWhiteBlueMap = [b', g', r'];

    figHandle = figure('Position', [50 50 1400 900], 'Color', 'w');

    row23_Tokens   = {'swap.*2.*3', 'swap.*3.*4', 'omit.*2', 'omit.*3', 'omit.*4'};
    displayNames   = {'Swap 2 3', 'Swap 3 4', 'Omit 2', 'Omit 3', 'Omit 4'};
    condColors     = {coolColors(1,:), coolColors(2,:), warmColors(1,:), warmColors(2,:), warmColors(3,:)};
    landmarkTarget = {[80 120], [120 160], 80, 120, 160};

    %% ROW 1 — Absolute heatmaps
    row1_Indices = [2, 3, 4, 5, 6, 7];
    row1_Titles  = {'Base Even', 'Swap 2 3', 'Swap 3 4', 'Omit 2', 'Omit 3', 'Omit 4'};
    row1_Colors  = {[0 0 0], coolColors(1,:), coolColors(2,:), warmColors(1,:), warmColors(2,:), warmColors(3,:)};

    for i = 1:nRow1Cols
        axAbs    = axes('Position', [colLefts(i), row1Bot, colW, row1H]); 
        storeIdx = row1_Indices(i);

        if storeIdx > length(condDataStore) || ~condDataStore(storeIdx).hasData
            text(axAbs, 0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', ...
                'FontName', targetFont, 'FontSize', 12, 'Color', [0.5 0.5 0.5]);
            title(axAbs, row1_Titles{i}, 'Color', row1_Colors{i}, 'FontName', targetFont, 'FontSize', 11, 'FontWeight', 'bold');
            axis(axAbs, 'off'); continue;
        end

        localMatrix       = condDataStore(storeIdx).rawMatrices(isolatedGlobalIndices, :);
        normalizedDisplay = (localMatrix - minIsolateOdd) ./ rangeIsolateOdd;
        normalizedDisplay(isnan(normalizedDisplay)) = 0;

        imagesc(axAbs, [0 200], [1 numN], normalizedDisplay(localSortIdx,:));
        colormap(axAbs, flipud(gray));
        set(axAbs, 'CLim', [0.25 0.75], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
        set(axAbs, 'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'}, 'FontName', targetFont, 'FontSize', 10);
        title(axAbs, row1_Titles{i}, 'Color', row1_Colors{i}, 'FontName', targetFont, 'FontSize', 11, 'FontWeight', 'bold');

        if i == 1
            ylabel(axAbs, sprintf('Sorted %s\n(n=%d)', cellType, numN), 'FontName', targetFont, 'FontSize', 11);
        else
            set(axAbs, 'YColor', 'none', 'YTick', []);
        end

        hold(axAbs, 'on');
        for lm = [40 80 120 160]
            xline(axAbs, lm, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        end
    end

    % Colorbar row 1
    axAbsCB = axes('Position', [cbLeft, row1Bot, cbW, row1H]);
    set(axAbsCB, 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
    colormap(axAbsCB, flipud(gray));
    cbAbs = colorbar(axAbsCB, 'Location', 'eastoutside');
    cbAbs.Position      = [cbLeft, row1Bot, cbW, row1H];
    cbAbs.Ticks         = [0.25 0.50 0.75]; cbAbs.TickLabels = {'0.25','0.50','0.75'};
    cbAbs.TickDirection = 'out'; cbAbs.Box = 'off';
    cbAbs.FontName      = targetFont; cbAbs.FontSize = 9;
    cbAbs.Label.String  = 'Activity (norm.)';
    cbAbs.Label.FontName = targetFont; cbAbs.Label.FontSize = 11;

    %% ROWS 2 & 3
    for i = 1:nCondCols
        matchIdx  = find(~cellfun(@isempty, regexpi({condDataStore.conditionNames}, row23_Tokens{i})), 1);
        thisColor = condColors{i};
        thisLeft  = colLefts(i+1);

        %% Row 2 — Delta heatmap
        axHM = axes('Position', [thisLeft, row2Bot, colW, row2H]); 

        if isempty(matchIdx) || ~condDataStore(matchIdx).hasData
            text(axHM, 0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', ...
                'FontName', targetFont, 'FontSize', 12, 'Color', [0.5 0.5 0.5]);
            title(axHM, ['\Delta ' displayNames{i}], 'Color', thisColor, 'FontName', targetFont, 'FontSize', 11, 'FontWeight', 'bold');
            axis(axHM, 'off');

            % Still draw empty row 3 panel
            axLP = axes('Position', [thisLeft, row3Bot, colW, row3H]); 
            text(axLP, 0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', ...
                'FontName', targetFont, 'FontSize', 12, 'Color', [0.5 0.5 0.5]);
            axis(axLP, 'off');
        else
            localMatrix = condDataStore(matchIdx).rawMatrices(isolatedGlobalIndices, :);
            normCond    = (localMatrix - minIsolateOdd) ./ rangeIsolateOdd;
            normCond(isnan(normCond)) = 0;
            diffMat     = normCond - normBaseEven;

            imagesc(axHM, [0 200], [1 numN], diffMat(localSortIdx,:));
            colormap(axHM, flipud(redWhiteBlueMap));
            set(axHM, 'CLim', [-0.4 0.4], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
            set(axHM, 'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'}, 'FontName', targetFont, 'FontSize', 10);
            title(axHM, ['\Delta ' displayNames{i}], 'Color', thisColor, 'FontName', targetFont, 'FontSize', 11, 'FontWeight', 'bold');

            if i == 1
                ylabel(axHM, sprintf('Sorted %s\n(n=%d)', cellType, numN), 'FontName', targetFont, 'FontSize', 11);
            else
                set(axHM, 'YColor', 'none', 'YTick', []);
            end

            hold(axHM, 'on');
            targets = landmarkTarget{i};
            for lm = [40 80 120 160]
                if ismember(lm, targets)
                    xline(axHM, lm, '--', 'Color', thisColor, 'LineWidth', 2.0);
                else
                    xline(axHM, lm, ':',  'Color', [0.6 0.6 0.6], 'LineWidth', 1);
                end
            end

            %% Row 3 — Delta mean trace
            axLP = axes('Position', [thisLeft, row3Bot, colW, row3H]); 
            hold(axLP, 'on');

            muDiff  = mean(diffMat, 1, 'omitnan');
            semDiff = std(diffMat, 0, 1, 'omitnan') / sqrt(numN);

            fill(axLP, [cmPositions, fliplr(cmPositions)], ...
                [muDiff+semDiff, fliplr(muDiff-semDiff)], ...
                thisColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            plot(axLP, cmPositions, muDiff, 'Color', thisColor, 'LineWidth', 1.5);
            yline(axLP, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);

            for lm = [40 80 120 160]
                if ismember(lm, targets)
                    xline(axLP, lm, '--', 'Color', thisColor, 'LineWidth', 1.5);
                else
                    xline(axLP, lm, ':',  'Color', [0.6 0.6 0.6], 'LineWidth', 1);
                end
            end

            set(axLP, 'Box', 'off', 'TickDir', 'out', 'XLim', [0 200], 'YLim', [-0.4 0.4], ...
                'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'}, ...
                'YTick', [-0.3 0 0.3], 'FontName', targetFont, 'FontSize', 10);
            xlabel(axLP, 'Position (cm)', 'FontName', targetFont, 'FontSize', 11);

            if i == 1
                ylabel(axLP, '\Delta Activity', 'FontName', targetFont, 'FontSize', 11);
            else
                set(axLP, 'YColor', 'none', 'YTick', []);
            end
        end
    end

    % Colorbar row 2
    axDiffCB = axes('Position', [cbLeft, row2Bot, cbW, row2H]); 
    set(axDiffCB, 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
    colormap(axDiffCB, flipud(redWhiteBlueMap));
    cbDiff = colorbar(axDiffCB, 'Location', 'eastoutside');
    cbDiff.Position      = [cbLeft, row2Bot, cbW, row2H];
    cbDiff.Ticks         = [-0.4 0 0.4];
    cbDiff.TickLabels    = {'Dec.', '0', 'Inc.'};
    cbDiff.TickDirection = 'out'; cbDiff.Box = 'off';
    cbDiff.FontName      = targetFont; cbDiff.FontSize = 9;
    cbDiff.Label.String  = '\Delta Activity (Cond - Base Even)';
    cbDiff.Label.FontName = targetFont; cbDiff.Label.FontSize = 11;

    if ~isempty(p.Results.SavePath)
        saveFigureFormats(figHandle, p.Results.SavePath);
    end
end

% function figHandle = plotBackgroundROIs_DifferenceOnly(allData, targetArea, varargin)
% % plotBackgroundROIs_DifferenceOnly: Ultimate Combined Publication Layout Engine.
% % Row 1: Absolute Normalized Heatmaps (Grayscale, 4 Columns)
% % Row 2: Relative Delta Change Heatmaps (Red/White/Blue, 3 Columns)
% % Row 3: Single Y-Axis Delta Mean Line Graphs (3 Columns)
% 
%     p = inputParser;
%     addRequired(p, 'allData', @isstruct);
%     addRequired(p, 'targetArea', @(x) ischar(x) || isstring(x));
%     
%     addParameter(p, 'DaysToPlot', [1, 2, 3, 4, 5], @isnumeric);
%     addParameter(p, 'TypeToPlot', 'Boutons', @(x) ischar(x) || isstring(x)); 
%     addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x)); 
%     addParameter(p, 'ApplySmoothing', true, @islogical);
%     addParameter(p, 'FontName', 'Arial', @(x) ischar(x) || isstring(x)); 
%     
%     parse(p, allData, targetArea, varargin{:});
%     targetFont = p.Results.FontName;
%     cellType = p.Results.TypeToPlot;
%     
%     warmColors = [0.541, 0.012, 0.012; 0.824, 0.016, 0.176];  
%     coolColors = [0.000, 0.400, 1.000];
%     
%     if all(cellfun(@isempty, {allData.TargetArea})), [allData.TargetArea] = deal(char(targetArea)); end
%     if isfield(allData, 'Type') && ~isfield(allData, 'TypeImaged')
%         [allData.TypeImaged] = allData.Type;
%     elseif all(cellfun(@isempty, {allData.TypeImaged}))
%         [allData.TypeImaged] = deal(char(cellType));
%     end
%     
%     daysToPlot = p.Results.DaysToPlot;
%     daysToPlot(daysToPlot == 200) = []; 
%     nDays = length(daysToPlot);
%     
%     condDataStore = struct('rawMatrices', {}, 'conditionNames', {}, 'condTypes', {}, 'titleColors', {});
%     
%     for d = 1:nDays
%         day = daysToPlot(d);
%         if day == 5
%             dayMask = ([allData.Day] == 5 | [allData.Day] == 200);
%         else
%             dayMask = ([allData.Day] == day);
%         end
%         
%         rawSessions = allData(dayMask & ...
%                               strcmpi(string({allData.TargetArea}), string(targetArea)) & ...
%                               strcmpi(string({allData.TypeImaged}), string(cellType)));
%         
%         if isempty(rawSessions), continue; end
%         
%         allCondsInDay = {};
%         for s = 1:length(rawSessions)
%             if isfield(rawSessions(s), 'ConditionData') && ~isempty(rawSessions(s).ConditionData)
%                 allCondsInDay = unique([allCondsInDay; fieldnames(rawSessions(s).ConditionData)]);
%             end
%         end
%         baseIdx = find(contains(lower(allCondsInDay), 'baseline') | contains(lower(allCondsInDay), 'default'), 1);
%         if isempty(baseIdx), baseIdx = 1; end
%         baseName = allCondsInDay{baseIdx};
%         
%         targetSwaps = {'Swap_2_3', 'Swap_3_4'};
%         targetOmits = {'Omit_2', 'Omit_3', 'Omit_4'}; 
%         requiredConditions = [{baseName}, targetSwaps, targetOmits];
%         
%         validSessionMask = false(1, length(rawSessions));
%         for s = 1:length(rawSessions)
%             if isfield(rawSessions(s), 'ConditionData') && ~isempty(rawSessions(s).ConditionData)
%                 sessConds = fieldnames(rawSessions(s).ConditionData);
%                 if all(ismember(requiredConditions, sessConds))
%                     validSessionMask(s) = true;
%                 end
%             end
%         end
%         daySessions = rawSessions(validSessionMask);
%         if isempty(daySessions), continue; end
%         
%         orderedConds = requiredConditions;
%         omitCount = 1; swapCount = 1;
%         
%         for c = 1:length(orderedConds)
%             cName = orderedConds{c};
%             isBaseBlock = (c == 1);
%             nameLow = lower(cName);
%             
%             if isBaseBlock
%                 cTypeOdd = 'baseline_odd'; cTypeEven = 'baseline_even';
%                 colColorOdd = [0 0 0]; colColorEven = [0 0 0];
%             elseif contains(nameLow, 'swap')
%                 cTypeOdd = 'swap'; cTypeEven = 'swap';
%                 colColorOdd = coolColors(mod(swapCount-1, size(coolColors,1))+1, :);
%                 colColorEven = colColorOdd; swapCount = swapCount + 1;
%             elseif contains(nameLow, 'omit')
%                 cTypeOdd = 'omit'; cTypeEven = 'omit';
%                 colColorOdd = warmColors(mod(omitCount-1, size(warmColors,1))+1, :);
%                 colColorEven = colColorOdd; omitCount = omitCount + 1;
%             end
%             
%             dayMatrixOdd = []; dayMatrixEven = [];
%             
%             for s = 1:length(daySessions)
%                 thisSess = daySessions(s);
%                 if isfield(thisSess, 'FilteredROIs') && ~isempty(thisSess.FilteredROIs)
%                     idx = thisSess.FilteredROIs;
%                 else
%                     idx = 1:size(thisSess.ConditionData.(baseName).LapActivity, 1);
%                 end
%                 
%                 lapActivity = thisSess.ConditionData.(cName).LapActivity(idx, :, :);
%                 if p.Results.ApplySmoothing, lapActivity = smoothLapActivity(lapActivity); end
%                 
%                 nTotalLaps = size(lapActivity, 2);
%                 
%                 meanOddVals = squeeze(mean(lapActivity(:, 1:2:nTotalLaps, :), 2, 'omitnan'));
%                 meanEvenVals = squeeze(mean(lapActivity(:, 2:2:nTotalLaps, :), 2, 'omitnan'));
%                 
%                 if size(meanOddVals, 2) == 1, meanOddVals = meanOddVals'; meanEvenVals = meanEvenVals'; end
%                 
%                 dayMatrixOdd  = vertcat(dayMatrixOdd, meanOddVals); %#ok<AGROW>
%                 dayMatrixEven = vertcat(dayMatrixEven, meanEvenVals); %#ok<AGROW>
%             end
%             
%             if isBaseBlock
%                 qEnd = length(condDataStore) + 1;
%                 condDataStore(qEnd).rawMatrices = dayMatrixOdd;
%                 condDataStore(qEnd).conditionNames = 'Baseline Odd';
%                 condDataStore(qEnd).condTypes = cTypeOdd;
%                 condDataStore(qEnd).titleColors = colColorOdd;
%                 
%                 qEnd = length(condDataStore) + 1;
%                 condDataStore(qEnd).rawMatrices = dayMatrixEven;
%                 condDataStore(qEnd).conditionNames = 'Baseline Even';
%                 condDataStore(qEnd).condTypes = cTypeEven;
%                 condDataStore(qEnd).titleColors = colColorEven;
%             else
%                 combinedOtherLaps = mean(cat(3, dayMatrixOdd, dayMatrixEven), 3, 'omitnan');
%                 qEnd = length(condDataStore) + 1;
%                 condDataStore(qEnd).rawMatrices = combinedOtherLaps;
%                 condDataStore(qEnd).conditionNames = strrep(cName, '_', ' ');
%                 condDataStore(qEnd).condTypes = cTypeOdd;
%                 condDataStore(qEnd).titleColors = colColorOdd;
%             end
%         end
%     end
%     
%     if isempty(condDataStore), figHandle = []; disp('No sessions matched.'); return; end
%     
%     %% --- MASTER SORTING & SUBSET ISOLATION ---
%     fullBaselineOdd = condDataStore(1).rawMatrices;
%     minFullOdd = min(fullBaselineOdd, [], 2, 'omitnan'); 
%     maxFullOdd = max(fullBaselineOdd, [], 2, 'omitnan'); 
%     rangeFullOdd = maxFullOdd - minFullOdd; rangeFullOdd(rangeFullOdd == 0) = 1;
%     
%     normFullOddRef = (fullBaselineOdd - minFullOdd) ./ rangeFullOdd;
%     normFullOddRef(isnan(normFullOddRef)) = 0;
%     [~, fullPeaks] = max(normFullOddRef, [], 2);
%     [~, fullSortIdx] = sort(fullPeaks);
%     totalSomas = size(fullBaselineOdd, 1);
%     
%     topRows = fullSortIdx(max(4300, 1) : totalSomas);%1800
%     candidateIndices = topRows;
%     
%     meanVal = mean(fullBaselineOdd(candidateIndices, :), 2, 'omitnan');
%     maxVal  = max(fullBaselineOdd(candidateIndices, :), [], 2, 'omitnan');
%     modulationDepth = meanVal ./ maxVal; 
%     
%     isTonicBackground = (modulationDepth > 0.40); 
%     isolatedGlobalIndices = candidateIndices(isTonicBackground);
%     
%     if isempty(isolatedGlobalIndices)
%         error('Sparsity thresholds are too aggressive. No un-tuned cells survived.');
%     end
%     
%     minIsolateOdd = minFullOdd(isolatedGlobalIndices);
%     maxIsolateOdd = maxFullOdd(isolatedGlobalIndices);
%     rangeIsolateOdd = maxIsolateOdd - minIsolateOdd; 
%     rangeIsolateOdd(rangeIsolateOdd == 0) = 1;
%     
%     isolateBaselineOdd = fullBaselineOdd(isolatedGlobalIndices, :);
%     smoothedForPeak = movmean(isolateBaselineOdd, 3, 2); 
%     [~, peakBins] = max(smoothedForPeak, [], 2);
%     [~, localSortIdx] = sort(peakBins, 'ascend');
%     
%     numN = length(isolatedGlobalIndices);
%     numBins = size(isolateBaselineOdd, 2);
%     cmPositions = linspace(0, 200, numBins);
%     
%     %% --- EXTRACT BASELINE ARRAYS FOR EXTRACTION ---
%     baseEvenIdx = find(contains(lower({condDataStore.conditionNames}), 'baseline even'), 1);
%     rawBaseEven = condDataStore(baseEvenIdx).rawMatrices(isolatedGlobalIndices, :);
%     normBaseEven = (rawBaseEven - minIsolateOdd) ./ rangeIsolateOdd;
%     normBaseEven(isnan(normBaseEven)) = 0;
%     
%     %% --- EXPLICIT GEOMETRIC LAYOUT GRID MAPPING ---
%     % Custom footprints structurally space rows 1, 2, and 3 cleanly.
%     % Format: [Left, Bottom, Width, Height]
%     row1_Positions = {[0.10, 0.70, 0.14, 0.23], [0.28, 0.70, 0.14, 0.23], [0.46, 0.70, 0.14, 0.23], [0.64, 0.70, 0.14, 0.23]}; % Absolute (4 Panels)
%     row2_Positions = {[0.28, 0.38, 0.14, 0.23], [0.46, 0.38, 0.14, 0.23], [0.64, 0.38, 0.14, 0.23]};                         % Delta (3 Panels)
%     row3_Positions = {[0.28, 0.10, 0.14, 0.16], [0.46, 0.10, 0.14, 0.16], [0.64, 0.10, 0.14, 0.16]};                         % Curves (3 Panels)
%     
%     figHandle = figure('Position', [50 50 1200 900], 'Color', 'w');
%     
%     % Symmetrical Divergent Red-White-Blue Map Setup
%     cMapLen = 256;
%     b = [linspace(1, 1, cMapLen/2), linspace(1, 0, cMapLen/2)];
%     g = [linspace(0, 1, cMapLen/2), linspace(1, 0, cMapLen/2)];
%     r = [linspace(0, 1, cMapLen/2), linspace(1, 1, cMapLen/2)];
%     redWhiteBlueMap = [b', g', r']; 
%     
%     row1_Indices = [2, 3, 4, 5]; % Base Even, Swap 2 3, Omit 2, Omit 3
%     row23_Tokens = {'swap', 'omit.*2', 'omit.*3'};
%     displayNames = {'Swap 2 3', 'Omit 2', 'Omit 3'};
%     colors = {coolColors(1,:), warmColors(1,:), warmColors(2,:)};
%     
%     %% ====================================================================
%     %% --- ROW 1: ABSOLUTE NORMALIZED HEATMAP PANELS (GRAYSCALE) ---
%     %% ====================================================================
%     for i = 1:4
%         axAbs = axes('Position', row1_Positions{i}); %#ok<LAXES>
%         currentCond = condDataStore(row1_Indices(i));
%         
%         localIsolateMatrix = currentCond.rawMatrices(isolatedGlobalIndices, :);
%         normalizedDisplay = (localIsolateMatrix - minIsolateOdd) ./ rangeIsolateOdd;
%         normalizedDisplay(isnan(normalizedDisplay)) = 0;
%         
%         imagesc(axAbs, [0, 200], [1, numN], normalizedDisplay(localSortIdx, :));
%         colormap(axAbs, flipud(gray));
%         set(axAbs, 'CLim', [0.25 0.75], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
%         set(axAbs, 'XTick', [40 80 120 160], 'XTickLabel', {'40', '80', '120', '160'}, 'FontName', targetFont, 'FontSize', 10);
%         
%         if i == 1
%             title(axAbs, 'Base Even', 'Color', [0 0 0], 'FontSize', 12, 'FontName', targetFont, 'FontWeight', 'bold');
%             ylabel(axAbs, sprintf('Sorted Background %s\n(n = %d)', cellType, numN), 'FontName', targetFont, 'FontSize', 11);
%             text(axAbs, -15, numN, sprintf('%d background %s', numN, lower(cellType)), ...
%                 'Rotation', 90, 'FontName', targetFont, 'FontSize', 11, 'FontAngle', 'italic', ...
%                 'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Clipping', 'off');
%         else
%             set(axAbs, 'YColor', 'none', 'YTick', []);
%             % Parse exact titles for manipulation tracks
%             cleanTitle = strrep(currentCond.conditionNames, '_', ' ');
%             title(axAbs, cleanTitle, 'Color', currentCond.titleColors, 'FontSize', 12, 'FontName', targetFont, 'FontWeight', 'bold');
%         end
%         
%         % Plot static gray landmark line trackers
%         hold(axAbs, 'on');
%         for posBin = [40 80 120 160]
%             xline(axAbs, posBin, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
%         end
%     end
%     
%     % Isolated Grayscale Colorbar Generation
%     axAbsCB = axes('Position', [0.79, 0.70, 0.015, 0.23]); 
%     set(axAbsCB, 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
%     colormap(axAbsCB, flipud(gray));
%     cbAbs = colorbar(axAbsCB, 'Location', 'eastoutside');
%     cbAbs.Position = [0.79, 0.70, 0.015, 0.23];
%     cbAbs.Ticks = [0.25 0.50 0.75]; cbAbs.TickLabels = {'0.25', '0.50', '0.75'};
%     cbAbs.TickDirection = 'out'; cbAbs.Box = 'off'; cbAbs.FontName = targetFont; cbAbs.FontSize = 9;
%     cbAbs.Label.String = 'Activity (norm.)'; cbAbs.Label.FontName = targetFont; cbAbs.Label.FontSize = 11;
% 
%     %% ====================================================================
%     %% --- ROW 2 & ROW 3: DELTA HEATMAPS & SINGLE-AXIS DELTA MEAN LINES ---
%     %% ====================================================================
%     for i = 1:3
%         matchIdx = find(~cellfun(@isempty, regexpi({condDataStore.conditionNames}, row23_Tokens{i})), 1);
%         
%         if isempty(matchIdx), continue; end
%         currentCond = condDataStore(matchIdx);
%         
%         localIsolateMatrix = currentCond.rawMatrices(isolatedGlobalIndices, :);
%         normCond = (localIsolateMatrix - minIsolateOdd) ./ rangeIsolateOdd;
%         normCond(isnan(normCond)) = 0; 
%         diffMat = normCond - normBaseEven;
%         
%         %% --- ROW 2: DELTA CHANGE HEATMAPS ---
%         axHM = axes('Position', row2_Positions{i}); %#ok<LAXES>
%         imagesc(axHM, [0, 200], [1, numN], diffMat(localSortIdx, :));
%         colormap(axHM, flipud(redWhiteBlueMap));
%         set(axHM, 'CLim', [-0.4 0.4], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
%         set(axHM, 'XTick', [40 80 120 160], 'XTickLabel', {'40', '80', '120', '160'}, 'FontName', targetFont, 'FontSize', 10);
%         title(axHM, ['\Delta ' displayNames{i}], 'Color', colors{i}, 'FontSize', 12, 'FontName', targetFont, 'FontWeight', 'bold');
%         
%         % Map condition milestone lines dynamically
%         landmarks = [40 80 120 160];
%         currentNameLow = lower(currentCond.conditionNames);
%         hold(axHM, 'on');
%         for lIdx = 1:4
%             posBin = landmarks(lIdx);
%             isTargetLandmark = false; 
%             if contains(currentNameLow, num2str(lIdx)), isTargetLandmark = true; end
%             
%             if isTargetLandmark
%                 xline(axHM, posBin, '--', 'Color', colors{i}, 'LineWidth', 2.0);
%             else
%                 xline(axHM, posBin, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
%             end
%         end
%         
%         if i == 1
%             ylabel(axHM, sprintf('Sorted Background %s\n(n = %d)', cellType, numN), 'FontName', targetFont, 'FontSize', 11);
%         else
%             set(axHM, 'YColor', 'none', 'YTick', []);
%         end
%         
%         %% --- ROW 3: COMPACT SINGLE Y-AXIS DELTA MEAN PROFILES ---
%         axLP = axes('Position', row3_Positions{i}); %#ok<LAXES>
%         hold(axLP, 'on');
%         
%         muDiff = mean(diffMat, 1, 'omitnan');
%         semDiff = std(diffMat, 0, 1, 'omitnan') / sqrt(numN);
%         fill(axLP, [cmPositions, fliplr(cmPositions)], [muDiff + semDiff, fliplr(muDiff - semDiff)], ...
%              colors{i}, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
%         plot(axLP, cmPositions, muDiff, 'Color', colors{i}, 'LineWidth', 1.5, 'DisplayName', '\Delta Mean');
%         
%         % Symmetric flat dash zero baseline reference line
%         yline(axLP, 0, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.2, 'DisplayName', 'No Change');
%         
%         set(axLP, 'Box', 'off', 'TickDir', 'out', 'XLim', [0 200], 'YLim', [-0.4 0.4], ...
%                 'XTick', [40 80 120 160], 'XTickLabel', {'40', '80', '120', '160'}, ...
%                 'YTick', [-0.3 0 0.3], 'FontName', targetFont, 'FontSize', 10);
%         xlabel(axLP, 'Position (cm)', 'FontName', targetFont, 'FontSize', 11);
%         
%         legend(axLP, 'show', 'Orientation', 'horizontal', 'Location', 'northoutside', 'Box', 'off', ...
%                'FontSize', 8, 'FontName', targetFont);
%         
%         if i == 1
%             ylabel(axLP, '\Delta Activity (Cond - Base)', 'FontName', targetFont, 'FontSize', 11);
%         else
%             set(axLP, 'YColor', 'none', 'YTick', []);
%         end
%     end
%     
%     %% --- SECTOR 6: RE-ALIGNED DELTA COLORBAR FLANK ---
%     axCB = axes('Position', [0.79, 0.38, 0.015, 0.23]); 
%     set(axCB, 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
%     colormap(axCB, flipud(redWhiteBlueMap));
%     
%     cbDiff = colorbar(axCB, 'Location', 'eastoutside');
%     cbDiff.Position = [0.79, 0.38, 0.015, 0.23];
%     cbDiff.Ticks = [-0.4 0 0.4]; cbDiff.TickLabels = {'Decrease (-)', 'No Change (0)', 'Increase (+)'};
%     cbDiff.TickDirection = 'out'; cbDiff.Box = 'off'; cbDiff.FontName = targetFont; cbDiff.FontSize = 9;
%     cbDiff.Label.String = '\Delta Activity (Condition - Baseline Even)';
%     cbDiff.Label.FontName = targetFont; cbDiff.Label.FontSize = 11;
%     
%     if ~isempty(p.Results.SavePath), saveFigureFormats(figHandle, p.Results.SavePath); end
% end