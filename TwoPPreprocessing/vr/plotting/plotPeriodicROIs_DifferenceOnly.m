function figHandle = plotPeriodicROIs_DifferenceOnly(allData, targetArea, varargin)
% plotPeriodicROIs_DifferenceOnly: Same heatmap/delta/mean-trace layout as
% plotBackgroundROIs_DifferenceOnly, but the cell-isolation step now finds
% PERIODICALLY-RESPONDING cells (e.g. ~7-8 repeating fields across the
% track) instead of tonic background cells.
%
% Periodicity is quantified per cell from its baseline (odd-lap) tuning
% curve using the fraction of spectral power falling within a
% [MinCycles, MaxCycles] band (cycles across the full track length).
% Significance is assessed against a null built by independently
% permuting each cell's own position bins (NShuffles times), which
% destroys periodic structure while preserving the value distribution.
%
% Cells are then sorted for display by the PHASE of their tuning curve at
% a shared reference spatial frequency (rather than by single-peak
% position), which aligns repeating fields into clean diagonal stripes.

    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    addRequired(p, 'targetArea', @(x) ischar(x) || isstring(x));
    addParameter(p, 'DaysToPlot', [1, 2, 3, 4, 5], @isnumeric);
    addParameter(p, 'TypeToPlot', 'Boutons', @(x) ischar(x) || isstring(x));
    addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ApplySmoothing', true, @islogical);
    addParameter(p, 'FontName', 'Arial', @(x) ischar(x) || isstring(x));
    addParameter(p, 'RequiredConditions', {}, @iscell);
    addParameter(p, 'MinFieldSpacingCm', 20, @isnumeric); % lower bound of candidate field-spacing search (~true spacing 24 cm)
    addParameter(p, 'MaxFieldSpacingCm', 28, @isnumeric); % upper bound of candidate field-spacing search (~true spacing 24 cm)
    addParameter(p, 'NHarmonics', 8, @isnumeric);         % number of harmonics (P, 2P, ..., 8P) averaged into the score — 200cm/24cm ≈ 8 repeats expected
    addParameter(p, 'MinOverlapFraction', 0.15, @isnumeric); % a harmonic lag is only used if at least this fraction of the track still overlaps (avoids noisy near-edge estimates)
    addParameter(p, 'NShuffles', 200, @isnumeric);        % permutations for the null distribution
    addParameter(p, 'PeriodicAlpha', 0.01, @isnumeric);   % significance threshold (uncorrected)
    addParameter(p, 'LandmarkPositions', [40 80 120 160], @isnumeric); % known landmark positions (cm) to exclude
    addParameter(p, 'LandmarkToleranceCm', 8, @isnumeric);             % how close (cm) a field must be to count as landmark-aligned
    addParameter(p, 'MaxLandmarkFraction', 0.5, @isnumeric);           % exclude cell if more than this fraction of its fields are landmark-aligned
    parse(p, allData, targetArea, varargin{:});

    targetFont   = p.Results.FontName;
    cellType     = p.Results.TypeToPlot;
    minSpacingCm = p.Results.MinFieldSpacingCm;
    maxSpacingCm = p.Results.MaxFieldSpacingCm;
    nHarmonics   = p.Results.NHarmonics;
    minOverlapFrac = p.Results.MinOverlapFraction;
    nShuffles    = p.Results.NShuffles;
    alphaThr     = p.Results.PeriodicAlpha;
    landmarkPositions = p.Results.LandmarkPositions;
    landmarkTolCm     = p.Results.LandmarkToleranceCm;
    maxLandmarkFrac   = p.Results.MaxLandmarkFraction;

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

        requiredForSelection = p.Results.RequiredConditions;
        if isempty(requiredForSelection)
            requiredForSelection = {baseName};
        else
            requiredForSelection = unique([{baseName}, requiredForSelection(:)']);
        end

        fprintf('  Requiring conditions: %s\n', strjoin(requiredForSelection, ', '));

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

    %% Sorting and normalisation (identical to background version)
    fullBaselineOdd  = condDataStore(1).rawMatrices;
    minFullOdd       = min(fullBaselineOdd, [], 2, 'omitnan');
    maxFullOdd       = max(fullBaselineOdd, [], 2, 'omitnan');
    rangeFullOdd     = maxFullOdd - minFullOdd; rangeFullOdd(rangeFullOdd == 0) = 1;

    normFullOddRef   = (fullBaselineOdd - minFullOdd) ./ rangeFullOdd;
    normFullOddRef(isnan(normFullOddRef)) = 0;
    totalSomas       = size(fullBaselineOdd, 1);

    %% --- PERIODICITY-BASED CELL ISOLATION (replaces tonic-background isolation) ---
    % Uses linear autocorrelation of each cell's baseline tuning curve. For
    % each cell, the candidate field spacing P is the lag (within
    % [MinFieldSpacingCm, MaxFieldSpacingCm]) with maximal autocorrelation.
    % The cell's score is the mean autocorrelation at P, 2P, 3P, ... (up to
    % NHarmonics, capped by track length) - a genuinely periodic cell keeps
    % recurring at multiples of its spacing, whereas a cell with only 1-2
    % real fields (e.g. landmark-tuned) will not. Significance is assessed
    % against a null built by independently permuting each cell's own bins
    % (autocorrelation, like the FFT, is shift-invariant, so plain circular
    % shifts would NOT work as a null here).
    fprintf('  Testing %d cells for periodicity via autocorrelation (spacing %g-%g cm, %d harmonics, %d shuffles)...\n', ...
        totalSomas, minSpacingCm, maxSpacingCm, nHarmonics, nShuffles);

    numBinsFull = size(fullBaselineOdd, 2);
    cmPerBin    = 200 / numBinsFull;
    minLagBins  = max(1, round(minSpacingCm / cmPerBin));
    maxLagBins  = max(minLagBins, round(maxSpacingCm / cmPerBin));
    candidateRange = minLagBins:maxLagBins;

    traceMat = movmean(normFullOddRef, 3, 2);
    minOverlapBins = round(minOverlapFrac * numBinsFull);

    acfMat = computeAutocorrBatch(traceMat);
    [~, bestRelIdx] = max(acfMat(:, candidateRange), [], 2);
    bestLagBins     = candidateRange(bestRelIdx)';
    harmonicScore   = computeHarmonicScore(acfMat, bestLagBins, nHarmonics, minOverlapBins);

    nullScore = nan(totalSomas, nShuffles);
    for sh = 1:nShuffles
        shuffledMat = permuteRowsIndependently(traceMat);
        acfShuffle  = computeAutocorrBatch(shuffledMat);
        [~, bestRelIdxSh] = max(acfShuffle(:, candidateRange), [], 2);
        bestLagBinsSh      = candidateRange(bestRelIdxSh)';
        nullScore(:, sh)   = computeHarmonicScore(acfShuffle, bestLagBinsSh, nHarmonics, minOverlapBins);
    end

    pValuePeriodic = (sum(nullScore >= harmonicScore, 2) + 1) / (nShuffles + 1);
    isPeriodic     = pValuePeriodic < alphaThr;
    autocorrPassIndices = find(isPeriodic);

    fprintf('  %d/%d cells passed the autocorrelation periodicity threshold (p < %.3g)\n', ...
        numel(autocorrPassIndices), totalSomas, alphaThr);

    %% --- LANDMARK-ALIGNMENT EXCLUSION ---
    % Autocorrelation alone can't fully separate "sustained fine-grained
    % periodicity" from "coarse landmark tuning" (e.g. 4 fields locked to
    % the 4 landmarks), since broad landmark responses can still leak some
    % correlational structure into the candidate spacing window. Here we
    % explicitly find each candidate cell's field locations and reject it
    % if most of its fields sit on top of a known landmark position rather
    % than being spread across the track.
    posCmFull = linspace(0, 200, numBinsFull);
    keepMask  = true(size(autocorrPassIndices));

    for ii = 1:numel(autocorrPassIndices)
        n = autocorrPassIndices(ii);
        [~, pkLocsBins] = findpeaks(traceMat(n, :), 'MinPeakProminence', 0.1);
        if isempty(pkLocsBins)
            keepMask(ii) = false;
            continue;
        end
        pkPositionsCm = posCmFull(pkLocsBins)';
        distToLandmark = min(abs(pkPositionsCm - landmarkPositions), [], 2);
        fracAligned = mean(distToLandmark <= landmarkTolCm);
        if fracAligned > maxLandmarkFrac
            keepMask(ii) = false;
        end
    end

    isolatedGlobalIndices = autocorrPassIndices(keepMask);

    fprintf('  %d/%d autocorrelation-passing cells kept after excluding landmark-aligned fields (>%.0f%% of fields within %g cm of a landmark)\n', ...
        numel(isolatedGlobalIndices), numel(autocorrPassIndices), maxLandmarkFrac*100, landmarkTolCm);

    if isempty(isolatedGlobalIndices)
        error('No periodic cells survived both the autocorrelation test and the landmark-exclusion filter — try relaxing PeriodicAlpha, MaxLandmarkFraction, or LandmarkToleranceCm.');
    end

    minIsolateOdd   = minFullOdd(isolatedGlobalIndices);
    maxIsolateOdd   = maxFullOdd(isolatedGlobalIndices);
    rangeIsolateOdd = maxIsolateOdd - minIsolateOdd; rangeIsolateOdd(rangeIsolateOdd == 0) = 1;

    isolateBaselineOdd = fullBaselineOdd(isolatedGlobalIndices, :);
    smoothedForPeak    = movmean(isolateBaselineOdd, 3, 2);

    % --- PHASE-BASED SORT (replaces peak-position sort for periodic cells) ---
    % Peak-position sorting is a poor match for periodic responses: noise
    % can make any one of several roughly-equal repeating fields appear
    % tallest, which scrambles cells that are otherwise well-aligned into a
    % noisy rather than diagonal ordering. Instead, a single representative
    % spatial period is taken as the median of each cell's own best-fit
    % field spacing (bestLagBins, from the autocorrelation test above), and
    % each cell's phase relative to this shared period is computed via a
    % discrete Fourier transform of its smoothed tuning curve at that
    % frequency. Sorting by this phase aligns repeating fields into clean
    % diagonal stripes, since phase reflects the position of the whole
    % repeating pattern rather than any single field's amplitude.
%     repPeriodBins    = median(bestLagBins(isolatedGlobalIndices));
%     freqCyclesPerBin = 1 / repPeriodBins;
% 
%     nColsIso = size(isolateBaselineOdd, 2);
%     binIdx   = (0:nColsIso-1);
%     complexPhase = zeros(length(isolatedGlobalIndices), 1);
%     for n = 1:length(isolatedGlobalIndices)
%         trace = smoothedForPeak(n, :);
%         coeff = sum(trace .* exp(-1i * 2*pi * freqCyclesPerBin * binIdx));
%         complexPhase(n) = angle(coeff);
%     end
% 
%     wrappedPhase = mod(complexPhase, 2*pi); % wrap to [0, 2*pi) for a continuous, non-negative ordering
%     [~, localSortIdx] = sort(wrappedPhase, 'descend');   % was 'ascend'

    repPeriodBins = median(bestLagBins(isolatedGlobalIndices));

    [~, peakBinRaw] = max(smoothedForPeak, [], 2);
    wrappedPhaseBins = mod(peakBinRaw - 1, repPeriodBins); % -1/+1 handles 1-based MATLAB indexing
    [~, localSortIdx] = sort(wrappedPhaseBins, 'ascend');

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
            ylabel(axAbs, sprintf('Sorted Periodic %s\n(n=%d)', cellType, numN), 'FontName', targetFont, 'FontSize', 11);
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
                ylabel(axHM, sprintf('Sorted Periodic %s\n(n=%d)', cellType, numN), 'FontName', targetFont, 'FontSize', 11);
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

            %% Row 3 — Clean mean delta trace (single axis, no dotted lines)
            axLP = axes('Position', [thisLeft, row3Bot, colW, row3H]); 
            hold(axLP, 'on');

            muDiff  = mean(diffMat, 1, 'omitnan');
            semDiff = std(diffMat, 0, 1, 'omitnan') / sqrt(numN);

            fill(axLP, [cmPositions, fliplr(cmPositions)], ...
                [muDiff+semDiff, fliplr(muDiff-semDiff)], ...
                thisColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            plot(axLP, cmPositions, muDiff, 'Color', thisColor, 'LineWidth', 1.5);
            yline(axLP, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);

            for lm = targets
                xline(axLP, lm, '--', 'Color', thisColor, 'LineWidth', 1.5);
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

%% ---- Local helper functions ----

function acfMat = computeAutocorrBatch(mat)
% Vectorized linear (non-circular) autocorrelation for every row of mat,
% via zero-padded FFT. Uses an UNBIASED normalization: each lag's raw sum
% is divided by its own number of overlapping bins (nCols - lag), not
% just by the lag-0 value. Without this correction, the estimate decays
% toward zero at long lags purely because fewer bins overlap there --
% not because periodicity is actually weaker -- which unfairly penalizes
% high-order harmonics near the edge of the track.
    [nRows, nCols] = size(mat);
    nfft = 2^nextpow2(2*nCols);
    X    = fft(mat, nfft, 2);
    S    = X .* conj(X);
    acfFull = real(ifft(S, [], 2));
    rawAcf  = acfFull(:, 1:nCols);

    overlapCounts = nCols - (0:nCols-1); % number of overlapping bins at each lag
    unbiasedAcf   = rawAcf ./ overlapCounts;

    lag0 = unbiasedAcf(:, 1);
    lag0(lag0 == 0) = eps;
    acfMat = unbiasedAcf ./ lag0;
end

function score = computeHarmonicScore(acfMat, bestLagBins, nHarmonics, minOverlapBins)
% For each row (cell), averages the autocorrelation value at lag =
% bestLagBins*k for k = 1..nHarmonics, skipping any harmonic whose lag
% exceeds the available number of columns (track length), OR whose
% remaining overlap is below minOverlapBins (too close to the edge to
% trust the estimate).
    [nRows, nCols] = size(acfMat);
    scoreSum    = zeros(nRows, 1);
    countValid  = zeros(nRows, 1);
    for k = 1:nHarmonics
        lagVal      = bestLagBins * k;
        colIdx      = lagVal + 1; % +1: lag 0 -> column 1
        overlapAtLag = nCols - lagVal;
        validMask   = (colIdx <= nCols) & (overlapAtLag >= minOverlapBins);
        rowsValid   = find(validMask);
        linIdx      = sub2ind([nRows, nCols], rowsValid, colIdx(rowsValid));
        scoreSum(rowsValid)   = scoreSum(rowsValid) + acfMat(linIdx);
        countValid(rowsValid) = countValid(rowsValid) + 1;
    end
    countValid(countValid == 0) = 1;
    score = scoreSum ./ countValid;
end

function out = permuteRowsIndependently(mat)
% Independently randomly permutes the columns of each row of mat.
% Used to build a null distribution that destroys periodic structure
% while preserving each cell's own value distribution.
    [nRows, nCols] = size(mat);
    [~, idx] = sort(rand(nRows, nCols), 2);
    linIdx   = sub2ind([nRows, nCols], repmat((1:nRows)', 1, nCols), idx);
    out      = mat(linIdx);
end