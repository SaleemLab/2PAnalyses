function figHandle = plotPeriodicROIs_BaselineOnly(allData, targetArea, varargin)
% plotPeriodicROIs_BaselineOnly: Detects PERIODICALLY-RESPONDING cells
% (e.g. ~7-8 repeating fields across the track) from baseline data only.
% This is a simplified version of plotPeriodicROIs_DifferenceOnly with
% the swap/omit delta-heatmap and delta-trace rows removed; it just
% identifies and displays periodic cells present in the baseline
% condition.
%
% Periodicity is quantified per cell from its baseline (odd-lap) tuning
% curve using autocorrelation at a candidate field spacing and its
% harmonics. Significance is assessed against a null built by
% independently permuting each cell's own position bins (NShuffles
% times), which destroys periodic structure while preserving the value
% distribution.
%
% Cells are sorted for display by the position of their peak within one
% repeat period (wrapped phase), which aligns repeating fields into
% clean diagonal stripes.

    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    addRequired(p, 'targetArea', @(x) ischar(x) || isstring(x));
    addParameter(p, 'DaysToPlot', [1, 2, 3, 4, 5], @isnumeric);
    addParameter(p, 'TypeToPlot', 'Boutons', @(x) ischar(x) || isstring(x));
    addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ApplySmoothing', true, @islogical);
    addParameter(p, 'FontName', 'Arial', @(x) ischar(x) || isstring(x));
    addParameter(p, 'RequiredConditions', {}, @iscell);
    addParameter(p, 'MinFieldSpacingCm', 24, @isnumeric); % lower bound of candidate field-spacing search (~true spacing 24 cm)
    addParameter(p, 'MaxFieldSpacingCm', 24, @isnumeric); % upper bound of candidate field-spacing search (~true spacing 24 cm)
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

    if all(cellfun(@isempty, {allData.TargetArea})), [allData.TargetArea] = deal(char(targetArea)); end
    if isfield(allData, 'Type') && ~isfield(allData, 'TypeImaged')
        [allData.TypeImaged] = allData.Type;
    elseif all(cellfun(@isempty, {allData.TypeImaged}))
        [allData.TypeImaged] = deal(char(cellType));
    end

    daysToPlot = p.Results.DaysToPlot;
    daysToPlot(daysToPlot == 200) = [];
    nDays = length(daysToPlot);

    % CHANGED: only baseline odd/even matrices are accumulated now — no
    % swap/omit condition handling at all.
    dayMatrixOdd = []; dayMatrixEven = [];

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
        if isempty(baseIdx), continue; end
        baseName = allCondsInDay{baseIdx};

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

        for s = 1:length(daySessions)
            thisSess = daySessions(s);

            if isfield(thisSess, 'FilteredROIs') && ~isempty(thisSess.FilteredROIs)
                idx = thisSess.FilteredROIs;
            else
                idx = 1:size(thisSess.ConditionData.(baseName).LapActivity, 1);
            end

            lapActivity = thisSess.ConditionData.(baseName).LapActivity(idx, :, :);
            if p.Results.ApplySmoothing, lapActivity = smoothLapActivity(lapActivity); end

            nTotalLaps   = size(lapActivity, 2);
            meanOddVals  = squeeze(mean(lapActivity(:, 1:2:nTotalLaps, :), 2, 'omitnan'));
            meanEvenVals = squeeze(mean(lapActivity(:, 2:2:nTotalLaps, :), 2, 'omitnan'));
            if size(meanOddVals,2) == 1, meanOddVals = meanOddVals'; meanEvenVals = meanEvenVals'; end

            dayMatrixOdd  = vertcat(dayMatrixOdd,  meanOddVals);
            dayMatrixEven = vertcat(dayMatrixEven, meanEvenVals);
        end
    end

    if isempty(dayMatrixOdd), figHandle = []; disp('No sessions matched.'); return; end

    %% Sorting and normalisation
    fullBaselineOdd  = dayMatrixOdd;
    minFullOdd       = min(fullBaselineOdd, [], 2, 'omitnan');
    maxFullOdd       = max(fullBaselineOdd, [], 2, 'omitnan');
    rangeFullOdd     = maxFullOdd - minFullOdd; rangeFullOdd(rangeFullOdd == 0) = 1;

    normFullOddRef   = (fullBaselineOdd - minFullOdd) ./ rangeFullOdd;
    normFullOddRef(isnan(normFullOddRef)) = 0;
    totalSomas       = size(fullBaselineOdd, 1);

    %% --- PERIODICITY-BASED CELL ISOLATION ---
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

    % --- PHASE-BASED SORT ---
    repPeriodBins = median(bestLagBins(isolatedGlobalIndices));

    [~, peakBinRaw] = max(smoothedForPeak, [], 2);
    wrappedPhaseBins = mod(peakBinRaw - 1, repPeriodBins); % -1/+1 handles 1-based MATLAB indexing
    [~, localSortIdx] = sort(wrappedPhaseBins, 'ascend');

    numN        = length(isolatedGlobalIndices);
    numBins     = size(isolateBaselineOdd, 2);
    cmPositions = linspace(0, 200, numBins);

    rawBaseEven  = dayMatrixEven(isolatedGlobalIndices, :);
    normBaseEven = (rawBaseEven - minIsolateOdd) ./ rangeIsolateOdd;
    normBaseEven(isnan(normBaseEven)) = 0;

    normBaseOdd  = (isolateBaselineOdd - minIsolateOdd) ./ rangeIsolateOdd;
    normBaseOdd(isnan(normBaseOdd)) = 0;

    %% Layout — simplified: two heatmaps (odd/even) plus a mean trace
    leftMargin  = 0.10;
    rightMargin = 0.12;
    colGap      = 0.03;
    cbW         = 0.02;
    cbGap       = 0.015;
    nCols       = 2;
    availW      = 1 - leftMargin - rightMargin - (nCols-1)*colGap;
    colW        = availW / nCols;
    heatmapBot  = 0.42; heatmapH = 0.50;
    traceBot    = 0.10; traceH = 0.22;
    colLefts    = leftMargin + (0:nCols-1) * (colW + colGap);
    cbLeft      = 1 - rightMargin + cbGap;

    figHandle = figure('Position', [50 50 900 800], 'Color', 'w');

    heatmapTitles = {'Baseline Odd', 'Baseline Even'};
    heatmapData   = {normBaseOdd, normBaseEven};

    for i = 1:nCols
        axAbs = axes('Position', [colLefts(i), heatmapBot, colW, heatmapH]);
        imagesc(axAbs, [0 200], [1 numN], heatmapData{i}(localSortIdx,:));
        colormap(axAbs, flipud(gray));
        set(axAbs, 'CLim', [0.25 0.75], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
        set(axAbs, 'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'}, 'FontName', targetFont, 'FontSize', 10);
        title(axAbs, heatmapTitles{i}, 'FontName', targetFont, 'FontSize', 12, 'FontWeight', 'bold');

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

    % Colorbar
    axAbsCB = axes('Position', [cbLeft, heatmapBot, cbW, heatmapH]);
    set(axAbsCB, 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
    colormap(axAbsCB, flipud(gray));
    cbAbs = colorbar(axAbsCB, 'Location', 'eastoutside');
    cbAbs.Position      = [cbLeft, heatmapBot, cbW, heatmapH];
    cbAbs.Ticks         = [0.25 0.50 0.75]; cbAbs.TickLabels = {'0.25','0.50','0.75'};
    cbAbs.TickDirection = 'out'; cbAbs.Box = 'off';
    cbAbs.FontName      = targetFont; cbAbs.FontSize = 9;
    cbAbs.Label.String  = 'Activity (norm.)';
    cbAbs.Label.FontName = targetFont; cbAbs.Label.FontSize = 11;

    % Mean trace (average across odd and even, single reference trace)
    axTrace = axes('Position', [colLefts(1), traceBot, availW, traceH]);
    hold(axTrace, 'on');

    combinedMean = mean(cat(3, normBaseOdd, normBaseEven), 3);
    muTrace  = mean(combinedMean, 1, 'omitnan');
    semTrace = std(combinedMean, 0, 1, 'omitnan') / sqrt(numN);

    fill(axTrace, [cmPositions, fliplr(cmPositions)], ...
        [muTrace+semTrace, fliplr(muTrace-semTrace)], ...
        [0 0 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(axTrace, cmPositions, muTrace, 'Color', [0 0 0], 'LineWidth', 1.5);

    for lm = [40 80 120 160]
        xline(axTrace, lm, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    end

    set(axTrace, 'Box', 'off', 'TickDir', 'out', 'XLim', [0 200], ...
        'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'}, ...
        'FontName', targetFont, 'FontSize', 10);
    xlabel(axTrace, 'Position (cm)', 'FontName', targetFont, 'FontSize', 11);
    ylabel(axTrace, 'Activity (norm.)', 'FontName', targetFont, 'FontSize', 11);

    if ~isempty(p.Results.SavePath)
        saveFigureFormats(figHandle, p.Results.SavePath);
    end
end

%% ---- Local helper functions (unchanged) ----

function acfMat = computeAutocorrBatch(mat)
    [nRows, nCols] = size(mat);
    nfft = 2^nextpow2(2*nCols);
    X    = fft(mat, nfft, 2);
    S    = X .* conj(X);
    acfFull = real(ifft(S, [], 2));
    rawAcf  = acfFull(:, 1:nCols);

    overlapCounts = nCols - (0:nCols-1);
    unbiasedAcf   = rawAcf ./ overlapCounts;

    lag0 = unbiasedAcf(:, 1);
    lag0(lag0 == 0) = eps;
    acfMat = unbiasedAcf ./ lag0;
end

function score = computeHarmonicScore(acfMat, bestLagBins, nHarmonics, minOverlapBins)
    [nRows, nCols] = size(acfMat);
    scoreSum    = zeros(nRows, 1);
    countValid  = zeros(nRows, 1);
    for k = 1:nHarmonics
        lagVal      = bestLagBins * k;
        colIdx      = lagVal + 1;
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
    [nRows, nCols] = size(mat);
    [~, idx] = sort(rand(nRows, nCols), 2);
    linIdx   = sub2ind([nRows, nCols], repmat((1:nRows)', 1, nCols), idx);
    out      = mat(linIdx);
end