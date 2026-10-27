function plotExampleLandmarkTraces(RegionData, landmarkCentres, gratingLandmarks, plaidLandmarks, tolerance, windowRadius, regionLabel, targetSI, nPerCategory)
    % Selects nPerCategory example ROIs closest to each target Selectivity
    % Index (SI) value (default: 3 examples each at the 25th/50th/75th
    % percentiles of the population) and plots, for each one:
    %   Col 1 - the normalized odd/even trace used to compute SI
    %   Col 2 - smoothed per-lap dF/F, single laps overlaid
    %   Col 3 - trial-averaged smoothed dF/F (mean +/- SEM across laps)
    %
    % SI = 0 means the response at the preferred landmark equals the response
    % at the alternate-type landmark (two comparable bumps). SI = 1 means the
    % alternate-type landmark shows ~0 response (one dominant bump only).
    %
    % targetSI, if provided, is treated as literal SI values to match (not
    % percentiles) -- pass e.g. [0, 0.5, 1] to go back to fixed low/mid/high
    % targets instead of data-driven percentiles.
    %
    % Candidates are picked without replacement across categories (a cell
    % used for one category's examples won't also appear in another's).
    
    if nargin < 7 || isempty(regionLabel),   regionLabel   = 'Region'; end
    if nargin < 9 || isempty(nPerCategory),  nPerCategory  = 3;        end
    usePercentileDefault = (nargin < 8 || isempty(targetSI));
    
    %% --- Extraction pass: same logic as extractSplitMetrics, but keep session/ROI IDs and full traces ---
    w_space = gausswin(15); w_space = w_space / sum(w_space);
    
    gratingBins = []; for g = gratingLandmarks, gratingBins = [gratingBins, (g - windowRadius):(g + windowRadius)]; end
    plaidBins   = []; for p = plaidLandmarks,   plaidBins   = [plaidBins,   (p - windowRadius):(p + windowRadius)]; end
    
    allowedLandmarkBins = [];
    for c = landmarkCentres
        allowedLandmarkBins = [allowedLandmarkBins, (c - tolerance):(c + tolerance)];
    end
    allowedLandmarkBins = unique(allowedLandmarkBins);
    
    SIvals        = [];
    sessionIdxArr = [];
    roiIdxArr     = [];
    normTrainArr  = {};
    normTestArr   = {};
    primaryLmkArr = [];
    
    for s = 1:length(RegionData)
        sess = RegionData(s);
        if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
           ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs), continue; end
        
        lapActivity = sess.ConditionData.Baseline.LapActivity;
        [numROIsTotal, numLaps, numBins] = size(lapActivity);
        
        validSearchBins = allowedLandmarkBins(allowedLandmarkBins >= 1 & allowedLandmarkBins <= numBins);
        
        smoothedActivity = lapActivity;
        for iCell = 1:numROIsTotal
            for iLap = 1:numLaps
                trace = squeeze(lapActivity(iCell, iLap, :));
                if all(isnan(trace)), continue; end
                nanMask = isnan(trace); trace(nanMask) = 0;
                smoothed = filtfilt(w_space, 1, trace); smoothed(nanMask) = NaN;
                smoothedActivity(iCell, iLap, :) = smoothed;
            end
        end
        
        roisToAnalyze = sess.FilteredROIs;
        roiActivity   = smoothedActivity(roisToAnalyze, :, :);
        numROIs       = length(roisToAnalyze);
        if numROIs == 0, continue; end
        
        meanOdd  = squeeze(mean(roiActivity(:, 1:2:end, :), 2, 'omitnan'));
        meanEven = squeeze(mean(roiActivity(:, 2:2:end, :), 2, 'omitnan'));
        
        hasSavedExclusions = isfield(sess, 'SMI') && isfield(sess.SMI, 'ExcludeEdgePeakCells') && ...
                             length(sess.SMI.ExcludeEdgePeakCells) == numROIsTotal;
        
        for i = 1:numROIs
            if hasSavedExclusions
                if sess.SMI.ExcludeEdgePeakCells(roisToAnalyze(i)), continue; end
            end
            
            trainTrace = meanOdd(i, :);
            testTrace  = meanEven(i, :);
            if all(isnan(trainTrace)) || all(isnan(testTrace)), continue; end
            
            minOdd = min(trainTrace, [], 'omitnan');
            maxOdd = max(trainTrace, [], 'omitnan');
            rangeOdd = maxOdd - minOdd; if rangeOdd == 0, rangeOdd = 1; end
            
            normTrainOdd  = (trainTrace - minOdd) ./ rangeOdd;
            normTestEven  = (testTrace  - minOdd) ./ rangeOdd;
            
            [~, maxIdxInSearch] = max(normTrainOdd(validSearchBins));
            prefBin = validSearchBins(maxIdxInSearch);
            
            [~, closestLandmarkIdx] = min(abs(landmarkCentres - prefBin));
            primaryLandmark = landmarkCentres(closestLandmarkIdx);
            
            actG    = mean(max(0, normTestEven(gratingBins)), 'omitnan');
            actP    = mean(max(0, normTestEven(plaidBins)),   'omitnan');
            actPeak = max(0, normTestEven(prefBin));
            
            if ismember(primaryLandmark, gratingLandmarks)
                actAlternative = actP;
            else
                actAlternative = actG;
            end
            
            if (actPeak + actAlternative) > 0
                valSI = (actPeak - actAlternative) / (actPeak + actAlternative);
                if ~isnan(valSI)
                    valSI = max(0, min(1, valSI));
                    SIvals(end+1)        = valSI;                        %#ok<AGROW>
                    sessionIdxArr(end+1) = s;                            %#ok<AGROW>
                    roiIdxArr(end+1)     = roisToAnalyze(i);             %#ok<AGROW>
                    normTrainArr{end+1}  = normTrainOdd;                 %#ok<AGROW>
                    normTestArr{end+1}   = normTestEven;                 %#ok<AGROW>
                    primaryLmkArr(end+1) = primaryLandmark;              %#ok<AGROW>
                end
            end
        end
    end
    
    if isempty(SIvals)
        error('No cells passed filtering, cannot select examples.');
    end
    
    percentileTargets = [25, 50, 75];
    if usePercentileDefault
        targetSI    = prctile(SIvals, percentileTargets);
        categoryTag = {'25th pct', '50th pct (median)', '75th pct'};
    else
        categoryTag = arrayfun(@(x) sprintf('Target %.2f', x), targetSI, 'UniformOutput', false);
    end
    
    %% --- Select nPerCategory example ROIs closest to each target, without reusing a cell across categories ---
    usedMask   = false(1, length(SIvals));
    rowCatIdx  = [];  % which category each output row belongs to
    rowRank    = [];  % 1..nPerCategory within that category
    rowCellIdx = [];  % index into SIvals/sessionIdxArr/etc.
    
    for t = 1:length(targetSI)
        diffs = abs(SIvals - targetSI(t));
        diffs(usedMask) = Inf;
        [~, sortedIdx] = sort(diffs);
        nAvailable = sum(~isinf(diffs));
        nTake = min(nPerCategory, nAvailable);
        chosen = sortedIdx(1:nTake);
        usedMask(chosen) = true;
        
        rowCatIdx  = [rowCatIdx,  repmat(t, 1, nTake)];   %#ok<AGROW>
        rowRank    = [rowRank,    1:nTake];               %#ok<AGROW>
        rowCellIdx = [rowCellIdx, chosen];                %#ok<AGROW>
    end
    
    numRows = length(rowCellIdx);
    
    %% --- Print a candidate summary table for picking a final example ---
    fprintf('\n=== %s: Example candidates (%d per category) ===\n', regionLabel, nPerCategory);
    fprintf('%-22s %-6s %-10s %-8s %-10s\n', 'Category', 'Rank', 'Session', 'ROI', 'Actual SI');
    for r = 1:numRows
        k = rowCellIdx(r);
        fprintf('%-22s %-6d %-10d %-8d %-10.3f\n', ...
            categoryTag{rowCatIdx(r)}, rowRank(r), sessionIdxArr(k), roiIdxArr(k), SIvals(k));
    end
    
    %% --- Plot ---
    figHeight = min(1400, 230 * numRows / 3 + 150);  % scale height with row count, cap to keep it manageable
    fig = figure('Color', 'w', 'Position', [80 60 1100 figHeight]);
    tl = tiledlayout(fig, numRows, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s: Example ROIs Across Selectivity Index (%d per category)', regionLabel, nPerCategory), ...
        'FontName', 'Arial', 'FontSize', 13);

    for row = 1:numRows
        k = rowCellIdx(row);
        sessNum = sessionIdxArr(k);
        roiNum  = roiIdxArr(k);
        actualSI = SIvals(k);
        
        sess = RegionData(sessNum);
        rawTrace = squeeze(sess.ConditionData.Baseline.LapActivity(roiNum, :, :)); % [numLaps x numBins]
        numBins = size(rawTrace, 2);
        xvals = 1:numBins;
        
        % Per-lap Gaussian smoothing (same w_space window used in the extraction pass above)
        numLapsRaw = size(rawTrace, 1);
        for lap = 1:numLapsRaw
            trace = rawTrace(lap, :);
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace);
            trace(nanMask) = 0;
            smoothed = filtfilt(w_space, 1, trace);
            smoothed(nanMask) = NaN;
            rawTrace(lap, :) = smoothed;
        end
        
        rowTag = sprintf('%s #%d\nROI %d, session %d\nactual SI = %.3f', ...
            categoryTag{rowCatIdx(row)}, rowRank(row), roiNum, sessNum, actualSI);
        
        %% Col 1: normalized odd/even trace
        ax1 = nexttile(tl, (row-1)*3 + 1); hold(ax1, 'on');
        plot(ax1, xvals, normTrainArr{k}, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.3, 'DisplayName', 'train (odd laps)');
        plot(ax1, xvals, normTestArr{k},  '-',  'Color', 'k',           'LineWidth', 1.6, 'DisplayName', 'test (even laps)');
        ylim(ax1, [-0.1, 1.1]);
        addLandmarkShading(ax1, gratingLandmarks, plaidLandmarks, windowRadius);
        ylabel(ax1, rowTag, 'FontName', 'Arial', 'FontSize', 8);
        if row == 1, title(ax1, 'Normalized odd/even trace', 'FontWeight', 'normal', 'FontName', 'Arial'); end
        if row == 1, legend(ax1, 'Location', 'best', 'Box', 'off', 'FontSize', 8); end
        if row == numRows, xlabel(ax1, 'Track position (cm)', 'FontName', 'Arial'); end
        box(ax1, 'off');
        xticks([40, 80, 120, 160]);
        xticklabels({'40', '80', '120', '160'});
        
        %% Col 2: smoothed per-lap traces overlaid
        ax2 = nexttile(tl, (row-1)*3 + 2); hold(ax2, 'on');
        numLaps = size(rawTrace, 1);
        for lap = 1:numLaps
            plot(ax2, xvals, rawTrace(lap, :), 'Color', [0.5 0.5 0.5, 0.35], 'LineWidth', 0.7, 'HandleVisibility', 'off');
        end
        addLandmarkShading(ax2, gratingLandmarks, plaidLandmarks, windowRadius);
        if row == 1, title(ax2, sprintf('Smoothed per-lap dF/F (n=%d laps)', numLaps), 'FontWeight', 'normal', 'FontName', 'Arial'); end
        if row == numRows, xlabel(ax2, 'Track position (cm)', 'FontName', 'Arial'); end
        ylabel(ax2, 'dF/F', 'FontName', 'Arial', 'FontSize', 8);
        box(ax2, 'off');
        xticks([40, 80, 120, 160]);
        xticklabels({'40', '80', '120', '160'});
        
        %% Col 3: trial-averaged smoothed trace with shading
        ax3 = nexttile(tl, (row-1)*3 + 3); hold(ax3, 'on');
        meanTrace = mean(rawTrace, 1, 'omitnan');
        nLapsValid = sum(~isnan(rawTrace), 1);
        nLapsValid(nLapsValid == 0) = 1; % avoid divide-by-zero on all-NaN bins
        semTrace  = std(rawTrace, 0, 1, 'omitnan') ./ sqrt(nLapsValid);
        fill(ax3, [xvals, fliplr(xvals)], [meanTrace + semTrace, fliplr(meanTrace - semTrace)], ...
            [0.3 0.3 0.3], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax3, xvals, meanTrace, 'Color', 'k', 'LineWidth', 1.6);
        addLandmarkShading(ax3, gratingLandmarks, plaidLandmarks, windowRadius);
        if row == 1, title(ax3, 'Trial-averaged smoothed dF/F (mean +/- SEM)', 'FontWeight', 'normal', 'FontName', 'Arial'); end
        if row == numRows, xlabel(ax3, 'Track position (cm)', 'FontName', 'Arial'); end
        ylabel(ax3, 'dF/F', 'FontName', 'Arial', 'FontSize', 8);
        box(ax3, 'off');
        xticks([40, 80, 120, 160]);
        xticklabels({'40', '80', '120', '160'});
        
        for ax = [ax1, ax2, ax3]
            if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax); end
            if exist('offsetAxes', 'file') == 2, offsetAxes(ax); end
        end
    end
%     
    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\feature_selective_rsp_visp';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    saveFigureFormats(fig, fullfile(outputDir, sprintf('%s_example_landmark_traces_candidates', lower(regionLabel))));
end

%% Shade grating (blue) and plaid (orange) landmark windows behind the traces
function addLandmarkShading(ax, gratingLandmarks, plaidLandmarks, windowRadius)
    yl = ylim(ax);
    for g = gratingLandmarks
        ph = patch(ax, [g-windowRadius, g+windowRadius, g+windowRadius, g-windowRadius], ...
            [yl(1) yl(1) yl(2) yl(2)], [0.7 0.8 1], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        uistack(ph, 'bottom');
    end
    for p = plaidLandmarks
        ph = patch(ax, [p-windowRadius, p+windowRadius, p+windowRadius, p-windowRadius], ...
            [yl(1) yl(1) yl(2) yl(2)], [1 0.8 0.65], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        uistack(ph, 'bottom');
    end
    ylim(ax, yl);
end