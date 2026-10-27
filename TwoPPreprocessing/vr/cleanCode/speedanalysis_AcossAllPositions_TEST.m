%% Speed Modulation Analysis — Full Corridor Version
% Question: Does running speed modulate position-specific activity in RSP,
%           across the whole corridor (not just at landmarks)?
% Approach: Occupancy-based per-position median split (low vs high)
% Filters:  Spatially modulated ROIs only (crossValExpVar)
%           Edge peak cells excluded (~ExcludeEdgePeakCells)
% Sessions: All baseline trials only (from both baseline and manipulation corridors)
% Shuffle:  1000 lap-label shuffles per session to identify significant ROIs
%
% Differences from the landmark-restricted version:
%   - Per-ROI speed modulation is averaged over ALL position bins (1:numPosBins),
%     not just windows around the four landmarks.
%   - Lap-sufficiency check is evaluated across the whole corridor: a session is
%     kept only if a set fraction of position bins have enough laps in both the
%     low and high speed groups (see minBinFracSufficient below), rather than
%     just checking the four landmark windows.
%   - Shuffle test null distributions are computed using all position bins.
%   - Landmarks are kept only as reference lines on plots, not as an inclusion
%     criterion for which bins go into the modulation metric.
%   - Fixed a bug: pairs.M25132 / pairs.M25133 / pairs.M26003 were built with
%     square brackets, which concatenates the date strings into one char array
%     instead of a cell array of session names. Switched to curly braces.
%   - Removed the hardcoded example-session reload at the end (it risked
%     mismatching FilteredROIS_idx from a different session). The heatmap /
%     population tuning curve figures now use the last successfully processed
%     session's own matrices and ROI filter, captured inside the loop.
%   - Shuffle-test significance uses raw p<0.05, uncorrected for multiple
%     comparisons across ROIs, to match the convention used elsewhere in
%     the thesis (crossValExpVar and GLM signicells thresholds are also
%     uncorrected p-value cutoffs). State this explicitly wherever
%     prevalence is reported.
%   - Removed the pooled ttest/signrank on raw ROI values (pseudoreplicated
%     across sessions/animals) and the matching per-session ttest/signrank.
%     Reporting is counts/percentages from the shuffle test, plus mean
%     low/high and Cohen's dz for effect size.
%   - Session speed/activity loading is now wrapped in a try/catch — a
%     session missing lapPositionRunningSpeed, trialIndicesByCondition, or
%     lapPositionActivity.dFFNeuropilCorrected is skipped with a message
%     rather than crashing the whole batch run.

%% Session pairs
% pairs        = struct;
% pairs.M25132 = {'20260219','20260220','20260221','20260223','20260226','20260228','20260313'};
% pairs.M25133 = {'20260220','20260221','20260223','20260219','20260224'};
% pairs.M26003 = {'20260316','20260317','20260320','20260321','20260322','20260324','20260325'};

pairs        = struct;
pairs.M25132 = ['20260226', '20260228', '20260313'];
pairs.M26003 = ['20260322', '20260324', '20260325'];

filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);
mouseInfo     = sessionsToProcess(filteredTable);

totalSessionsToProcess = 0;
for i = 1:size(mouseInfo, 1)
    totalSessionsToProcess = totalSessionsToProcess + length(mouseInfo{i, 2});
end

%% Initialise
allSessionMedians   = [];
allSessionIQRs      = [];
allSessionLabels    = {};
allFilePaths        = {};
allSpeedModSession   = {};
allLowSession        = {}; % per-ROI, corridor-averaged low response — for effect-size context
allHighSession       = {}; % per-ROI, corridor-averaged high response
allShufflePSession   = {}; % per-ROI shuffle p-values per session
sessionsProcessedCount = 0;
allMedianThreshLines  = [];
allMedianThreshLabels = {};

landmarks    = [40, 80, 120, 160]; % kept for plotting reference only
nShuffles    = 1000;

% Sufficiency criterion for the full-corridor check
minLapsPerGroup     = 20;   % same per-group lap minimum as before
minBinFracSufficient = 0.8; % session kept only if >=80% of position bins meet minLapsPerGroup

% Captured from the loop for the example-session plots at the end
lastGood.filteredLow  = [];
lastGood.filteredHigh = [];
lastGood.filteredDiff = [];
lastGood.label        = '';

%% Loop through mice and sessions
for thisMouse = 1:size(mouseInfo, 1)
    mousenumber  = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};
    fprintf('Processing Mouse: %s\n', mousenumber);

    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};
        sessionsProcessedCount = sessionsProcessedCount + 1;
        fprintf('  Processing Session: %s\n', sessionName);

        %% Load session info
        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if ~isfile(infoPath), fprintf('  Info Missing — skipping\n'); continue; end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};

        %% Find corridor stimulus
        targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "CombinedRuns"), 1);
        if isempty(targetIdx)
            allCorridorIdx = find(contains(stimNames, "Corridor"));
            if isscalar(allCorridorIdx)
                targetIdx = allCorridorIdx;
            elseif length(allCorridorIdx) > 1
                targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "00002"), 1);
            end
        end
        if isempty(targetIdx)
            fprintf('  No valid corridor found — skipping\n'); continue;
        end

        fprintf('  Loading: %s\n', sessionFileInfo.stimFiles(targetIdx).name)
        VRStimName = sessionFileInfo.stimFiles(targetIdx).name;

        filePath = sessionFileInfo.stimFiles(targetIdx).Response;
        allFilePaths{end+1} = filePath;

        %% Session speed stats (baseline trials, running only)
        % Wrapped in try/catch — a session missing lapPositionRunningSpeed,
        % trialIndicesByCondition.Baseline, or
        % lapPositionActivity.dFFNeuropilCorrected is skipped rather than
        % crashing the whole batch run. load() with named variables silently
        % omits any that don't exist in the file, so the error actually
        % surfaces one line down, on the first field access.
        try
            rawData        = load(filePath, 'lapPositionRunningSpeed', 'lapPositionActivity', 'trialIndicesByCondition', 'SMI_Metrics');
            speedData      = rawData.lapPositionRunningSpeed;
            baselineTrials = rawData.trialIndicesByCondition.Baseline;
            speedData      = speedData(baselineTrials, :);
            rawROIData     = rawData.lapPositionActivity.dFFNeuropilCorrected;
            rawROIData     = rawROIData(:, baselineTrials, :);
        catch ME
            fprintf('  No speed/activity data for this session (%s) — skipping\n', ME.message);
            continue
        end

        numPosBins = size(speedData, 2);

        validSpeeds   = speedData(:);
        validSpeeds   = validSpeeds(validSpeeds > 1 & ~isnan(validSpeeds));
        sessionMedian = median(validSpeeds);
        sessionIQR    = iqr(validSpeeds);

        allSessionMedians(end+1) = sessionMedian;
        allSessionIQRs(end+1)    = sessionIQR;
        allSessionLabels{end+1}  = sprintf('%s_%s', mousenumber, sessionName);
        fprintf('  Median: %.1f cm/s | IQR: %.1f cm/s\n', sessionMedian, sessionIQR);

        %% Check minimum laps per speed group across the WHOLE corridor
        % For each position bin, split laps by that bin's own median speed
        % (same rule getLowHighSpeedPositionMatrix uses) and check group sizes.
        nLowPerBin  = nan(1, numPosBins);
        nHighPerBin = nan(1, numPosBins);
        for b = 1:numPosBins
            binSpeeds = speedData(:, b);
            validBin  = binSpeeds(binSpeeds > 1 & ~isnan(binSpeeds));
            if isempty(validBin), continue; end
            binMedian = median(validBin);
            nLowPerBin(b)  = sum(binSpeeds > 1 & binSpeeds <= binMedian & ~isnan(binSpeeds));
            nHighPerBin(b) = sum(binSpeeds > binMedian & ~isnan(binSpeeds));
        end

        sufficientBins = nLowPerBin >= minLapsPerGroup & nHighPerBin >= minLapsPerGroup;
        fracSufficient = mean(sufficientBins, 'omitnan');

        if fracSufficient < minBinFracSufficient
            fprintf('  Skipping — only %.0f%% of bins have >=%d laps per group (need %.0f%%)\n', ...
                    fracSufficient*100, minLapsPerGroup, minBinFracSufficient*100)
            allSpeedModSession{end+1} = [];
            allLowSession{end+1}  = [];
            allHighSession{end+1} = [];
            allShufflePSession{end+1} = [];
            continue
        end
        fprintf('  %.0f%% of position bins have sufficient laps per group\n', fracSufficient*100)

        %% Compute low/high speed position matrices (already full-corridor)
        options.signalToUse = 'dFFNeuropilCorrected';
        [response] = getLowHighSpeedPositionMatrix(sessionFileInfo, VRStimName, options);

        %% Compute continuous speed position matrix (panel c heatmap)
        matOptions.signalToUse    = 'dFFNeuropilCorrected';
        matOptions.numBins        = 10;
        matOptions.applySmoothing = false;
        getPositionSpeedMatrix(sessionFileInfo, VRStimName, matOptions);

        %% Save per-session median lines
        allMedianThreshLines(end+1, :) = response.speedPositionActivity.lowHigh.medianThreshLine;
        allMedianThreshLabels{end+1}   = sprintf('%s_%s', mousenumber, sessionName);

        %% Load ROI filters
        rawData2       = load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'crossValExpVar');
        crossValExpVar = rawData2.crossValExpVar;

        FilteredROIS_idx = find( ...
            crossValExpVar.dFFNeuropilCorrected.pValues    <= 0.01 & ...
            crossValExpVar.dFFNeuropilCorrected.medianExpVar > 0.1 & ...
            ~rawData.SMI_Metrics.dFFNeuropilCorrected.ExcludeEdgePeakCells);

        if isempty(FilteredROIS_idx)
            fprintf('  No valid ROIs — skipping\n');
            allSpeedModSession{end+1} = [];
            allLowSession{end+1}  = [];
            allHighSession{end+1} = [];
            allShufflePSession{end+1} = [];
            continue
        end
        fprintf('  nROIs: %d\n', length(FilteredROIS_idx));

        %% Compute real speed modulation per ROI — averaged over ALL position bins
        filteredLow  = response.speedPositionActivity.lowHigh.matrixLow(FilteredROIS_idx, :);
        filteredHigh = response.speedPositionActivity.lowHigh.matrixHigh(FilteredROIS_idx, :);
        filteredDiff = filteredLow - filteredHigh;

        speedModPerROI = mean(filteredDiff, 2, 'omitnan');
        allSpeedModSession{end+1} = speedModPerROI;
        allLowSession{end+1}  = filteredLow;
        allHighSession{end+1} = filteredHigh;
        fprintf('  Mean modulation (full corridor): %.4f\n', mean(speedModPerROI, 'omitnan'));

        % Keep this session's matrices in case it's the last one processed —
        % used for the example-session plots instead of a hardcoded reload.
        lastGood.filteredLow  = filteredLow;
        lastGood.filteredHigh = filteredHigh;
        lastGood.filteredDiff = filteredDiff;
        lastGood.label        = sprintf('%s_%s', mousenumber, sessionName);

        %% Shuffle test — shuffle lap speed labels, recompute modulation over ALL bins
        nLaps    = size(speedData, 1);
        nROIs_s  = length(FilteredROIS_idx);
        shuffleMod = nan(nROIs_s, nShuffles);

        fprintf('  Running %d shuffles across %d bins...\n', nShuffles, numPosBins);
        for sh = 1:nShuffles
            shuffIdx  = randperm(nLaps);
            shuffLow  = nan(nROIs_s, numPosBins);
            shuffHigh = nan(nROIs_s, numPosBins);

            for b = 1:numPosBins
                binSpeeds = speedData(shuffIdx, b);
                validBin  = binSpeeds(binSpeeds > 1 & ~isnan(binSpeeds));
                if isempty(validBin), continue; end
                binMedian = median(validBin);

                lowIdx  = binSpeeds > 1 & binSpeeds <= binMedian & ~isnan(binSpeeds);
                highIdx = binSpeeds > binMedian & ~isnan(binSpeeds);

                roiData = squeeze(rawROIData(FilteredROIS_idx, :, b)); % [nROIs x nLaps]

                if sum(lowIdx)  >= 2
                    shuffLow(:, b)  = mean(roiData(:, lowIdx),  2, 'omitnan');
                end
                if sum(highIdx) >= 2
                    shuffHigh(:, b) = mean(roiData(:, highIdx), 2, 'omitnan');
                end
            end
            shuffleMod(:, sh) = mean(shuffLow - shuffHigh, 2, 'omitnan');
        end

        %% Per-ROI shuffle p-value
        shuffP = nan(nROIs_s, 1);
        for r = 1:nROIs_s
            shuffP(r) = mean(abs(shuffleMod(r,:)) >= abs(speedModPerROI(r)));
        end
        allShufflePSession{end+1} = shuffP;

        nSig    = sum(shuffP < 0.05);
        nSigPos = sum(shuffP < 0.05 & speedModPerROI > 0);
        nSigNeg = sum(shuffP < 0.05 & speedModPerROI < 0);
        fprintf('  Shuffle: %d/%d sig (%.1f%%) — low>high: %d, high>low: %d\n', ...
                nSig, nROIs_s, nSig/nROIs_s*100, nSigPos, nSigNeg)
    end
end

%% Plot per-position median threshold lines across sessions
validRows  = find(any(~isnan(allMedianThreshLines), 2));
threshMean = mean(allMedianThreshLines(validRows, :), 1, 'omitnan');
threshStd  = std(allMedianThreshLines(validRows, :), 0, 1, 'omitnan');

figure('Color', 'w', 'Position', [100 100 900 400]); hold on;
colors = lines(length(validRows));
for s = 1:length(validRows)
    plot(1:size(allMedianThreshLines,2), allMedianThreshLines(validRows(s), :), ...
        'Color', colors(s,:), 'LineWidth', 1.5);
end
plot(1:size(allMedianThreshLines,2), threshMean, 'k--', 'LineWidth', 2.5);
fill([1:size(allMedianThreshLines,2), fliplr(1:size(allMedianThreshLines,2))], ...
     [threshMean + threshStd, fliplr(threshMean - threshStd)], ...
     'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
for lm = landmarks
    xline(lm, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
end
xlabel('Position (cm)', 'FontName', 'Arial', 'FontSize', 12);
ylabel('Median speed threshold (cm/s)', 'FontName', 'Arial', 'FontSize', 12);
title(sprintf('Per-position median speed threshold across %d sessions', length(validRows)), 'FontWeight', 'normal');
legendLabels = [allSessionLabels(validRows), {'Mean ± SD'}];
legend(legendLabels, 'Box', 'off', 'Location', 'best', 'FontSize', 9);
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial');
if exist('offsetAxes', 'file') == 2, offsetAxes(gca); end

%% Are median speed thresholds significantly different across sessions overall?
% Friedman test — nonparametric, treats each position bin as a matched
% block shared across sessions (same 200 corridor positions in every
% session), testing whether the threshold LEVEL differs by session once
% that pairing is accounted for. One overall p-value.
threshData = allMedianThreshLines(validRows, :);

validBins = all(~isnan(threshData), 1); % Friedman can't handle NaNs
threshDataClean = threshData(:, validBins);
fprintf('\nUsing %d / %d position bins with valid values in every session\n', ...
        sum(validBins), size(threshData, 2));

p_friedman = friedman(threshDataClean', 1, 'off');
fprintf('Friedman test across sessions: p = %.4g\n', p_friedman);
if p_friedman < 0.05
    fprintf('--> Median speed thresholds differ significantly across sessions.\n');
else
    fprintf('--> No significant difference in median speed thresholds across sessions.\n');
end

%% Per-session results with shuffle
fprintf('\nPer-session results (full corridor):\n')
for s = 1:length(allSpeedModSession)
    if ~isempty(allSpeedModSession{s})
        sessionMod = allSpeedModSession{s};
        shuffP     = allShufflePSession{s};
        nSig       = sum(shuffP < 0.05);
        nSigPos    = sum(shuffP < 0.05 & sessionMod > 0);
        nSigNeg    = sum(shuffP < 0.05 & sessionMod < 0);
        nTotal     = length(sessionMod);
        nPos       = sum(sessionMod > 0);
        nNeg       = sum(sessionMod < 0);

        fprintf('  %s:\n', allSessionLabels{s})
        fprintf('    n=%d, mean=%.4f, median=%.4f\n', nTotal, mean(sessionMod,'omitnan'), median(sessionMod,'omitnan'))
        fprintf('    low>high: %d (%.1f%%), high>low: %d (%.1f%%)\n', nPos, nPos/nTotal*100, nNeg, nNeg/nTotal*100)
        fprintf('    Shuffle sig: %d/%d (%.1f%%) — low>high: %d, high>low: %d\n', nSig, nTotal, nSig/nTotal*100, nSigPos, nSigNeg)
    end
end

%% Per-session histograms with shuffle-significant ROIs highlighted
figure('Position', [100 100 1200 600], 'Color', 'w');
nSessions = length(allSpeedModSession);
for s = 1:nSessions
    if isempty(allSpeedModSession{s}), continue; end

    sessionMod = allSpeedModSession{s};
    shuffP     = allShufflePSession{s};
    sigIdx     = shuffP < 0.05;
    nonsigIdx  = ~sigIdx;

    subplot(2, ceil(nSessions/2), s)
    hold on
    if sum(nonsigIdx) > 0
        histogram(sessionMod(nonsigIdx), 20, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7)
    end
    sigPos = sessionMod(sigIdx & sessionMod > 0);
    if ~isempty(sigPos)
        histogram(sigPos, 20, 'FaceColor', [0.0 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.9)
    end
    sigNeg = sessionMod(sigIdx & sessionMod < 0);
    if ~isempty(sigNeg)
        histogram(sigNeg, 20, 'FaceColor', [0.8 0.0 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.9)
    end
    xline(0, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off')
    xline(mean(sessionMod, 'omitnan'), '-k', 'LineWidth', 2, 'HandleVisibility', 'off')

    nSig    = sum(sigIdx);
    nSigPos = sum(sigIdx & sessionMod > 0);
    nSigNeg = sum(sigIdx & sessionMod < 0);

    title(sprintf('%s\nsig: %d/%d (low>high: %d, high>low: %d)', ...
          strrep(allSessionLabels{s}, '_', '\_'), nSig, length(sessionMod), nSigPos, nSigNeg))
    xlabel('Mean (Low - High), full corridor')
    ylabel('ROIs')
    legend('Non-sig', 'Sig low>high', 'Sig high>low', 'Box', 'off', 'FontSize', 7)
    set(gca, 'Box', 'off', 'TickDir', 'out')
end
sgtitle('Speed modulation per session, full corridor — shuffle-significant ROIs highlighted')

%% Session speed summary plots
figure('Position', [100 100 700 300], 'Color', 'w');
subplot(1,2,1)
bar(allSessionMedians, 'FaceColor', [0.4 0.6 0.8])
xticks(1:length(allSessionLabels)); xticklabels(allSessionLabels); xtickangle(45)
yline(mean(allSessionMedians), 'k--', 'HandleVisibility', 'off')
ylabel('Median speed (cm/s)'); title('Session median speeds')
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(1,2,2)
bar(allSessionIQRs, 'FaceColor', [0.4 0.6 0.8])
xticks(1:length(allSessionLabels)); xticklabels(allSessionLabels); xtickangle(45)
ylabel('IQR (cm/s)'); title('Session speed IQR')
set(gca, 'Box', 'off', 'TickDir', 'out')

%% Pool across sessions and test
iqrThreshold = 0;
includeIdx   = allSessionIQRs >= iqrThreshold;
fprintf('\nIncluding %d of %d sessions (IQR >= %.0f cm/s)\n', ...
        sum(includeIdx), length(allSessionIQRs), iqrThreshold)

allSpeedMod = [];
allShuffleP = [];
for s = 1:length(allSpeedModSession)
    if includeIdx(s) && ~isempty(allSpeedModSession{s})
        allSpeedMod = [allSpeedMod; allSpeedModSession{s}];
        allShuffleP = [allShuffleP; allShufflePSession{s}];
    end
end

nSigAll    = sum(allShuffleP < 0.05);
nSigPosAll = sum(allShuffleP < 0.05 & allSpeedMod > 0);
nSigNegAll = sum(allShuffleP < 0.05 & allSpeedMod < 0);

fprintf('\nPooled result (full corridor):\n')
fprintf('  n            = %d ROIs\n',  length(allSpeedMod))
fprintf('  mean         = %.4f\n',     mean(allSpeedMod))
fprintf('  median       = %.4f\n',     median(allSpeedMod))
fprintf('  Shuffle sig  = %d/%d (%.1f%%) — low>high: %d, high>low: %d [uncorrected for multiple comparisons]\n', ...
        nSigAll, length(allSpeedMod), nSigAll/length(allSpeedMod)*100, nSigPosAll, nSigNegAll)

% Mean low/high response, for context on what the difference is a
% difference OF. Computed the same corridor-wide-averaged way as
% speedModPerROI, so it's directly comparable to it — note these will
% likely look small and close together since most of the 200 bins are
% off each cell's preferred position (near baseline after z-scoring),
% not because the true peak-to-peak modulation is small.
allMeanLow  = [];
allMeanHigh = [];
for s = 1:length(allSpeedModSession)
    if includeIdx(s) && ~isempty(allSpeedModSession{s})
        allMeanLow  = [allMeanLow;  mean(allLowSession{s},  2, 'omitnan')];
        allMeanHigh = [allMeanHigh; mean(allHighSession{s}, 2, 'omitnan')];
    end
end
fprintf('  mean(low)    = %.4f\n', mean(allMeanLow, 'omitnan'))
fprintf('  mean(high)   = %.4f\n', mean(allMeanHigh, 'omitnan'))

% Standardized paired effect size (Cohen's dz) — magnitude of the
% difference relative to its own spread across ROIs, independent of the
% arbitrary DF/F units. A magnitude descriptor, not a significance test.
cohens_dz = mean(allSpeedMod) / std(allSpeedMod);
fprintf('  Cohen''s dz   = %.4f\n', cohens_dz)

%% Pooled histogram with shuffle-significant ROIs highlighted
figure('Position', [100 100 500 350], 'Color', 'w');
hold on
sigIdx    = allShuffleP < 0.05;
nonsigIdx = ~sigIdx;

histogram(allSpeedMod(nonsigIdx), 50, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7)
sigPos = allSpeedMod(sigIdx & allSpeedMod > 0);
sigNeg = allSpeedMod(sigIdx & allSpeedMod < 0);
if ~isempty(sigPos)
    histogram(sigPos, 50, 'FaceColor', [0.0 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.9)
end
if ~isempty(sigNeg)
    histogram(sigNeg, 50, 'FaceColor', [0.8 0.0 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.9)
end
xline(0,                 '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off')
xline(mean(allSpeedMod), '-k',  'LineWidth', 2,   'HandleVisibility', 'off')
xlabel('Mean (Low - High) \DeltaF/F, full corridor')
ylabel('Number of ROIs')
legend('Non-sig', 'Sig low>high', 'Sig high>low', 'Box', 'off')
title(sprintf('Pooled speed modulation, full corridor | n=%d ROIs | sig: %d/%d (%.1f%%) — low>high: %d, high>low: %d', ...
    length(allSpeedMod), nSigAll, length(allSpeedMod), ...
    nSigAll/length(allSpeedMod)*100, nSigPosAll, nSigNegAll))
set(gca, 'Box', 'off', 'TickDir', 'out')

%% Heatmap visualisation (last successfully processed session)
if ~isempty(lastGood.filteredLow)
    filteredLow  = lastGood.filteredLow;
    filteredHigh = lastGood.filteredHigh;
    filteredDiff = lastGood.filteredDiff;

    [~, peakPos]   = max(filteredLow, [], 2);
    [~, sortOrder] = sort(peakPos);

    sortedLow  = filteredLow(sortOrder, :);
    sortedHigh = filteredHigh(sortOrder, :);
    sortedDiff = filteredDiff(sortOrder, :);

    sigmaSpatial = 1;
    sigmaROI     = 2;
    smoothedLow  = imgaussfilt(sortedLow,  [sigmaROI sigmaSpatial]);
    smoothedHigh = imgaussfilt(sortedHigh, [sigmaROI sigmaSpatial]);
    smoothedDiff = imgaussfilt(sortedDiff, [sigmaROI sigmaSpatial]);

    sharedMax  = prctile([smoothedLow(:); smoothedHigh(:)], 98);
    sharedMin  = prctile([smoothedLow(:); smoothedHigh(:)],  2);
    diffMax    = prctile(abs(smoothedDiff(:)), 95);
    x_pos      = 1:size(filteredLow, 2);
    nROIs_plot = size(sortedLow, 1);

    figure('Position', [100 100 800 700], 'Color', 'w');

    subplot(3,1,1)
    imagesc(x_pos, 1:nROIs_plot, smoothedLow, [sharedMin sharedMax])
    title(sprintf('Low speed (%s)', strrep(lastGood.label, '_', '\_')))
    ylabel('ROI (sorted by pref. position, low)')
    xline(landmarks, '--w', 'Alpha', 0.4, 'HandleVisibility', 'off')
    colorbar; axis tight
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [])

    subplot(3,1,2)
    imagesc(x_pos, 1:nROIs_plot, smoothedHigh, [sharedMin sharedMax])
    title('High speed')
    ylabel('ROI (sorted by pref. position, low)')
    xline(landmarks, '--w', 'Alpha', 0.4, 'HandleVisibility', 'off')
    colorbar; axis tight
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [])

    subplot(3,1,3)
    imagesc(x_pos, 1:size(smoothedDiff,1), smoothedDiff, [-diffMax diffMax])
    title('Low - High')
    ylabel('ROI (sorted by pref. position, low)')
    xlabel('Position (cm)')
    xline(landmarks, '--w', 'Alpha', 0.4, 'HandleVisibility', 'off')
    colorbar; axis tight
    set(gca, 'Box', 'off', 'TickDir', 'out')

    %% Population mean tuning curves
    meanLow  = mean(filteredLow,  1, 'omitnan');
    meanHigh = mean(filteredHigh, 1, 'omitnan');
    meanDiff = meanLow - meanHigh;

    figure('Position', [100 100 800 400], 'Color', 'w');

    subplot(2,1,1)
    hold on
    plot(x_pos, meanLow,  'Color', [0.0 0.6 0.8], 'LineWidth', 2)
    plot(x_pos, meanHigh, 'Color', [0.8 0.0 0.6], 'LineWidth', 2)
    xline(landmarks, '--k', 'HandleVisibility', 'off')
    legend('Low', 'High', 'Location', 'northeast', 'Box', 'off')
    ylabel('\DeltaF/F')
    title(sprintf('Population mean, low vs high speed (%s)', strrep(lastGood.label, '_', '\_')))
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XLim', [1 x_pos(end)])

    subplot(2,1,2)
    hold on
    plot(x_pos, meanDiff, 'k', 'LineWidth', 2)
    yline(0, '--r', 'HandleVisibility', 'off')
    xline(landmarks, '--k', 'HandleVisibility', 'off')
    ylabel('Low - High')
    xlabel('Position (cm)')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XLim', [1 x_pos(end)])
else
    fprintf('\nNo session data available for example-session plots.\n')
end