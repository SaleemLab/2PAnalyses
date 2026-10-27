%% Speed Modulation Analysis — Full Corridor Version: in thesis!!
% Question: Does running speed modulate position-specific activity in RSP,
%           across the whole corridor (not just at landmarks)?
% Approach: Occupancy-based per-position median split (low vs high)
% Filters:  Spatially modulated ROIs only (crossValExpVar)
%           Edge peak cells excluded (~ExcludeEdgePeakCells)
% Sessions: All baseline trials only (from both baseline and manipulation corridors)
% Shuffle:  1000 lap-label shuffles per session to identify significant ROIs
%
% This is the agreed final version:
%   - Per-ROI speed modulation averaged over ALL position bins (full corridor,
%     not landmark-restricted).
%   - Lap-sufficiency check evaluated across the whole corridor.
%   - Raw (low - high) difference as the modulation metric — no ratio or
%     z-score normalization (the underlying DF/F is already per-ROI
%     z-scored across the full recording, so no further normalization is
%     needed; see prior discussion).
%   - Shuffle-test significance uses raw p<0.05, uncorrected for multiple
%     comparisons, matching the convention used elsewhere in the thesis
%     (crossValExpVar and GLM signicells are also uncorrected p-value
%     cutoffs). Stated explicitly wherever prevalence is reported.
%   - No pooled ttest/signrank (pseudoreplicated across sessions/animals).
%     Reporting is counts/percentages from the shuffle test, plus mean
%     low/high and Cohen's dz for effect size.
%   - Session speed/activity loading wrapped in try/catch — a session
%     missing lapPositionRunningSpeed, trialIndicesByCondition, or
%     lapPositionActivity.dFFNeuropilCorrected is skipped with a message
%     rather than crashing the whole batch run.
%   - Friedman test (bins as matched blocks) for whether median speed
%     thresholds differ significantly across sessions overall.
%   - Figures saved to outputDir, two subfolders: 'summary' (median-lines-
%     across-sessions plot and the pooled low/high histogram) and 'ref'
%     (per-session histograms, session speed summary bars, and the
%     example-session heatmap/tuning-curve figures).
%   - Fixed the pairs struct bug (again): square brackets concatenate
%     date strings into one char array instead of a cell array of
%     session names, and a single un-celled string can't be indexed with
%     {} at all. All entries now use curly braces.
%   - FIXED: allSessionMedians/allSessionLabels/allSessionIQRs were
%     appended for every session that loaded speed data, BEFORE the
%     lap-sufficiency check — including sessions later dropped. Session
%     speed summary plot/stats now filtered to only sessions actually
%     used in the pooled analysis (usedInPooling, based on
%     allSpeedModSession).
%   - FIXED: legend labels on the per-position median threshold plot were
%     pulled from allSessionLabels(validRows), which is indexed
%     differently from allMedianThreshLines — a silent misalignment bug.
%     Now correctly uses allMedianThreshLabels(validRows), which is built
%     in lockstep with allMedianThreshLines.
%   - ADDED: tracks the actual speed VALUES within each session's low- and
%     high-speed groups (not just the per-position median threshold), so
%     you can check whether "low speed" and "high speed" mean roughly the
%     same thing across sessions, not just whether the threshold differs.

%% Output folders
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section3_Fig3.6\runBeh_medianAcrossPos_withlapspeedpos';
summaryDir = fullfile(outputDir, 'summary');
refDir     = fullfile(outputDir, 'ref');
if ~exist(outputDir,  'dir'), mkdir(outputDir);  end
if ~exist(summaryDir, 'dir'), mkdir(summaryDir); end
if ~exist(refDir,     'dir'), mkdir(refDir);     end

%% Session pairs
pairs        = struct;
pairs.M25132 = {'20260226', '20260228', '20260313'};
pairs.M25133 = {'20260224'};
pairs.M26003 = {'20260322', '20260324', '20260325'};

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
allLowSpeedValuesBySession  = {}; % raw speed values within the low-speed group, kept per session
allHighSpeedValuesBySession = {}; % same for high-speed group

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

        %% Track actual speed values within low/high groups, not just the threshold
        lowSpeedValues  = [];
        highSpeedValues = [];
        for b = 1:numPosBins
            binSpeeds = speedData(:, b);
            validBin  = binSpeeds(binSpeeds > 1 & ~isnan(binSpeeds));
            if isempty(validBin), continue; end
            binMedian = median(validBin);
            lowSpeedValues  = [lowSpeedValues;  validBin(validBin <= binMedian)];
            highSpeedValues = [highSpeedValues; validBin(validBin > binMedian)];
        end
        allLowSpeedValuesBySession{end+1}  = lowSpeedValues;
        allHighSpeedValuesBySession{end+1} = highSpeedValues;

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

figMedianLines = figure('Color', 'w', 'Position', [100 100 900 400]); hold on;
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
legendLabels = [allMedianThreshLabels(validRows), {'Mean ± SD'}];
legend(legendLabels, 'Box', 'off', 'Location', 'best', 'FontSize', 9);
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial');
if exist('offsetAxes', 'file') == 2, offsetAxes(gca); end

if exist('saveFigureFormats', 'file') == 2
    saveFigureFormats(figMedianLines, fullfile(summaryDir, 'medianSpeedThreshold_acrossSessions'));
else
    saveas(figMedianLines, fullfile(summaryDir, 'medianSpeedThreshold_acrossSessions.png'));
end

%% Are median speed thresholds significantly different across sessions overall?
threshData = allMedianThreshLines(validRows, :);
validBins = all(~isnan(threshData), 1);
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

%% Compare mean speed within low/high groups across sessions
lowValuesUsed  = allLowSpeedValuesBySession(validRows);
highValuesUsed = allHighSpeedValuesBySession(validRows);

fprintf('\nMean speed within low/high groups (sessions used in pooled analysis only):\n');
lowMeans  = nan(length(lowValuesUsed), 1);
highMeans = nan(length(highValuesUsed), 1);
for s = 1:length(lowValuesUsed)
    lowMeans(s)  = mean(lowValuesUsed{s});
    highMeans(s) = mean(highValuesUsed{s});
    fprintf('  %s: low=%.1f cm/s, high=%.1f cm/s\n', allMedianThreshLabels{validRows(s)}, lowMeans(s), highMeans(s));
end
fprintf('Low-group mean range: %.1f-%.1f cm/s\n', min(lowMeans), max(lowMeans));
fprintf('High-group mean range: %.1f-%.1f cm/s\n', min(highMeans), max(highMeans));

lowPooled = []; lowGroupLabels = [];
for s = 1:length(lowValuesUsed)
    lowPooled = [lowPooled; lowValuesUsed{s}];
    lowGroupLabels = [lowGroupLabels; repmat(s, length(lowValuesUsed{s}), 1)];
end
pLow = kruskalwallis(lowPooled, lowGroupLabels, 'off');
fprintf('Kruskal-Wallis, low-speed values across sessions: p = %.4g\n', pLow);

highPooled = []; highGroupLabels = [];
for s = 1:length(highValuesUsed)
    highPooled = [highPooled; highValuesUsed{s}];
    highGroupLabels = [highGroupLabels; repmat(s, length(highValuesUsed{s}), 1)];
end
pHigh = kruskalwallis(highPooled, highGroupLabels, 'off');
fprintf('Kruskal-Wallis, high-speed values across sessions: p = %.4g\n', pHigh);

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
figPerSessionHist = figure('Position', [100 100 1200 600], 'Color', 'w');
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

if exist('saveFigureFormats', 'file') == 2
    saveFigureFormats(figPerSessionHist, fullfile(refDir, 'ref_perSessionHistograms'));
else
    saveas(figPerSessionHist, fullfile(refDir, 'ref_perSessionHistograms.png'));
end

%% Session speed summary plots
usedInPooling = ~cellfun(@isempty, allSpeedModSession);
fprintf('\n%d / %d sessions that loaded speed data were actually used in the pooled analysis\n', ...
        sum(usedInPooling), length(usedInPooling));

sessionMediansUsed = allSessionMedians(usedInPooling);
sessionIQRsUsed    = allSessionIQRs(usedInPooling);
sessionLabelsUsed  = allSessionLabels(usedInPooling);

fprintf('Session median running speed (sessions used in pooled analysis only): range %.1f-%.1f cm/s\n', ...
        min(sessionMediansUsed), max(sessionMediansUsed));

figSpeedSummary = figure('Position', [100 100 700 300], 'Color', 'w');
subplot(1,2,1)
bar(sessionMediansUsed, 'FaceColor', [0.4 0.6 0.8])
xticks(1:length(sessionLabelsUsed)); xticklabels(sessionLabelsUsed); xtickangle(45)
yline(mean(sessionMediansUsed), 'k--', 'HandleVisibility', 'off')
ylabel('Median speed (cm/s)'); title('Session median speeds (pooled sessions only)')
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(1,2,2)
bar(sessionIQRsUsed, 'FaceColor', [0.4 0.6 0.8])
xticks(1:length(sessionLabelsUsed)); xticklabels(sessionLabelsUsed); xtickangle(45)
ylabel('IQR (cm/s)'); title('Session speed IQR (pooled sessions only)')
set(gca, 'Box', 'off', 'TickDir', 'out')

if exist('saveFigureFormats', 'file') == 2
    saveFigureFormats(figSpeedSummary, fullfile(refDir, 'ref_sessionSpeedSummary'));
else
    saveas(figSpeedSummary, fullfile(refDir, 'ref_sessionSpeedSummary.png'));
end

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

cohens_dz = mean(allSpeedMod) / std(allSpeedMod);
fprintf('  Cohen''s dz   = %.4f\n', cohens_dz)

%% Pooled histogram with shuffle-significant ROIs highlighted
figPooledHist = figure('Position', [100 100 500 350], 'Color', 'w');
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

if exist('saveFigureFormats', 'file') == 2
    saveFigureFormats(figPooledHist, fullfile(summaryDir, 'pooledSpeedModulation_lowHigh_fullCorridor'));
else
    saveas(figPooledHist, fullfile(summaryDir, 'pooledSpeedModulation_lowHigh_fullCorridor.png'));
end

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

    figHeatmap = figure('Position', [100 100 800 700], 'Color', 'w');

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

    if exist('saveFigureFormats', 'file') == 2
        saveFigureFormats(figHeatmap, fullfile(refDir, sprintf('ref_exampleSessionHeatmap_%s', lastGood.label)));
    else
        saveas(figHeatmap, fullfile(refDir, sprintf('ref_exampleSessionHeatmap_%s.png', lastGood.label)));
    end

    %% Population mean tuning curves
    meanLow  = mean(filteredLow,  1, 'omitnan');
    meanHigh = mean(filteredHigh, 1, 'omitnan');
    meanDiff = meanLow - meanHigh;

    figTuningCurves = figure('Position', [100 100 800 400], 'Color', 'w');

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

    if exist('saveFigureFormats', 'file') == 2
        saveFigureFormats(figTuningCurves, fullfile(refDir, sprintf('ref_exampleSessionTuningCurves_%s', lastGood.label)));
    else
        saveas(figTuningCurves, fullfile(refDir, sprintf('ref_exampleSessionTuningCurves_%s.png', lastGood.label)));
    end
else
    fprintf('\nNo session data available for example-session plots.\n')
end

fprintf('\nMedian speed within low/high groups (sessions used in pooled analysis only):\n');
lowMedians  = nan(length(lowValuesUsed), 1);
highMedians = nan(length(highValuesUsed), 1);
for s = 1:length(lowValuesUsed)
    lowMedians(s)  = median(lowValuesUsed{s});
    highMedians(s) = median(highValuesUsed{s});
    fprintf('  %s: low=%.1f cm/s, high=%.1f cm/s\n', allMedianThreshLabels{validRows(s)}, lowMedians(s), highMedians(s));
end

lowOverallMedian  = median(lowMedians);
lowOverallMAD     = mad(lowMedians, 1);   % flag=1 -> median absolute deviation, matching your SMI convention
highOverallMedian = median(highMedians);
highOverallMAD    = mad(highMedians, 1);

fprintf('Low-group:  %.1f +/- %.1f cm/s (median +/- MAD across sessions)\n', lowOverallMedian, lowOverallMAD);
fprintf('High-group: %.1f +/- %.1f cm/s (median +/- MAD across sessions)\n', highOverallMedian, highOverallMAD);
