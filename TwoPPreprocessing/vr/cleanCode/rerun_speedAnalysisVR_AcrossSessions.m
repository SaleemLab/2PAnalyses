%% Speed Modulation Analysis
% Question: Does running speed modulate position-specific activity in RSP?
% Approach: Occupancy-based per-position median split (low vs high)
% Filters:  Spatially modulated ROIs only (crossValExpVar)
%           Edge peak cells excluded (~ExcludeEdgePeakCells)
% Sessions: All baseline trials only (from both baseline and manipulation corridors)
% Shuffle:  1000 lap-label shuffles per session to identify significant ROIs

%% Session pairs
pairs        = struct;
pairs.M25132 = ['20260226', '20260228', '20260313'];
pairs.M26003 = ['20260322', '20260324', '20260325'];

filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);
mouseInfo     = sessionsToProcess(filteredTable);

totalSessionsToProcess = 0;
sessionsProcessedCount = 0;
for i = 1:size(mouseInfo, 1)
    totalSessionsToProcess = totalSessionsToProcess + length(mouseInfo{i, 2});
end

%% Initialise
allSessionMedians  = [];
allSessionIQRs     = [];
allSessionLabels   = {};
allFilePaths       = {};
allSpeedModSession = {};
allShufflePSession = {}; % per-ROI shuffle p-values per session
sessionsProcessedCount = 0;
allMedianThreshLines = [];
allMedianThreshLabels = {};

landmarkBins = [30:55, 65:95, 105:135, 145:175];
landmarks    = [40, 80, 120, 160];
nShuffles    = 1000;

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

        % run get speed position matrix 

        filePath   = sessionFileInfo.stimFiles(targetIdx).Response;
        allFilePaths{end+1} = filePath;

        %% Session speed stats (baseline trials, running only)
        rawData        = load(filePath, 'lapPositionRunningSpeed', 'lapPositionActivity', 'trialIndicesByCondition', 'SMI_Metrics');
        speedData      = rawData.lapPositionRunningSpeed;
        baselineTrials = rawData.trialIndicesByCondition.Baseline;
        speedData      = speedData(baselineTrials, :);
        rawROIData     = rawData.lapPositionActivity.dFFNeuropilCorrected;
        rawROIData     = rawROIData(:, baselineTrials, :);

        validSpeeds   = speedData(:);
        validSpeeds   = validSpeeds(validSpeeds > 1 & ~isnan(validSpeeds));
        sessionMedian = median(validSpeeds);
        sessionIQR    = iqr(validSpeeds);

        allSessionMedians(end+1) = sessionMedian;
        allSessionIQRs(end+1)    = sessionIQR;
        allSessionLabels{end+1}  = sprintf('%s_%s', mousenumber, sessionName);
        fprintf('  Median: %.1f cm/s | IQR: %.1f cm/s\n', sessionMedian, sessionIQR);


        %% Check minimum laps per speed group at landmark positions
        minLapsPerGroup = 20;
        sufficientLaps  = true;

        for lm = 1:length(landmarks)
            window   = landmarks(lm)-15 : landmarks(lm)+15;
            lmSpeeds = mean(speedData(:, window), 2);
            lmMedian = median(lmSpeeds(lmSpeeds > 1 & ~isnan(lmSpeeds)));
            nLow     = sum(lmSpeeds > 1 & lmSpeeds <= lmMedian);
            nHigh    = sum(lmSpeeds > lmMedian);

            if nLow < minLapsPerGroup || nHigh < minLapsPerGroup
                sufficientLaps = false;
                fprintf('  Skipping — insufficient laps at LM %d cm: nLow=%d, nHigh=%d\n', landmarks(lm), nLow, nHigh)
                break
            end
        end

        if ~sufficientLaps
            allSpeedModSession{end+1} = [];
            allShufflePSession{end+1} = [];
            continue
        end

        %% Compute low/high speed position matrices
        options.signalToUse = 'dFFNeuropilCorrected';
        [response] = getLowHighSpeedPositionMatrix(sessionFileInfo, VRStimName, options);

        %% 
        %% Compute continuous speed position matrix (panel c heatmap)
        matOptions.signalToUse    = 'dFFNeuropilCorrected';
        matOptions.numBins        = 10;
        matOptions.applySmoothing = false;
        getPositionSpeedMatrix(sessionFileInfo, VRStimName, matOptions);


        %% save per sessions median lines

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
            allShufflePSession{end+1} = [];
            continue
        end
        fprintf('  nROIs: %d\n', length(FilteredROIS_idx));

        %% Compute real speed modulation per ROI
        filteredLow  = response.speedPositionActivity.lowHigh.matrixLow(FilteredROIS_idx, :);
        filteredHigh = response.speedPositionActivity.lowHigh.matrixHigh(FilteredROIS_idx, :);

        filteredLow_lm  = filteredLow(:,  landmarkBins);
        filteredHigh_lm = filteredHigh(:, landmarkBins);

        speedModPerROI = mean(filteredLow_lm - filteredHigh_lm, 2, 'omitnan');
        allSpeedModSession{end+1} = speedModPerROI;
        fprintf('  Mean modulation: %.4f\n', mean(speedModPerROI, 'omitnan'));

        %% Shuffle test — shuffle lap speed labels, recompute modulation
        nLaps    = size(speedData, 1);
        nROIs_s  = length(FilteredROIS_idx);
        shuffleMod = nan(nROIs_s, nShuffles);

        fprintf('  Running %d shuffles...\n', nShuffles);
        for sh = 1:nShuffles
            shuffIdx  = randperm(nLaps);
            shuffLow  = nan(nROIs_s, length(landmarkBins));
            shuffHigh = nan(nROIs_s, length(landmarkBins));

            for bi = 1:length(landmarkBins)
                b         = landmarkBins(bi);
                binSpeeds = speedData(shuffIdx, b);
                lmMedian  = median(binSpeeds(binSpeeds > 1 & ~isnan(binSpeeds)));

                lowIdx  = binSpeeds > 1 & binSpeeds <= lmMedian & ~isnan(binSpeeds);
                highIdx = binSpeeds > lmMedian & ~isnan(binSpeeds);

                roiData = squeeze(rawROIData(FilteredROIS_idx, :, b)); % [nROIs x nLaps]

                if sum(lowIdx)  >= 2
                    shuffLow(:, bi)  = mean(roiData(:, lowIdx),  2, 'omitnan');
                end
                if sum(highIdx) >= 2
                    shuffHigh(:, bi) = mean(roiData(:, highIdx), 2, 'omitnan');
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



%% Plot per-position median lines across sessions 

%% Plot per-position median threshold lines across sessions

validRows = find(any(~isnan(allMedianThreshLines), 2));
threshMean = mean(allMedianThreshLines(validRows, :), 1, 'omitnan');
threshStd  = std(allMedianThreshLines(validRows, :), 0, 1, 'omitnan');

figure('Color', 'w', 'Position', [100 100 900 400]); hold on;

colors = lines(length(validRows));
for s = 1:length(validRows)
    plot(1:200, allMedianThreshLines(validRows(s), :), ...
        'Color', colors(s,:), 'LineWidth', 1.5);
end
plot(1:200, threshMean, 'k--', 'LineWidth', 2.5);

% Shade mean +/- std
fill([1:200, fliplr(1:200)], ...
     [threshMean + threshStd, fliplr(threshMean - threshStd)], ...
     'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Add landmark lines
for lm = [40, 80, 120, 160]
    xline(lm, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
end

xlabel('Position (cm)', 'FontName', 'Arial', 'FontSize', 12);
ylabel('Median speed threshold (cm/s)', 'FontName', 'Arial', 'FontSize', 12);
title(sprintf('Per-position median speed threshold across %d sessions', length(validRows)), ...
    'FontWeight', 'normal');

% Use allSessionLabels which is built correctly across all mice/sessions
legendLabels = [allSessionLabels(validRows), {'Mean ± SD'}];
legend(legendLabels, 'Box', 'off', 'Location', 'best', 'FontSize', 9);

set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial');
if exist('offsetAxes', 'file') == 2, offsetAxes(gca); end

fprintf('\nMedian threshold consistency across sessions:\n');
fprintf('  Mean CV across position bins: %.3f\n', ...
    mean(std(allMedianThreshLines(validRows,:), 0, 1, 'omitnan') ./ ...
         mean(allMedianThreshLines(validRows,:), 1, 'omitnan'), 'omitnan'));

%%
%% Test consistency of median speed threshold across sessions

validRows = find(any(~isnan(allMedianThreshLines), 2));
threshData = allMedianThreshLines(validRows, :);
nSessions  = length(validRows);

% 1 — Coefficient of variation per position bin
% Low CV = consistent across sessions
threshMean = mean(threshData, 1, 'omitnan');
threshStd  = std(threshData,  0, 1, 'omitnan');
threshCV   = threshStd ./ threshMean;

fprintf('\n=== Median Speed Threshold Consistency ===\n');
fprintf('  Mean CV across position bins: %.3f\n', mean(threshCV, 'omitnan'));
fprintf('  Max CV across position bins:  %.3f\n', max(threshCV, [], 'omitnan'));
fprintf('  (CV < 0.3 generally indicates acceptable consistency for pooling)\n');

% 2 — Kruskal-Wallis across position bins (are sessions different?)
% This tests: at each bin, do sessions have different medians?
pVals = nan(1, 200);
for b = 1:200
    binData = threshData(:, b);
    if sum(~isnan(binData)) >= 2
        groupLabel = (1:nSessions)';
        pVals(b) = kruskalwallis(binData, groupLabel, 'off');
    end
end

nSigBins = sum(pVals < 0.05, 'omitnan');
fprintf('  Significant position bins (KW): %d/200 (%.1f%%)\n', nSigBins, nSigBins/200*100);
fprintf('  (Expected ~5%% by chance if sessions are consistent)\n');

% 3 — ICC-style summary: correlation between sessions across position
% High correlation = sessions track same spatial speed profile
fprintf('\n  Pairwise correlations between session threshold lines:\n');
corrMat = corr(threshData', 'rows', 'pairwise');
offDiag = corrMat(tril(true(nSessions), -1));
fprintf('  Mean pairwise r = %.3f (range: %.3f to %.3f)\n', ...
    mean(offDiag), min(offDiag), max(offDiag));
fprintf('  (r > 0.8 suggests sessions are tracking the same speed profile)\n');

% 4 — Plot CV across position to see if inconsistency is localised
figure('Color', 'w', 'Position', [100 100 900 300]); hold on;
plot(1:200, threshCV, 'k-', 'LineWidth', 1.5);
yline(0.3, 'r--', 'LineWidth', 1.5);
for lm = [40, 80, 120, 160]
    xline(lm, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
end
xlabel('Position (cm)', 'FontName', 'Arial', 'FontSize', 12);
ylabel('Coefficient of Variation', 'FontName', 'Arial', 'FontSize', 12);
title('Median speed threshold consistency across sessions', 'FontWeight', 'normal');
set(gca, 'Box', 'off', 'TickDir', 'out');
if exist('offsetAxes', 'file') == 2, offsetAxes(gca); end
%% Per-session results with shuffle
fprintf('\nPer-session results:\n')
for s = 1:length(allSpeedModSession)
    if ~isempty(allSpeedModSession{s})
        sessionMod = allSpeedModSession{s};
        shuffP     = allShufflePSession{s};
        [~, p_t]   = ttest(sessionMod);
        [p_sr, ~]  = signrank(sessionMod);
        nSig       = sum(shuffP < 0.05);
        nSigPos    = sum(shuffP < 0.05 & sessionMod > 0);
        nSigNeg    = sum(shuffP < 0.05 & sessionMod < 0);
        nTotal     = length(sessionMod);
        nPos       = sum(sessionMod > 0);
        nNeg       = sum(sessionMod < 0);

        fprintf('  %s:\n', allSessionLabels{s})
        fprintf('    n=%d, mean=%.4f, median=%.4f\n', nTotal, mean(sessionMod,'omitnan'), median(sessionMod,'omitnan'))
        fprintf('    p(ttest)=%.4f, p(signrank)=%.4f\n', p_t, p_sr)
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

    % non-significant ROIs
    if sum(nonsigIdx) > 0
        histogram(sessionMod(nonsigIdx), 20, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7)
    end

    % significant ROIs — low > high
    sigPos = sessionMod(sigIdx & sessionMod > 0);
    if ~isempty(sigPos)
        histogram(sigPos, 20, 'FaceColor', [0.0 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.9)
    end

    % significant ROIs — high > low
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
    xlabel('Mean (Low - High)')
    ylabel('ROIs')
    legend('Non-sig', 'Sig low>high', 'Sig high>low', 'Box', 'off', 'FontSize', 7)
    set(gca, 'Box', 'off', 'TickDir', 'out')
end
sgtitle('Speed modulation per session — shuffle-significant ROIs highlighted')

%% Speed distributions per session
fprintf('\nLandmark speed distributions:\n')
for s = 1:length(allSessionLabels)
    rawData        = load(allFilePaths{s}, 'lapPositionRunningSpeed', 'trialIndicesByCondition');
    speedData      = rawData.lapPositionRunningSpeed;
    baselineTrials = rawData.trialIndicesByCondition.Baseline;
    speedData      = speedData(baselineTrials, :);
    fprintf('\n  %s:\n', allSessionLabels{s})
    for lm = 1:length(landmarks)
        window    = landmarks(lm)-15 : landmarks(lm)+15;
        lmSpeeds  = mean(speedData(:, window), 2);
        lmMedian  = median(lmSpeeds(lmSpeeds > 1 & ~isnan(lmSpeeds)));
        lowSpeeds  = lmSpeeds(lmSpeeds > 1 & lmSpeeds <= lmMedian);
        highSpeeds = lmSpeeds(lmSpeeds > lmMedian);
        fprintf('    LM %d cm: nLow=%d (%.1f cm/s), nHigh=%d (%.1f cm/s), contrast=%.1f cm/s\n', ...
                landmarks(lm), length(lowSpeeds), mean(lowSpeeds), ...
                length(highSpeeds), mean(highSpeeds), mean(highSpeeds)-mean(lowSpeeds))
    end
end

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

%% Speed ranges per session
fprintf('\nSpeed ranges per session:\n')
for s = 1:length(allSessionLabels)
    rawData        = load(allFilePaths{s}, 'lapPositionRunningSpeed', 'trialIndicesByCondition');
    speedData      = rawData.lapPositionRunningSpeed;
    baselineTrials = rawData.trialIndicesByCondition.Baseline;
    speedData      = speedData(baselineTrials, :);
    validSpeeds    = speedData(:);
    validSpeeds    = validSpeeds(validSpeeds > 1 & ~isnan(validSpeeds));
    meanLow        = mean(validSpeeds(validSpeeds <= allSessionMedians(s)));
    meanHigh       = mean(validSpeeds(validSpeeds >  allSessionMedians(s)));
    fprintf('  %s: low=%.1f cm/s, high=%.1f cm/s, MeanDiff=%.1f cm/s\n', ...
            allSessionLabels{s}, meanLow, meanHigh, meanHigh-meanLow)
end

%% Pool across sessions and test
iqrThreshold = 0;
includeIdx   = allSessionIQRs >= iqrThreshold;
fprintf('\nIncluding %d of %d sessions (IQR >= %.0f cm/s)\n', ...
        sum(includeIdx), length(allSessionIQRs), iqrThreshold)

allSpeedMod  = [];
allShuffleP  = [];
for s = 1:length(allSpeedModSession)
    if includeIdx(s) && ~isempty(allSpeedModSession{s})
        allSpeedMod  = [allSpeedMod;  allSpeedModSession{s}];
        allShuffleP  = [allShuffleP;  allShufflePSession{s}];
    end
end

[~, p_t, ~, stats] = ttest(allSpeedMod);
[p_sr, ~]          = signrank(allSpeedMod);
nSigAll            = sum(allShuffleP < 0.05);
nSigPosAll         = sum(allShuffleP < 0.05 & allSpeedMod > 0);
nSigNegAll         = sum(allShuffleP < 0.05 & allSpeedMod < 0);

fprintf('\nPooled result:\n')
fprintf('  n            = %d ROIs\n',  length(allSpeedMod))
fprintf('  mean         = %.4f\n',     mean(allSpeedMod))
fprintf('  median       = %.4f\n',     median(allSpeedMod))
fprintf('  t            = %.3f\n',     stats.tstat)
fprintf('  p (ttest)    = %.4f\n',     p_t)
fprintf('  p (signrank) = %.4f\n',     p_sr)
fprintf('  Shuffle sig  = %d/%d (%.1f%%) — low>high: %d, high>low: %d\n', ...
        nSigAll, length(allSpeedMod), nSigAll/length(allSpeedMod)*100, nSigPosAll, nSigNegAll)

%% Pooled histogram with shuffle-significant ROIs highlighted
figure('Position', [100 100 500 350], 'Color', 'w');
hold on

sigIdx   = allShuffleP < 0.05;
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

xline(0,                  '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off')
xline(mean(allSpeedMod),  '-k',  'LineWidth', 2,   'HandleVisibility', 'off')
xlabel('Mean (Low - High) \DeltaF/F')
ylabel('Number of ROIs')
legend('Non-sig', 'Sig low>high', 'Sig high>low', 'Box', 'off')
% title(sprintf('Pooled speed modulation\np(ttest)=%.4f, p(signrank)=%.4f, n=%d ROIs\nShuffle sig: %d/%d (low>high: %d, high>low: %d)', ...
      % p_t, p_sr, length(allSpeedMod), nSigAll, length(allSpeedMod), nSigPosAll, nSigNegAll))

title(sprintf('Pooled speed modulation | n=%d ROIs total | sig: %d/%d (%.1f%%) — low>high: %d, high>low: %d', ...
    length(allSpeedMod), nSigAll, length(allSpeedMod), ...
    nSigAll/length(allSpeedMod)*100, nSigPosAll, nSigNegAll))
set(gca, 'Box', 'off', 'TickDir', 'out')

%% Heatmap visualisation (last session loaded)

response = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_Response_M26003_BaselineCorridor_20260322_00002.mat")


filteredLow  = response.speedPositionActivity.lowHigh.matrixLow(FilteredROIS_idx, :);
filteredHigh = response.speedPositionActivity.lowHigh.matrixHigh(FilteredROIS_idx, :);
filteredDiff = response.speedPositionActivity.lowHigh.matrixDiff(FilteredROIS_idx, :);

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

sharedMax = prctile([smoothedLow(:); smoothedHigh(:)], 98);
sharedMin = prctile([smoothedLow(:); smoothedHigh(:)],  2);
diffMax   = prctile(abs(smoothedDiff(:)), 95);
x_pos     = 1:200;
nROIs_plot = size(sortedLow, 1);

figure('Position', [100 100 800 700], 'Color', 'w');

subplot(3,1,1)
imagesc(x_pos, 1:nROIs_plot, smoothedLow, [sharedMin sharedMax])
title('Low speed')
ylabel('ROI (sorted by pref. position (low)')
xline(landmarks, '--w', 'Alpha', 0.4, 'HandleVisibility', 'off')
colorbar; axis tight
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [])

subplot(3,1,2)
imagesc(x_pos, 1:nROIs_plot, smoothedHigh, [sharedMin sharedMax])
title('High speed')
ylabel('ROI (sorted by pref. position (low)')
xline(landmarks, '--w', 'Alpha', 0.4, 'HandleVisibility', 'off')
colorbar; axis tight
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [])

subplot(3,1,3)
imagesc(x_pos, 1:size(smoothedDiff,1), smoothedDiff, [-diffMax diffMax])
title('Low - High')
ylabel('ROI (sorted by pref. position (low))')
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
title('Population mean — Low vs High speed')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XLim', [1 200])

subplot(2,1,2)
hold on
plot(x_pos, meanDiff, 'k', 'LineWidth', 2)
yline(0, '--r', 'HandleVisibility', 'off')
xline(landmarks, '--k', 'HandleVisibility', 'off')
ylabel('Low - High')
xlabel('Position (cm)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XLim', [1 200])