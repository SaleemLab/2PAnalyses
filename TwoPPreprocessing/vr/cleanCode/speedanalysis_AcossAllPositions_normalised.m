%% Speed Modulation Analysis — Full Corridor, Normalized Index
% Same as the full-corridor version, but the modulation metric is now
% normalized: (low - high) / (low + high), computed at the SUMMARY level
% (summed across all position bins per ROI) rather than per-bin, to avoid
% the division instability that occurs when DF/F approaches zero at
% individual off-peak bins.
%
% IMPORTANT: the shuffle-test null distribution is now built on this same
% normalized metric, not the raw difference. If you only swap the real
% statistic and leave the shuffle null on raw differences, the resulting
% p-values compare two different quantities and are meaningless.
%
% ROIs whose |sumLow + sumHigh| falls below denomFloor are flagged as
% unreliable rather than silently included — check unreliableCount before
% trusting a session's results.

%% Session pairs
pairs        = struct;
pairs.M25132 = {'20260226', '20260228', '20260313'};
pairs.M26003 = {'20260322', '20260324', '20260325'};

filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);
mouseInfo     = sessionsToProcess(filteredTable);

%% Initialise
allSessionMedians    = [];
allSessionIQRs       = [];
allSessionLabels     = {};
allFilePaths         = {};
allSpeedModSession    = {}; % now holds the NORMALIZED index per ROI
allUnreliableSession  = {}; % logical mask, per ROI, per session
allShufflePSession    = {};
allMedianThreshLines  = [];

landmarks    = [40, 80, 120, 160]; % reference lines only
nShuffles    = 1000;

minLapsPerGroup      = 20;
minBinFracSufficient = 0.8;

% Denominator floor for the normalized index — a session's sumLow+sumHigh
% below this magnitude is flagged unreliable rather than trusted at face
% value. TUNE THIS against your actual DF/F scale; there's no principled
% default. Starting point: look at the distribution of |sumLow+sumHigh|
% for a session and set this near its 5th percentile.
denomFloor = 1.0;

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
        fprintf('  Processing Session: %s\n', sessionName);

        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if ~isfile(infoPath), fprintf('  Info Missing — skipping\n'); continue; end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};

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
        filePath   = sessionFileInfo.stimFiles(targetIdx).Response;
        allFilePaths{end+1} = filePath;

        rawData        = load(filePath, 'lapPositionRunningSpeed', 'lapPositionActivity', 'trialIndicesByCondition', 'SMI_Metrics');
        speedData      = rawData.lapPositionRunningSpeed;
        baselineTrials = rawData.trialIndicesByCondition.Baseline;
        speedData      = speedData(baselineTrials, :);
        rawROIData     = rawData.lapPositionActivity.dFFNeuropilCorrected;
        rawROIData     = rawROIData(:, baselineTrials, :);

        numPosBins = size(speedData, 2);

        validSpeeds   = speedData(:);
        validSpeeds   = validSpeeds(validSpeeds > 1 & ~isnan(validSpeeds));
        sessionMedian = median(validSpeeds);
        sessionIQR    = iqr(validSpeeds);

        allSessionMedians(end+1) = sessionMedian;
        allSessionIQRs(end+1)    = sessionIQR;
        allSessionLabels{end+1}  = sprintf('%s_%s', mousenumber, sessionName);
        fprintf('  Median: %.1f cm/s | IQR: %.1f cm/s\n', sessionMedian, sessionIQR);

        %% Lap-sufficiency check across the whole corridor
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
            allSpeedModSession{end+1}   = [];
            allUnreliableSession{end+1} = [];
            allShufflePSession{end+1}   = [];
            continue
        end
        fprintf('  %.0f%% of position bins have sufficient laps per group\n', fracSufficient*100)

        %% Low/high speed position matrices
        options.signalToUse = 'dFFNeuropilCorrected';
        [response] = getLowHighSpeedPositionMatrix(sessionFileInfo, VRStimName, options);

        matOptions.signalToUse    = 'dFFNeuropilCorrected';
        matOptions.numBins        = 10;
        matOptions.applySmoothing = false;
        getPositionSpeedMatrix(sessionFileInfo, VRStimName, matOptions);

        allMedianThreshLines(end+1, :) = response.speedPositionActivity.lowHigh.medianThreshLine;

        %% ROI filters
        rawData2       = load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'crossValExpVar');
        crossValExpVar = rawData2.crossValExpVar;

        FilteredROIS_idx = find( ...
            crossValExpVar.dFFNeuropilCorrected.pValues    <= 0.01 & ...
            crossValExpVar.dFFNeuropilCorrected.medianExpVar > 0.1 & ...
            ~rawData.SMI_Metrics.dFFNeuropilCorrected.ExcludeEdgePeakCells);

        if isempty(FilteredROIS_idx)
            fprintf('  No valid ROIs — skipping\n');
            allSpeedModSession{end+1}   = [];
            allUnreliableSession{end+1} = [];
            allShufflePSession{end+1}   = [];
            continue
        end
        fprintf('  nROIs: %d\n', length(FilteredROIS_idx));

        filteredLow  = response.speedPositionActivity.lowHigh.matrixLow(FilteredROIS_idx, :);
        filteredHigh = response.speedPositionActivity.lowHigh.matrixHigh(FilteredROIS_idx, :);
        filteredDiff = filteredLow - filteredHigh;

        %% Normalized modulation index — summary level
        sumLow  = sum(filteredLow,  2, 'omitnan');
        sumHigh = sum(filteredHigh, 2, 'omitnan');
        denom   = sumLow + sumHigh;

        speedModPerROI = (sumLow - sumHigh) ./ denom;
        unreliableROI  = abs(denom) < denomFloor;

        allSpeedModSession{end+1}   = speedModPerROI;
        allUnreliableSession{end+1} = unreliableROI;

        nUnreliable = sum(unreliableROI);
        fprintf('  Normalized modulation (full corridor): mean=%.4f, median=%.4f | %d/%d ROIs flagged unreliable (|denom|<%.2f)\n', ...
                mean(speedModPerROI(~unreliableROI), 'omitnan'), ...
                median(speedModPerROI(~unreliableROI), 'omitnan'), ...
                nUnreliable, length(speedModPerROI), denomFloor);

        lastGood.filteredLow  = filteredLow;
        lastGood.filteredHigh = filteredHigh;
        lastGood.filteredDiff = filteredDiff;
        lastGood.label        = sprintf('%s_%s', mousenumber, sessionName);

        %% Shuffle test — null built on the SAME normalized metric
        nLaps    = size(speedData, 1);
        nROIs_s  = length(FilteredROIS_idx);
        shuffleMod = nan(nROIs_s, nShuffles);

        fprintf('  Running %d shuffles across %d bins (normalized metric)...\n', nShuffles, numPosBins);
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

                roiData = squeeze(rawROIData(FilteredROIS_idx, :, b));

                if sum(lowIdx)  >= 2
                    shuffLow(:, b)  = mean(roiData(:, lowIdx),  2, 'omitnan');
                end
                if sum(highIdx) >= 2
                    shuffHigh(:, b) = mean(roiData(:, highIdx), 2, 'omitnan');
                end
            end

            shuffSumLow  = sum(shuffLow,  2, 'omitnan');
            shuffSumHigh = sum(shuffHigh, 2, 'omitnan');
            shuffleMod(:, sh) = (shuffSumLow - shuffSumHigh) ./ (shuffSumLow + shuffSumHigh);
        end

        shuffP = nan(nROIs_s, 1);
        for r = 1:nROIs_s
            shuffP(r) = mean(abs(shuffleMod(r,:)) >= abs(speedModPerROI(r)));
        end
        % Unreliable ROIs get their p-value set to NaN rather than a
        % number that looks trustworthy — exclude them from sig counts.
        shuffP(unreliableROI) = NaN;
        allShufflePSession{end+1} = shuffP;

        nSig    = sum(shuffP < 0.05, 'omitnan');
        nSigPos = sum(shuffP < 0.05 & speedModPerROI > 0, 'omitnan');
        nSigNeg = sum(shuffP < 0.05 & speedModPerROI < 0, 'omitnan');
        fprintf('  Shuffle: %d/%d sig (%.1f%%) — low>high: %d, high>low: %d [excluding %d unreliable]\n', ...
                nSig, nROIs_s - nUnreliable, nSig/(nROIs_s - nUnreliable)*100, nSigPos, nSigNeg, nUnreliable)
    end
end

%% Per-session results
fprintf('\nPer-session results (normalized, full corridor):\n')
for s = 1:length(allSpeedModSession)
    if ~isempty(allSpeedModSession{s})
        sessionMod   = allSpeedModSession{s};
        unreliable   = allUnreliableSession{s};
        shuffP       = allShufflePSession{s};
        reliableMask = ~unreliable;

        [~, p_t]  = ttest(sessionMod(reliableMask));
        [p_sr, ~] = signrank(sessionMod(reliableMask));
        nSig      = sum(shuffP < 0.05, 'omitnan');
        nSigPos   = sum(shuffP < 0.05 & sessionMod > 0, 'omitnan');
        nSigNeg   = sum(shuffP < 0.05 & sessionMod < 0, 'omitnan');
        nTotal    = sum(reliableMask);
        nPos      = sum(sessionMod(reliableMask) > 0);
        nNeg      = sum(sessionMod(reliableMask) < 0);

        fprintf('  %s:\n', allSessionLabels{s})
        fprintf('    n=%d (reliable), %d flagged unreliable\n', nTotal, sum(unreliable))
        fprintf('    mean=%.4f, median=%.4f\n', mean(sessionMod(reliableMask),'omitnan'), median(sessionMod(reliableMask),'omitnan'))
        fprintf('    p(ttest)=%.4f, p(signrank)=%.4f\n', p_t, p_sr)
        fprintf('    low>high: %d (%.1f%%), high>low: %d (%.1f%%)\n', nPos, nPos/nTotal*100, nNeg, nNeg/nTotal*100)
        fprintf('    Shuffle sig: %d/%d (%.1f%%) — low>high: %d, high>low: %d\n', nSig, nTotal, nSig/nTotal*100, nSigPos, nSigNeg)
    end
end

%% Pool across sessions
allSpeedMod    = [];
allShuffleP    = [];
allUnreliable  = [];
for s = 1:length(allSpeedModSession)
    if ~isempty(allSpeedModSession{s})
        allSpeedMod   = [allSpeedMod;   allSpeedModSession{s}];
        allShuffleP   = [allShuffleP;   allShufflePSession{s}];
        allUnreliable = [allUnreliable; allUnreliableSession{s}];
    end
end
reliableMask = ~allUnreliable;

[~, p_t, ~, stats] = ttest(allSpeedMod(reliableMask));
[p_sr, ~]           = signrank(allSpeedMod(reliableMask));
nSigAll             = sum(allShuffleP < 0.05, 'omitnan');
nSigPosAll          = sum(allShuffleP < 0.05 & allSpeedMod > 0, 'omitnan');
nSigNegAll          = sum(allShuffleP < 0.05 & allSpeedMod < 0, 'omitnan');
nReliable           = sum(reliableMask);

fprintf('\nPooled result (normalized, full corridor):\n')
fprintf('  n            = %d ROIs (%d flagged unreliable, excluded)\n', nReliable, sum(allUnreliable))
fprintf('  mean         = %.4f\n',     mean(allSpeedMod(reliableMask)))
fprintf('  median       = %.4f\n',     median(allSpeedMod(reliableMask)))
fprintf('  t            = %.3f\n',     stats.tstat)
fprintf('  p (ttest)    = %.4f\n',     p_t)
fprintf('  p (signrank) = %.4f\n',     p_sr)
fprintf('  Shuffle sig  = %d/%d (%.1f%%) — low>high: %d, high>low: %d\n', ...
        nSigAll, nReliable, nSigAll/nReliable*100, nSigPosAll, nSigNegAll)

%% Pooled histogram
figure('Position', [100 100 500 350], 'Color', 'w');
hold on
sigIdx    = allShuffleP < 0.05 & reliableMask;
nonsigIdx = reliableMask & ~sigIdx;

histogram(allSpeedMod(nonsigIdx), 50, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7)
sigPos = allSpeedMod(sigIdx & allSpeedMod > 0);
sigNeg = allSpeedMod(sigIdx & allSpeedMod < 0);
if ~isempty(sigPos)
    histogram(sigPos, 50, 'FaceColor', [0.0 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.9)
end
if ~isempty(sigNeg)
    histogram(sigNeg, 50, 'FaceColor', [0.8 0.0 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.9)
end
xline(0, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off')
xlabel('(Low - High) / (Low + High), full corridor, summary-level')
ylabel('Number of ROIs')
legend('Non-sig', 'Sig low>high', 'Sig high>low', 'Box', 'off')
title(sprintf('Pooled normalized speed modulation | n=%d reliable ROIs | sig: %d/%d (%.1f%%)', ...
    nReliable, nSigAll, nReliable, nSigAll/nReliable*100))
set(gca, 'Box', 'off', 'TickDir', 'out')

%% Example-session heatmap (unchanged — still shows raw DF/F, not the index)
if ~isempty(lastGood.filteredLow)
    filteredLow  = lastGood.filteredLow;
    filteredHigh = lastGood.filteredHigh;
    filteredDiff = lastGood.filteredDiff;

    [~, peakPos]   = max(filteredLow, [], 2);
    [~, sortOrder] = sort(peakPos);

    sortedLow  = filteredLow(sortOrder, :);
    sortedHigh = filteredHigh(sortOrder, :);
    sortedDiff = filteredDiff(sortOrder, :);

    sigmaSpatial = 1; sigmaROI = 2;
    smoothedLow  = imgaussfilt(sortedLow,  [sigmaROI sigmaSpatial]);
    smoothedHigh = imgaussfilt(sortedHigh, [sigmaROI sigmaSpatial]);
    smoothedDiff = imgaussfilt(sortedDiff, [sigmaROI sigmaSpatial]);

    sharedMax = prctile([smoothedLow(:); smoothedHigh(:)], 98);
    sharedMin = prctile([smoothedLow(:); smoothedHigh(:)],  2);
    diffMax   = prctile(abs(smoothedDiff(:)), 95);
    x_pos     = 1:size(filteredLow, 2);
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
    title('Low - High (raw, unnormalized)')
    ylabel('ROI (sorted by pref. position, low)')
    xlabel('Position (cm)')
    xline(landmarks, '--w', 'Alpha', 0.4, 'HandleVisibility', 'off')
    colorbar; axis tight
    set(gca, 'Box', 'off', 'TickDir', 'out')
else
    fprintf('\nNo session data available for example-session plots.\n')
end