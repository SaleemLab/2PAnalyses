%% plotBaselineSubtractionDemo_script.m
% Plain script version -- edit the params below and just run it (F5).
% Loads one session's response, picks one ROI, and plots OLD (track-wide
% min/max) vs NEW (local-baseline-subtracted) normalization side by side.

%% ---- EDIT THESE ----
VRStimName          = 'M26004_LandManipCorridor_20260318_CombinedRuns';   % <-- set this
roiIdx              = [];              
excludeFlaggedLaps  = true;
applySmoothing      = true;
baselineMethod      = 'prctile5';       % 'prctile5' or 'median'
% sessionFileInfo must already exist in your workspace
% -----------------------

stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName ''%s'' not found.', VRStimName); end

response = load(sessionFileInfo.stimFiles(stimIdx).Response, ...
    'lapPositionActivity', 'trialIndicesByCondition', 'stimName', 'flaggedLaps');

partnerLandmarkPosition = 80;
excludeStart = 30; excludeEnd = 30;
landmarkCentres = [40, 80, 120, 160];
tolerance = 10;

signalNames = fieldnames(response.lapPositionActivity);
currentSignalName = signalNames{1};   % change index if you want a different signal type
fprintf('Using signal type: %s\n', currentSignalName);

lapActivity = response.lapPositionActivity.(currentSignalName);
baseTrials = response.trialIndicesByCondition.Baseline;

if excludeFlaggedLaps && isfield(response, 'flaggedLaps') && ~isempty(response.flaggedLaps)
    baseTrials = setdiff(baseTrials, response.flaggedLaps);
end

baseLapActivity = lapActivity(:, baseTrials, :);
[numROIs, numLaps, numBins] = size(baseLapActivity);

allowedLandmarkBins = [];
for c = landmarkCentres
    allowedLandmarkBins = [allowedLandmarkBins, (c - tolerance):(c + tolerance)]; %#ok<AGROW>
end
allowedLandmarkBins = unique(allowedLandmarkBins);
validInnerRange = (excludeStart + 1) : (numBins - excludeEnd);
allowedSearchBins = intersect(allowedLandmarkBins, validInnerRange);
baselineBins = setdiff(validInnerRange, allowedLandmarkBins);

if applySmoothing
    w = gausswin(15); w = w / sum(w);
    for iCell = 1:numROIs
        for iLap = 1:numLaps
            trace = squeeze(baseLapActivity(iCell, iLap, :));
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace); trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace); smoothed(nanMask) = NaN;
            baseLapActivity(iCell, iLap, :) = smoothed;
        end
    end
end

meanOdd  = squeeze(mean(baseLapActivity(:, 1:2:end, :), 2, 'omitnan'));
meanEven = squeeze(mean(baseLapActivity(:, 2:2:end, :), 2, 'omitnan'));

%% Auto-pick ROI if none given (biggest old-normalization peak-vs-baseline gap)
if isempty(roiIdx)
    bestGap = -Inf; bestROI = 1;
    for r = 1:numROIs
        trainTrace = meanOdd(r, :);
        if all(isnan(trainTrace)), continue; end
        minOdd = min(trainTrace, [], 'omitnan'); maxOdd = max(trainTrace, [], 'omitnan');
        rangeOdd = maxOdd - minOdd; if rangeOdd == 0, rangeOdd = 1; end
        normTrain = (trainTrace - minOdd) ./ rangeOdd;
        [~, mIdx] = max(normTrain(allowedSearchBins));
        pBin = allowedSearchBins(mIdx);
        gap = abs(normTrain(pBin) - median(normTrain(baselineBins), 'omitnan'));
        if gap > bestGap, bestGap = gap; bestROI = r; end
    end
    roiIdx = bestROI;
end
fprintf('Plotting ROI #%d\n', roiIdx);

trainTrace = meanOdd(roiIdx, :);
testTrace  = meanEven(roiIdx, :);

%% ---- OLD normalization: track-wide min/max ----
minOdd = min(trainTrace, [], 'omitnan');
maxOdd = max(trainTrace, [], 'omitnan');
rangeOdd = maxOdd - minOdd; if rangeOdd == 0, rangeOdd = 1; end
normTrain_old = (trainTrace - minOdd) ./ rangeOdd;
normTest_old  = (testTrace - minOdd) ./ rangeOdd;

[~, maxIdxInSearch] = max(normTrain_old(allowedSearchBins));
prefBin = allowedSearchBins(maxIdxInSearch);
if prefBin <= 100
    partnerBin = prefBin + partnerLandmarkPosition;
else
    partnerBin = prefBin - partnerLandmarkPosition;
end
partnerBin = min(max(partnerBin, 1), numBins);

Rp_old = normTest_old(prefBin);
Rn_old = normTest_old(partnerBin);
SMI_old = NaN;
if (Rp_old + Rn_old) ~= 0
    SMI_old = (Rp_old - Rn_old) / (Rp_old + Rn_old);
end

%% ---- NEW normalization: local-baseline-subtracted ----
if strcmp(baselineMethod, 'median')
    baselineVal = median(trainTrace(baselineBins), 'omitnan');
else
    baselineVal = prctile(trainTrace(baselineBins), 5);
end
trainTrace_bs = trainTrace - baselineVal;
testTrace_bs  = testTrace  - baselineVal;

minOdd_bs = min(trainTrace_bs, [], 'omitnan');
maxOdd_bs = max(trainTrace_bs, [], 'omitnan');
rangeOdd_bs = maxOdd_bs - minOdd_bs; if rangeOdd_bs == 0, rangeOdd_bs = 1; end
normTrain_new = (trainTrace_bs - minOdd_bs) ./ rangeOdd_bs;
normTest_new  = (testTrace_bs - minOdd_bs) ./ rangeOdd_bs;

[~, maxIdxInSearch_new] = max(normTrain_new(allowedSearchBins));
prefBin_new = allowedSearchBins(maxIdxInSearch_new);
if prefBin_new <= 100
    partnerBin_new = prefBin_new + partnerLandmarkPosition;
else
    partnerBin_new = prefBin_new - partnerLandmarkPosition;
end
partnerBin_new = min(max(partnerBin_new, 1), numBins);

Rp_new = normTest_new(prefBin_new);
Rn_new = normTest_new(partnerBin_new);
SMI_new = NaN;
if (Rp_new + Rn_new) ~= 0
    SMI_new = (Rp_new - Rn_new) / (Rp_new + Rn_new);
end

%% ---- Plot ----
figure('Name', sprintf('Baseline subtraction demo - ROI %d', roiIdx), 'Color', [1 1 1], 'Position', [100 100 1100 700]);

subplot(2,2,1); hold on;
plot(normTest_old, 'k', 'LineWidth', 1.5);
plot(prefBin, Rp_old, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(partnerBin, Rn_old, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
for c = landmarkCentres, xline(c, ':', 'Color', [0.6 0.6 0.6]); end
title(sprintf('OLD norm (track-wide min/max)\nSMI = %.3f', SMI_old));
xlabel('position bin'); ylabel('normalized activity (test/even)');
legend({'test trace', 'preferred (Rp)', 'partner/non-pref (Rn)'}, 'Location', 'best');
ylim([-0.2 1.2]); grid on;

subplot(2,2,2); hold on;
plot(normTest_new, 'k', 'LineWidth', 1.5);
plot(prefBin_new, Rp_new, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(partnerBin_new, Rn_new, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
for c = landmarkCentres, xline(c, ':', 'Color', [0.6 0.6 0.6]); end
title(sprintf('NEW norm (baseline-subtracted, %s)\nSMI = %.3f', baselineMethod, SMI_new));
xlabel('position bin'); ylabel('normalized activity (test/even)');
legend({'test trace', 'preferred (Rp)', 'partner/non-pref (Rn)'}, 'Location', 'best');
ylim([-0.2 1.2]); grid on;

subplot(2,2,3); hold on;
plot(trainTrace, 'b', 'LineWidth', 1.5, 'DisplayName', 'raw train trace (odd laps)');
yline(baselineVal, '--', 'Color', [0.8 0.3 0.3], 'DisplayName', sprintf('baseline (%s) = %.3f', baselineMethod, baselineVal));
for c = landmarkCentres, xline(c, ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off'); end
scatter(baselineBins, trainTrace(baselineBins), 8, [0.7 0.7 0.9], 'filled', 'DisplayName', 'inter-landmark bins used for baseline');
title('raw train trace + estimated local baseline');
xlabel('position bin'); ylabel('raw activity (train/odd)');
legend('Location', 'best'); grid on;

subplot(2,2,4); hold on;
bar([1 2], [SMI_old, SMI_new], 0.5, 'FaceColor', [0.5 0.5 0.5]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'old SMI', 'new SMI'});
ylabel('SMI value'); ylim([-1.1, 1.1]); yline(0, 'k-');
title('SMI comparison for this ROI'); grid on;

fprintf('\nROI %d summary:\n', roiIdx);
fprintf('  OLD: prefBin=%d, partnerBin=%d, Rp=%.3f, Rn=%.3f, SMI=%.3f\n', prefBin, partnerBin, Rp_old, Rn_old, SMI_old);
fprintf('  NEW: prefBin=%d, partnerBin=%d, Rp=%.3f, Rn=%.3f, SMI=%.3f, baseline=%.3f\n', prefBin_new, partnerBin_new, Rp_new, Rn_new, SMI_new, baselineVal);