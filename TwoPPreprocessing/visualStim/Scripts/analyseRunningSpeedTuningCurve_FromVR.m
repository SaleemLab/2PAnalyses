%% ===== TEST: running speed tuning in VR, one session, corridor runs combined (nothing saved) =====
mouse = 'M26003';  sess = '20260322';
signalName = 'dFFNeuropilCorrected';
useZScoredSignals = false;
p_thresh = 0.01;

sfi = load(findSessionFileInfoFilePath(mouse, sess), 'sessionFileInfo');  sessionFileInfo = sfi.sessionFileInfo;
stimNames = {sessionFileInfo.stimFiles.name};
corrIdx   = find(contains(stimNames, 'BaselineCorridor'));
fprintf('Corridor runs: %s\n', strjoin(stimNames(corrIdx), ', '));

sourceStructName = 'processedSignals';
if useZScoredSignals, sourceStructName = 'zScoredProcessedSignals'; end
tickToCm = 3.1415 * 20 / 1024;                       % same as getRunningSpeedTuningCurves

sigAll = [];  speedAll = [];
for k = corrIdx
    % Same files as getRunningSpeedTuningCurves
    data2P     = load(sessionFileInfo.stimFiles(k).processedMergedBonsaiSuite2pData, sourceStructName);
    dataPeriph = load(sessionFileInfo.stimFiles(k).processedPeripheralData, 'peripheralData');
    sigData = data2P.(sourceStructName).(signalName);

    % Wheel speed, exactly as in getRunningSpeedTuningCurves (computed per run)
    w = dataPeriph.peripheralData.Wheel;
    wheelSpeed = [0; diff(w.Value * tickToCm)] ./ [1; diff(w.sampleTimes)];
    wheelSpeed(abs(wheelSpeed) > 150) = NaN;
    if size(sigData, 2) ~= numel(wheelSpeed)
        error('%s: %d signal frames vs %d wheel samples', stimNames{k}, size(sigData, 2), numel(wheelSpeed));
    end

    % Keep within-lap frames only (from each lap's first to last imaging frame)
    resp = load(sessionFileInfo.stimFiles(k).Response, 'lapPosition2PFrameIdx');
    inLap = false(1, numel(wheelSpeed));
    for lap = 1:size(resp.lapPosition2PFrameIdx, 1)
        fr = cellfun(@(x) x(:), resp.lapPosition2PFrameIdx(lap, :), 'UniformOutput', false);
        fr = vertcat(fr{:});
        if ~isempty(fr), inLap(min(fr):max(fr)) = true; end
    end
    fprintf('  %s: %d of %d frames within laps\n', stimNames{k}, sum(inLap), numel(inLap));

    sigAll   = [sigAll,   sigData(:, inLap)];        %#ok<AGROW>
    speedAll = [speedAll; wheelSpeed(inLap)];        %#ok<AGROW>
end
clear data2P dataPeriph sigData

%% Tuning curves: core of getRunningSpeedTuningCurves, unchanged
response = struct();
response.stimName = 'VRSpeedTuning_test';
response.tuningCurve = computeSpeedTuningCurve(sigAll, speedAll, signalName, sourceStructName);

tc = response.tuningCurve;
fprintf('\nVR: %d speed bins, edges %s cm/s\n', numel(tc.speedBins) - 1, mat2str(round(tc.speedBins, 1)));
fprintf('Running time per bin (s): %s | stationary: %.0f s\n', mat2str(round(tc.occupancy.moving)), tc.occupancy.stationary);

%% Classify with v5 (no saving)
response = classifySpeedTuningFromCurves(sessionFileInfo, response, signalName, false);

pv  = response.tuningCurve.(signalName).pValMoving;
cls = response.tuningCurve.(signalName).classification;
sig = pv <= p_thresh;
fprintf('\nSpeed-tuned in VR (p <= %.2f): %d of %d (%.1f%%)\n', p_thresh, sum(sig), numel(sig), 100*mean(sig));
types = {'lowpass', 'highpass', 'bandpass', 'trough_inverted'};
for t = 1:4
    n = sum(strcmp(cls.tuningType(sig), types{t}));
    fprintf('  %-16s %4d (%.1f%%)\n', types{t}, n, 100*n/max(sum(sig), 1));
end

%% Example tuning curves and fits (your existing plotting function)
plotAllSpeedTuningCategories(response, signalName);



%%
%% Compare VR with darkness, same boutons
stimNames = {sessionFileInfo.stimFiles.name};
darkIdx = find(contains(stimNames, 'Darkness'), 1);
dk = load(sessionFileInfo.stimFiles(darkIdx).Response, 'response');  dk = dk.response;

pD = dk.tuningCurve.(signalName).pValMoving;  clsD = dk.tuningCurve.(signalName).classification;
pV = response.tuningCurve.(signalName).pValMoving;  clsV = response.tuningCurve.(signalName).classification;
if numel(pD) ~= numel(pV), error('Different numbers of boutons in darkness (%d) and VR (%d).', numel(pD), numel(pV)); end

sigD = pD <= p_thresh;  sigV = pV <= p_thresh;
fprintf('\nDarkness: %d of %d tuned (%.1f%%) | VR: %d (%.1f%%)\n', sum(sigD), numel(sigD), 100*mean(sigD), sum(sigV), 100*mean(sigV));
fprintf('Tuned in both: %d | darkness only: %d | VR only: %d | neither: %d\n', ...
    sum(sigD & sigV), sum(sigD & ~sigV), sum(~sigD & sigV), sum(~sigD & ~sigV));

% Shape agreement for boutons tuned in both
types = {'lowpass', 'highpass', 'bandpass', 'trough_inverted'};
both = find(sigD & sigV);
T = zeros(4);
for i = both'
    a = find(strcmp(types, clsD.tuningType{i}));  b = find(strcmp(types, clsV.tuningType{i}));
    if ~isempty(a) && ~isempty(b), T(a, b) = T(a, b) + 1; end
end
disp(array2table(T, 'RowNames', strcat('dark_', types), 'VariableNames', strcat('VR_', types)));
fprintf('Same shape in both: %d of %d (%.0f%%)\n', trace(T), sum(T(:)), 100*trace(T)/max(sum(T(:)), 1));



%%
%% Superimpose darkness and VR tuning curves for the same boutons
tcD = dk.tuningCurve;        tcV = response.tuningCurve;
eD  = tcD.speedBins(:)';     eV  = tcV.speedBins(:)';
cD  = eD(1:end-1) + diff(eD)/2;   cV = eV(1:end-1) + diff(eV)/2;     % bin centres (cm/s)
gfun = @(p, x) p(1) + p(2) * exp(-((x - p(3)).^2) / (2 * p(4)^2));

groups = { find(sigD & sigV), find(~sigD & sigV), find(sigD & ~sigV) };
groupNames = {'Tuned in both', 'Tuned in VR only', 'Tuned in darkness only'};
nPer = 6;
colD = [0.2 0.2 0.2];  colV = [0.85 0.33 0.10];

figEx = figure('Color', 'w', 'Position', [50 50 1300 650]);
for g = 1:3
    idx = groups{g};
    % Order by strength of tuning (smallest p in the condition defining the group)
    if g == 3, [~, o] = sort(pD(idx)); else, [~, o] = sort(pV(idx)); end
    idx = idx(o(1:min(nPer, numel(idx))));
    for k = 1:numel(idx)
        r  = idx(k);
        ax = subplot(3, nPer, (g-1)*nPer + k); hold(ax, 'on');

        % Darkness
        errorbar(ax, cD, tcD.(signalName).moveMean(r,:), tcD.(signalName).moveSEM(r,:), 'o', ...
            'Color', colD, 'MarkerFaceColor', colD, 'MarkerSize', 3, 'CapSize', 0);
        plot(ax, 0, tcD.(signalName).statMean(r), 's', 'Color', colD, 'MarkerFaceColor', 'w', 'MarkerSize', 4);
        pF = clsD.fitParams(r, :);
        if all(isfinite(pF))
            xx = linspace(cD(1), cD(end), 100);  plot(ax, xx, gfun(pF, xx), '-', 'Color', colD, 'LineWidth', 1.2);
        end

        % VR
        errorbar(ax, cV, tcV.(signalName).moveMean(r,:), tcV.(signalName).moveSEM(r,:), 'o', ...
            'Color', colV, 'MarkerFaceColor', colV, 'MarkerSize', 3, 'CapSize', 0);
        plot(ax, 0, tcV.(signalName).statMean(r), 's', 'Color', colV, 'MarkerFaceColor', 'w', 'MarkerSize', 4);
        pF = clsV.fitParams(r, :);
        if all(isfinite(pF))
            xx = linspace(cV(1), cV(end), 100);  plot(ax, xx, gfun(pF, xx), '-', 'Color', colV, 'LineWidth', 1.2);
        end

        title(ax, sprintf('ROI %d\nD: %s (p=%.3f)\nVR: %s (p=%.3f)', r, clsD.tuningType{r}, pD(r), ...
            clsV.tuningType{r}, pV(r)), 'FontSize', 7, 'FontWeight', 'normal', 'Interpreter', 'none');
        set(ax, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 7);
        if k == 1, ylabel(ax, sprintf('%s\n\\DeltaF/F', groupNames{g}), 'FontSize', 8); end
        if g == 3, xlabel(ax, 'Speed (cm/s)', 'FontSize', 8); end
    end
end
sgtitle(figEx, 'Darkness (black) vs VR (orange); squares = stationary', 'FontSize', 10);


%%
%% Fair comparison: VR tuning curves on the darkness bins, speed index scatter
eD = dk.tuningCurve.speedBins(:)';  nB = numel(eD) - 1;
mD = dk.tuningCurve.(signalName).moveMean;                    % darkness curves (ROIs x bins)

% VR tuning curves on the darkness bins (frames outside the darkness range ignored)
mV = nan(size(mD));
for b = 1:nB
    idx = speedAll >= eD(b) & speedAll < eD(b+1);
    mV(:, b) = mean(sigAll(:, idx), 2, 'omitnan');
    fprintf('Bin %d (%.1f-%.1f cm/s): %.0f s in VR\n', b, eD(b), eD(b+1), sum(idx)/60);
end

% Speed index: Spearman correlation between bin number and response
binRank = (1:nB)';
siD = corr(binRank, mD', 'type', 'Spearman')';
siV = corr(binRank, mV', 'type', 'Spearman')';

% Across boutons
[rAll, pAll] = corr(siD, siV, 'type', 'Spearman', 'rows', 'complete');
fprintf('\nSpeed index, darkness vs VR (all %d boutons): rho = %.2f, p = %.3g\n', numel(siD), rAll, pAll);
tD = find(sigD);
fprintf('Boutons tuned in darkness (n=%d): same sign in VR for %d (%.0f%%)\n', ...
    numel(tD), sum(sign(siD(tD)) == sign(siV(tD))), 100*mean(sign(siD(tD)) == sign(siV(tD))));

% Tuning-curve correlation: own bouton vs random pairing
ccOwn = arrayfun(@(r) corr(mD(r,:)', mV(r,:)'), (1:size(mD,1))');
rng(1);  perm = randperm(size(mD,1));
ccNull = arrayfun(@(r) corr(mD(r,:)', mV(perm(r),:)'), (1:size(mD,1))');
fprintf('Tuning-curve correlation: own bouton median %.2f vs random pairing %.2f (rank-sum p = %.3g)\n', ...
    median(ccOwn, 'omitnan'), median(ccNull, 'omitnan'), ranksum(ccOwn, ccNull));

%% Figure
figFair = figure('Color', 'w', 'Position', [100 100 700 300]);
ax1 = subplot(1,2,1); hold(ax1, 'on');
scatter(ax1, siD, siV, 12, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.5);
scatter(ax1, siD(tD), siV(tD), 25, [0.85 0.33 0.10], 'filled');
plot(ax1, [-1 1], [-1 1], 'k--');  xline(ax1, 0, ':');  yline(ax1, 0, ':');
axis(ax1, 'square');  xlim(ax1, [-1 1]);  ylim(ax1, [-1 1]);
xlabel(ax1, 'Speed index, darkness');  ylabel(ax1, 'Speed index, VR (darkness bins)');
title(ax1, sprintf('\\rho = %.2f', rAll), 'FontWeight', 'normal');
legend(ax1, {'All boutons', 'Tuned in darkness'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);

ax2 = subplot(1,2,2); hold(ax2, 'on');
[f1, x1] = ecdf(ccNull);  plot(ax2, x1, f1, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
[f2, x2] = ecdf(ccOwn);   plot(ax2, x2, f2, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
xlabel(ax2, 'Darkness vs VR tuning-curve correlation');  ylabel(ax2, 'Cumulative fraction');
legend(ax2, {'Random pairing', 'Same bouton'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
set([ax1 ax2], 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);



%% ---------- local function: copied from getRunningSpeedTuningCurves ('tuningCurve' framework) ----------
function tc = computeSpeedTuningCurve(sigData, wheelSpeed, fname, sourceStructName)
    fs = 60;
    mIdx = find(wheelSpeed > 1 & ~isnan(wheelSpeed));
    sIdx = find(wheelSpeed <= 1 & ~isnan(wheelSpeed));
    numROIs = size(sigData, 1);

    ptsPerBin = floor(0.07 * length(wheelSpeed));
    nB = floor(length(mIdx) / ptsPerBin);
    edges = quantile(wheelSpeed(mIdx), linspace(0, 1, max(1, nB) + 1));
    edges = unique(edges);
    numBins = length(edges) - 1;

    tc.source = sourceStructName;
    tc.speedBins = edges;
    tc.occupancy.stationary = length(sIdx) / fs;

    statMean = mean(sigData(:, sIdx), 2, 'omitnan');
    statSEM  = std(sigData(:, sIdx), 0, 2, 'omitnan') / sqrt(length(sIdx));
    mMean = zeros(numROIs, numBins);  mSEM = zeros(numROIs, numBins);  mCount = zeros(numROIs, numBins);
    for b = 1:numBins
        idx = find(wheelSpeed >= edges(b) & wheelSpeed < edges(b+1));
        mMean(:, b)  = mean(sigData(:, idx), 2, 'omitnan');
        mSEM(:, b)   = std(sigData(:, idx), 0, 2, 'omitnan') / sqrt(length(idx));
        mCount(:, b) = length(idx);
        tc.occupancy.moving(b) = length(idx) / fs;
    end
    realVarMovingOnly = var(mMean, 0, 2, 'omitnan');

    % Moving-only shuffle, as in getRunningSpeedTuningCurves
    numShifts = 1000;  minShift = 600;
    sigDataMovingOnly    = sigData(:, mIdx);
    wheelSpeedMovingOnly = wheelSpeed(mIdx);
    maxShiftMove = length(mIdx) - minShift;
    rng(1); shiftValuesMove = randi([minShift, maxShiftMove], [1, numShifts]);
    binMaskMove = false(numel(wheelSpeedMovingOnly), numBins);
    for b = 1:numBins
        binMaskMove(:, b) = wheelSpeedMovingOnly >= edges(b) & wheelSpeedMovingOnly < edges(b+1);
    end
    shuffVarsMoving = zeros(numROIs, numShifts);
    for iperm = 1:numShifts
        sigS_move = circshift(sigDataMovingOnly, shiftValuesMove(iperm), 2);
        mM = zeros(numROIs, numBins);
        for b = 1:numBins
            mM(:, b) = mean(sigS_move(:, binMaskMove(:, b)), 2, 'omitnan');
        end
        shuffVarsMoving(:, iperm) = var(mM, 0, 2, 'omitnan');
    end
    pValMoving = sum(shuffVarsMoving >= realVarMovingOnly, 2) / numShifts;
    pValMoving(realVarMovingOnly == 0) = 1.0;

    tc.(fname).statMean = statMean;   tc.(fname).statSEM  = statSEM;
    tc.(fname).statCount = repmat(length(sIdx), numROIs, 1);
    tc.(fname).moveMean = mMean;      tc.(fname).moveSEM  = mSEM;   tc.(fname).moveCount = mCount;
    tc.(fname).tuningVarMoving = realVarMovingOnly;
    tc.(fname).pValMoving = pValMoving;
end