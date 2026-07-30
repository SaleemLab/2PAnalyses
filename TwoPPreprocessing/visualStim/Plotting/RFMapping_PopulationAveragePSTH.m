% RFMapping_PopulationAveragePSTH.m
%
% Population-averaged PSTH for RF mapping, analogous to
% DotFields_PopulationAveragePSTH.m, to determine a data-driven response
% window (onset latency, peak, half-decay) instead of assuming one.
%
% RUN THIS AFTER analyse_RFBoutons_basic.m (up to and including the
% "Recompute responsiveness and RF centers" section) -- it reuses:
%   allRFMapping          (pooled bouton struct array)
%   timeVector            (shared time vector, from RFMappingMetadata)
%   uAz, uEl_plot          (grid azimuth/elevation values)
%   respIdxList            (indices of boutons classified as responsive)
%
% For each responsive bouton, pulls the full trial-averaged trace AT ITS
% OWN PREFERRED GRID POSITION (centerAz/centerEl, already computed in the
% main script), then pools across boutons to get population mean +/- SEM.

%%
minSustainSec = 0.1;  % onset must hold above threshold for this long
% NOTE: adjust upward (e.g. 0.2-0.3) if RF mapping's per-trial SNR is
% lower than dot fields' (fewer trials per grid position is common for
% RF mapping, so a shorter/noisier sustain window is more likely to
% trigger a false-positive onset -- sanity check against the plotted
% trace either way, same as we did for dot fields).

%%
if isempty(allRFMapping) || isempty(timeVector)
    error('allRFMapping / timeVector not found -- run analyse_RFBoutons_basic.m first.');
end
if ~exist('respIdxList', 'var') || isempty(respIdxList)
    respIdxList = find([allRFMapping.isResponsive]);
end

fprintf('Building population PSTH from %d RF-responsive boutons...\n', numel(respIdxList));

psthTraces = {};
skippedCount = 0;

for k = 1:numel(respIdxList)
    iROI = respIdxList(k);
    b = allRFMapping(iROI);

    centerAz = b.centerAz;
    centerEl = b.centerEl;
    if isnan(centerAz) || isnan(centerEl)
        skippedCount = skippedCount + 1;
        continue;
    end

    cIdx = find(uAz == centerAz, 1);
    rIdx = find(uEl_plot == centerEl, 1);
    if isempty(cIdx) || isempty(rIdx)
        skippedCount = skippedCount + 1;
        continue;
    end

    trialMatrix = b.baselineSubtracted{rIdx, cIdx}; % [Trials x Time]
    if isempty(trialMatrix)
        skippedCount = skippedCount + 1;
        continue;
    end

    meanTrace = mean(double(trialMatrix), 1, 'omitnan');
    if all(isnan(meanTrace))
        skippedCount = skippedCount + 1;
        continue;
    end
    psthTraces{end+1} = meanTrace; %#ok<SAGROW>
end

fprintf('Built %d population PSTH traces (%d skipped -- missing position/data).\n', ...
    numel(psthTraces), skippedCount);

matAll = cat(1, psthTraces{:});
meanAll = mean(matAll, 1, 'omitnan');
semAll  = std(matAll, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(matAll), 1));

%% plot population average +/- SEM
figRF = figure('Color','w','Position',[200 200 800 450]);
hold on;
fill([timeVector(:)', fliplr(timeVector(:)')], [meanAll+semAll, fliplr(meanAll-semAll)], ...
    [0.2 0.2 0.6], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');
plot(timeVector, meanAll, '-', 'Color', [0.2 0.2 0.6], 'LineWidth', 2, ...
    'DisplayName', sprintf('RF-responsive, preferred position (n=%d)', size(matAll,1)));
xline(0, 'k--', 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('\DeltaF/F (baseline-subtracted, mean across responsive boutons)');
title('Population-averaged PSTH (RF mapping, preferred position)');
legend('Location','best');
box off;

%% quantify peak + half-decay
[peakVal, peakIdx] = max(meanAll);
fprintf('\nPeak: %.4f at t=%.2fs\n', peakVal, timeVector(peakIdx));

baselineVal = mean(meanAll(timeVector < 0), 'omitnan');
halfVal = baselineVal + 0.5*(peakVal - baselineVal);
postPeakIdx = peakIdx:numel(timeVector);
decayIdx = postPeakIdx(find(meanAll(postPeakIdx) < halfVal, 1));
if ~isempty(decayIdx)
    fprintf('Decays to half-peak at t=%.2fs\n', timeVector(decayIdx));
else
    fprintf('Does not decay to half-peak within the trial window -- inspect trace manually.\n');
end

%% data-driven onset latency (first sustained crossing of baseline + 2SD)
dt = mean(diff(timeVector));
minSustainSamples = max(1, round(minSustainSec / dt));

baselineSD = std(meanAll(timeVector < 0), 'omitnan');
thresh = baselineVal + 2*baselineSD;

postStimIdx = find(timeVector >= 0);
above = meanAll(postStimIdx) > thresh;
onsetIdx = [];
for kk = 1:(numel(above) - minSustainSamples + 1)
    if all(above(kk:kk+minSustainSamples-1))
        onsetIdx = postStimIdx(kk);
        break;
    end
end

if ~isempty(onsetIdx)
    fprintf('Onset (first sustained crossing of baseline+2SD, held %.0fms): t=%.2fs (baseline=%.4f, threshold=%.4f)\n', ...
        minSustainSec*1000, timeVector(onsetIdx), baselineVal, thresh);
    xline(timeVector(onsetIdx), 'k:', 'Onset', 'HandleVisibility','off');
else
    fprintf('No sustained onset crossing found -- inspect trace manually / consider lowering threshold or minSustainSec.\n');
end

if ~isempty(onsetIdx) && ~isempty(decayIdx)
    fprintf('\nSuggested response window: [%.2f %.2f]\n', timeVector(onsetIdx), timeVector(decayIdx));
    fprintf('(Current script default respWin = [0.5 2] -- compare against this before deciding whether to update it.)\n');
end
