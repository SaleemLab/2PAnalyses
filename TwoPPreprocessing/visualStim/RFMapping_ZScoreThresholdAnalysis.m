% RFMapping_ZScoreThresholdAnalysis.m
%
% Rather than debating 1SD vs 2SD in the abstract, computes a continuous
% z-score per bouton: z = (prefVal - blankMean) / blankStd -- "how many
% SDs above blank is this bouton's preferred response, exactly" -- and
% uses ANOVA significance (a DIFFERENT, independent test: do conditions
% differ at all) as a rough validation signal for where genuinely-tuned
% boutons' z-scores actually fall.
%
% Produces:
%   1. Histogram of z, split into ANOVA-significant vs not.
%   2. Threshold-sweep: fraction of ANOVA-significant boutons retained
%      as the SD cutoff varies continuously from 0 to 3.
%   3. Visual spot-check: pulls up plotRFHeatmapWithTraces for a sample
%      of "borderline" boutons (ANOVA-significant, but z between 1 and
%      2 -- i.e. boutons the 2SD rule would reject but 1SD would keep)
%      so you can look at them directly.
%
% Requires allRFMapping already has pValANOVA computed (from the
% recomputation loop), and uAz/uEl_plot/timeVector in the workspace for
% the visual spot-check section.

ALPHA = 0.05;
numBoutons = numel(allRFMapping);

zScore    = nan(numBoutons, 1);
prefValAll  = nan(numBoutons, 1);
blankMeanAll = nan(numBoutons, 1);
blankStdAll  = nan(numBoutons, 1);

for iROI = 1:numBoutons
    bTrialsCorrected = allRFMapping(iROI).baselineSubtractedBlank;
    meanGridResponse = allRFMapping(iROI).meanGridResponse;

    if isempty(bTrialsCorrected) || isempty(meanGridResponse)
        continue;
    end
    bTrialsCorrected = double(bTrialsCorrected);
    meanGridResponse = double(meanGridResponse);

    % NOTE: respIdx must already be defined in the workspace (same as
    % used in the main recomputation loop) for this to match exactly.
    blankTrialMeans = mean(bTrialsCorrected(:, respIdx), 2, 'omitnan');
    blankMean = mean(blankTrialMeans, 'omitnan');
    blankStd  = std(blankTrialMeans, 'omitnan');
    prefVal   = max(meanGridResponse(:), [], 'omitnan');

    if blankStd > 0
        zScore(iROI) = (prefVal - blankMean) / blankStd;
    end
    prefValAll(iROI)   = prefVal;
    blankMeanAll(iROI) = blankMean;
    blankStdAll(iROI)  = blankStd;
end

for iROI = 1:numBoutons
    allRFMapping(iROI).zScore    = zScore(iROI);
    allRFMapping(iROI).prefVal   = prefValAll(iROI);
    allRFMapping(iROI).blankMean = blankMeanAll(iROI);
    allRFMapping(iROI).blankStd  = blankStdAll(iROI);
end

anovaSig = [allRFMapping.pValANOVA] < ALPHA;
validZ   = ~isnan(zScore);

fprintf('%d / %d boutons have a valid z-score.\n', sum(validZ), numBoutons);
fprintf('%d boutons pass ANOVA (p<%.2f).\n', sum(anovaSig & validZ'), ALPHA);

%% ===================== Figure 1: z-score histogram, split by ANOVA significance =====================
figure('Position', [100 100 600 450]);
edges = -2:0.2:8;
histogram(zScore(validZ & ~anovaSig'), edges, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.6, 'DisplayName', 'ANOVA not significant');
hold on;
histogram(zScore(validZ & anovaSig'), edges, 'FaceColor', [0.2 0.5 0.8], 'FaceAlpha', 0.7, 'DisplayName', 'ANOVA significant (p<0.05)');
xline(1, 'k--', 'LineWidth', 1.5, 'Label', '1 SD');
xline(2, 'k-',  'LineWidth', 1.5, 'Label', '2 SD');
xlabel('z-score: (prefVal - blankMean) / blankStd'); ylabel('Number of boutons');
legend('Location', 'best');
title('Distribution of preferred-response z-scores');

nAnovaSig = sum(anovaSig' & validZ);
n_above1  = sum(anovaSig' & validZ & zScore > 1);
n_above2  = sum(anovaSig' & validZ & zScore > 2);
n_between = sum(anovaSig' & validZ & zScore > 1 & zScore <= 2);
fprintf('\nAmong ANOVA-significant boutons (n=%d):\n', nAnovaSig);
fprintf('  z > 1: %d (%.1f%%)\n', n_above1, 100*n_above1/nAnovaSig);
fprintf('  z > 2: %d (%.1f%%)\n', n_above2, 100*n_above2/nAnovaSig);
fprintf('  Between 1 and 2 (kept by 1SD, rejected by 2SD): %d (%.1f%%)\n', n_between, 100*n_between/nAnovaSig);

%%  Figure 2: threshold sweep 
sdRange = 0:0.1:3;
fracRetained = nan(size(sdRange));
for i = 1:numel(sdRange)
    fracRetained(i) = sum(anovaSig' & validZ & zScore > sdRange(i)) / nAnovaSig;
end

figure('Position', [100 100 500 400]);
plot(sdRange, 100*fracRetained, 'k-', 'LineWidth', 2);
hold on;
xline(1, 'b--', 'LineWidth', 1.2, 'Label', '1 SD');
xline(2, 'r--', 'LineWidth', 1.2, 'Label', '2 SD');
xlabel('SD threshold'); ylabel('% of ANOVA-significant boutons retained');
title('How much does raising the SD threshold cost you?');
grid on;

%% visual spot-check: borderline boutons (ANOVA-sig, 1<z<=2) 
borderlineIdx = find(anovaSig(:) & validZ(:) & zScore > 1 & zScore <= 2);
fprintf('\n%d borderline boutons (ANOVA-significant, 1 < z <= 2) available for visual spot-check.\n', numel(borderlineIdx));

nToShow = min(6, numel(borderlineIdx));
if nToShow > 0
    exampleBorderline = borderlineIdx(round(linspace(1, numel(borderlineIdx), nToShow)));
    figure('Color', 'w', 'Position', [50 50 nToShow*300 550], 'Name', 'Borderline boutons (1<z<=2, ANOVA-significant)');
    for i = 1:nToShow
        b = exampleBorderline(i);
        axPanel = subplot(2, nToShow, i);
        plotRFHeatmapWithTraces(axPanel, allRFMapping(b), uAz, uEl_plot, timeVector, 'Smooth', true);
        title(axPanel, sprintf('Bouton %d\nz=%.2f, p=%.3g', b, zScore(b), allRFMapping(b).pValANOVA), 'FontSize', 9);

        axBlank = subplot(2, nToShow, nToShow + i);
        localPlotBlankVsPref(axBlank, allRFMapping(b), respIdx);
    end
else
    fprintf('No borderline boutons found -- nothing to plot.\n');
end

%% visual spot-check #2: LOWER slice (1.0 <= z < 1.2) 
% This is specifically the slice a 1SD cutoff would newly include beyond
% what's already been visually checked at 1.2-1.4. If these still look
% like coherent RFs (not noise), that supports using 1SD as the final
% cutoff. If they start looking noisy, a higher cutoff (e.g. 1.1-1.2)
% may be more appropriate than a full drop to 1.0.
lowerSliceIdx = find(anovaSig(:) & validZ(:) & zScore >= 1.0 & zScore < 1.2);
fprintf('\n%d boutons in the LOWER slice (ANOVA-significant, 1.0 <= z < 1.2) available for visual spot-check.\n', numel(lowerSliceIdx));

nToShow2 = min(6, numel(lowerSliceIdx));
if nToShow2 > 0
    exampleLowerSlice = lowerSliceIdx(round(linspace(1, numel(lowerSliceIdx), nToShow2)));
    figure('Color', 'w', 'Position', [50 50 nToShow2*300 550], 'Name', 'Lower-slice boutons (1.0<=z<1.2, ANOVA-significant)');
    for i = 1:nToShow2
        b = exampleLowerSlice(i);
        axPanel = subplot(2, nToShow2, i);
        plotRFHeatmapWithTraces(axPanel, allRFMapping(b), uAz, uEl_plot, timeVector, 'Smooth', true);
        title(axPanel, sprintf('Bouton %d\nz=%.2f, p=%.3g', b, zScore(b), allRFMapping(b).pValANOVA), 'FontSize', 9);

        axBlank = subplot(2, nToShow2, nToShow2 + i);
        localPlotBlankVsPref(axBlank, allRFMapping(b), respIdx);
    end
else
    fprintf('No lower-slice boutons found -- nothing to plot.\n');
end

%% local helper: blank trials vs preferred response
function localPlotBlankVsPref(ax, boutonData, respIdx)
    % Shows individual blank trial values (jittered scatter), the blank
    % mean +/- 1SD/2SD lines, and the preferred-position response marked
    % against them -- makes the "how many SDs above blank" story visible,
    % not just a number.
    bTrials = double(boutonData.baselineSubtractedBlank); % [Trials x Time]

    blankTrialMeans = mean(bTrials(:, respIdx), 2, 'omitnan');
    blankMean = mean(blankTrialMeans, 'omitnan');
    blankStd  = std(blankTrialMeans, 'omitnan');
    prefVal   = boutonData.prefVal;

    axes(ax); hold(ax, 'on');
    jitterX = 1 + 0.15*(rand(size(blankTrialMeans))-0.5);
    scatter(ax, jitterX, blankTrialMeans, 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Blank trials');

    yline(ax, blankMean, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Blank mean');
    yline(ax, blankMean + blankStd,   'b--', 'LineWidth', 1, 'DisplayName', '+1 SD');
    yline(ax, blankMean + 2*blankStd, 'r--', 'LineWidth', 1, 'DisplayName', '+2 SD');

    scatter(ax, 1, prefVal, 80, [0.9 0.3 0.1], 'filled', 'Marker', 'p', 'DisplayName', 'Preferred response');

    xlim(ax, [0.5 1.5]);
    set(ax, 'XTick', []);
    ylabel(ax, '\DeltaF/F');
    box(ax, 'off');
    if prefVal > blankMean, sign_str = '+'; else, sign_str = ''; end
    title(ax, sprintf('z=%s%.2f SD', sign_str, (prefVal-blankMean)/max(blankStd,eps)), 'FontSize', 8);
end