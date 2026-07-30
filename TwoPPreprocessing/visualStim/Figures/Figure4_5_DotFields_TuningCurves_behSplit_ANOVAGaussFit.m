% visual speed tuning - behaviour split
% ALTERNATE VERSION: uses ANOVA(speeds+blank)+2SD responsiveness
% (isResponsive_stat/isResponsive_run) instead of cross-validated R^2 to
% decide "tuned", and the non-cross-validated Gaussian fit goodness
% (gaussR2_stat/gaussR2_run) to decide which fits are trustworthy enough
% to report parameters/preferred-speed/shape from.
%
% KEY DIFFERENCE FROM THE R^2-BASED VERSION: neither of these criteria
% requires matched trial counts between states, so NO downsampling is
% needed anywhere in this script -- every bouton uses its full available
% trial count in each state, independently.
%
% WHAT THIS MEANS THE R^2 PANEL CANNOT DO ANYMORE: the paired
% cross-validated R^2 comparison (stat vs run reliability, tested against
% a shuffled null) doesn't exist in this pipeline -- gaussR2 is a
% same-data goodness-of-fit measure, not a held-out reliability measure,
% so it will tend to look optimistic (especially for low-trial boutons)
% and should NOT be described as "tuning reliability". It answers a
% different, narrower question: "is a Gaussian a reasonable description
% of this bouton's tuning curve", not "would this same curve reproduce on
% new trials". Keep this distinction explicit in any write-up.

%% run this to load the appropriate data into workspace
DotFields_TuningCurveAnalysis_compareStatesV2_2PData

%% Post-hoc session exclusion (FAST -- does NOT rerun ANOVA/Gaussian fits)
% Same mechanism as the R^2-based script, in case you still want to drop
% a specific session for reasons other than trial count (e.g. suspected
% data-quality issues) -- this pipeline doesn't NEED it for fairness the
% way the R^2 version did, since nothing here requires matched trials.
sessionLabels_all = {allDotUnits.sessionLabel}';
uniqueSessionsAll  = unique(sessionLabels_all, 'stable');

fprintf('\n--- Per-session bouton counts (for reference only -- not a fairness requirement here) ---\n');
for s = 1:numel(uniqueSessionsAll)
    thisSession = uniqueSessionsAll{s};
    fprintf('%-25s %6d boutons\n', thisSession, sum(strcmp(sessionLabels_all, thisSession)));
end

% EDIT THIS LIST to exclude sessions from all panels below.
sessionsToExclude = {};

excludeBoutonMask = ismember(sessionLabels_all, sessionsToExclude);
if any(excludeBoutonMask)
    fprintf('\nExcluding %d boutons from %d session(s): %s\n', ...
        sum(excludeBoutonMask), numel(sessionsToExclude), strjoin(sessionsToExclude, ', '));
else
    fprintf('\nNo sessions excluded -- using all %d boutons.\n', numel(allDotUnits));
end

%% Build the ANOVA+2SD / gaussR2 combined criterion
% "Tuned" for prevalence purposes = ANOVA(speeds+blank)+2SD-above-blank
% responsive (isResponsive_stat/_run), computed from FULL trial data, no
% downsampling. This is a real hypothesis test on its own.
%
% "Trustworthy for shape/preferred-speed/example-plot purposes" ADDS a
% goodness-of-fit gate on top: responsive AND gaussR2 above a threshold
% you should set by actually looking at the gaussR2 distribution first
% (not borrowed from the R^2>0.1 convention -- that was for a different,
% cross-validated statistic).

gaussR2Thresh = 0.3; % PLACEHOLDER -- inspect the histogram below before trusting this

isResponsive_stat_all = cat(1, allDotUnits.isResponsive_stat);
isResponsive_run_all  = cat(1, allDotUnits.isResponsive_run);
gaussR2_stat_all      = cat(1, allDotUnits.gaussR2_stat);
gaussR2_run_all       = cat(1, allDotUnits.gaussR2_run);

isResponsive_stat_all(excludeBoutonMask) = false;
isResponsive_run_all(excludeBoutonMask)  = false;
gaussR2_stat_all(excludeBoutonMask)      = NaN;
gaussR2_run_all(excludeBoutonMask)       = NaN;

figure('Color', 'w', 'Name', 'Gaussian fit R^2 distribution (non-cross-validated)', 'Position', [200, 200, 500, 400]);
histogram(gaussR2_stat_all(~isnan(gaussR2_stat_all)), 20, 'FaceColor', 'k', 'FaceAlpha', 0.5, 'DisplayName', 'Stationary');
hold on;
histogram(gaussR2_run_all(~isnan(gaussR2_run_all)), 20, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'Locomotion');
xline(gaussR2Thresh, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'DisplayName', sprintf('threshold=%.2f', gaussR2Thresh));
xlabel('Gaussian fit R^2 (non-cross-validated, same-data goodness of fit)');
ylabel('# boutons');
title('Check this BEFORE trusting gaussR2Thresh above');
legend('Location', 'best'); box off; hold off;

% "Tuned" (prevalence question) -- responsiveness only, no fit-quality gate
isTunedStat_alt = isResponsive_stat_all;
isTunedRun_alt  = isResponsive_run_all;

% "Trustworthy for shape/parameters" -- responsive AND well-fit
validIdx_stat_alt = find(isResponsive_stat_all & gaussR2_stat_all > gaussR2Thresh);
validIdx_run_alt  = find(isResponsive_run_all  & gaussR2_run_all  > gaussR2Thresh);
validIdx_both_alt = validIdx_stat_alt(ismember(validIdx_stat_alt, validIdx_run_alt));

fprintf('\n--- ANOVA+2SD / gaussR2 combined criterion ---\n');
fprintf('Responsive (stationary): %d / %d\n', sum(isTunedStat_alt), numel(allDotUnits));
fprintf('Responsive (locomotion): %d / %d\n', sum(isTunedRun_alt),  numel(allDotUnits));
fprintf('Responsive + well-fit (stationary): %d\n', numel(validIdx_stat_alt));
fprintf('Responsive + well-fit (locomotion): %d\n', numel(validIdx_run_alt));
fprintf('Responsive + well-fit in BOTH states: %d\n', numel(validIdx_both_alt));

%% PANEL: Goodness-of-fit paired comparison (signrank + scatter)
% NOT a reliability comparison -- see header note. Restricted to boutons
% responsive in both states (fair comparison of "how Gaussian-shaped is
% the curve" only among boutons that are actually driven by the stimulus
% in both states).
pairIdx_alt = isResponsive_stat_all & isResponsive_run_all & ~isnan(gaussR2_stat_all) & ~isnan(gaussR2_run_all);
gaussR2_stat_paired = gaussR2_stat_all(pairIdx_alt);
gaussR2_run_paired  = gaussR2_run_all(pairIdx_alt);

pval_gaussR2 = signrank(gaussR2_stat_paired, gaussR2_run_paired);

figure('Color', 'w', 'Name', 'Goodness-of-fit paired comparison', 'Position', [200, 200, 450, 450]);
scatter(gaussR2_stat_paired, gaussR2_run_paired, 20, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.4);
hold on;
minVal = min([gaussR2_stat_paired; gaussR2_run_paired]); maxVal = max([gaussR2_stat_paired; gaussR2_run_paired]);
plot([minVal maxVal], [minVal maxVal], 'r--', 'LineWidth', 1.2);
xlabel('Stationary gaussR^2 (non-cross-validated)'); ylabel('Locomotion gaussR^2 (non-cross-validated)');
title(sprintf('Gaussian fit quality, paired per bouton (n=%d)\nWilcoxon signed-rank p=%.3g\n(NOT a reliability comparison -- see header note)', ...
    numel(gaussR2_stat_paired), pval_gaussR2), 'FontSize', 10);
axis square; box off;
hold off;

%% PANEL: Prevalence of responsiveness, paired (McNemar test)
tunedStat_p = isTunedStat_alt;
tunedRun_p  = isTunedRun_alt;

nBothTuned    = sum(tunedStat_p & tunedRun_p);
nStatOnly     = sum(tunedStat_p & ~tunedRun_p);
nRunOnly      = sum(~tunedStat_p & tunedRun_p);
nNeitherTuned = sum(~tunedStat_p & ~tunedRun_p);

fprintf('\n--- Prevalence (ANOVA+2SD responsive), paired (n=%d boutons) ---\n', numel(tunedStat_p));
fprintf('%-20s %6d\n', 'Responsive in both', nBothTuned);
fprintf('%-20s %6d\n', 'Stationary only',    nStatOnly);
fprintf('%-20s %6d\n', 'Locomotion only',    nRunOnly);
fprintf('%-20s %6d\n', 'Neither',            nNeitherTuned);

% NOTE: chi-square continuity-corrected McNemar is only reliable when
% n12+n21 is reasonably large (rule of thumb >=10) -- check this before
% trusting the p-value.
n12 = nStatOnly; n21 = nRunOnly;
mcnemarChi2 = (abs(n12 - n21) - 1)^2 / (n12 + n21);
pval_mcnemar = 1 - chi2cdf(mcnemarChi2, 1);
fprintf('McNemar test (stat-only vs run-only): chi2=%.2f, p=%.3g (n12+n21=%d)\n', ...
    mcnemarChi2, pval_mcnemar, n12+n21);

figure('Color', 'w', 'Name', 'Prevalence, paired (ANOVA+2SD)', 'Position', [200, 200, 400, 400]);
bar([nBothTuned, nStatOnly, nRunOnly, nNeitherTuned], 'FaceColor', [0.3 0.5 0.7]);
set(gca, 'XTickLabel', {'Both', 'Stat only', 'Run only', 'Neither'});
ylabel('# boutons');
title(sprintf('Responsiveness prevalence, paired (n=%d)\nMcNemar p=%.3g', numel(tunedStat_p), pval_mcnemar), ...
    'FontSize', 11);
box off;

%% PANEL: Preferred speed agreement, paired (all shape categories)
if ~isempty(validIdx_both_alt)
    prefSpeed_stat_both = cat(1, allDotUnits(validIdx_both_alt).prefSpeed_stat);
    prefSpeed_run_both  = cat(1, allDotUnits(validIdx_both_alt).prefSpeed_run);

    pval_prefSpeed = signrank(prefSpeed_stat_both, prefSpeed_run_both);
    pctAgree = 100 * mean(prefSpeed_stat_both == prefSpeed_run_both);

    figure('Color', 'w', 'Name', 'Preferred speed agreement (ANOVA+gaussR2)', 'Position', [200, 200, 450, 400]);
    vals = histcounts2(prefSpeed_stat_both, prefSpeed_run_both, 0.5:1:6.5, 0.5:1:6.5);
    imagesc(vals'); axis xy; colorbar;
    hold on; plot([0.5 6.5], [0.5 6.5], 'r');
    xlabel('Preferred speed index (stationary)'); ylabel('Preferred speed index (locomotion)');
    title(sprintf('Preferred speed agreement (n=%d)\n%.1f%% exact match, signrank p=%.3g', ...
        numel(validIdx_both_alt), pctAgree, pval_prefSpeed), 'FontSize', 11);
else
    fprintf('\nvalidIdx_both_alt is empty -- no boutons responsive+well-fit in both states. Consider lowering gaussR2Thresh.\n');
end

%% PANEL: Tuning shape stability, paired (confusion matrix, category 1-4)
if ~isempty(validIdx_both_alt)
    gaussChar_stat_both = cat(1, allDotUnits(validIdx_both_alt).gaussChar_stat);
    gaussChar_run_both  = cat(1, allDotUnits(validIdx_both_alt).gaussChar_run);

    validShapeIdx = ~isnan(gaussChar_stat_both) & ~isnan(gaussChar_run_both);
    gcs = gaussChar_stat_both(validShapeIdx);
    gcr = gaussChar_run_both(validShapeIdx);

    shapeConfusion = histcounts2(gcs, gcr, 0.5:1:4.5, 0.5:1:4.5);
    pctShapeAgree = 100 * sum(diag(shapeConfusion)) / sum(shapeConfusion(:));

    categoryNames = {'Lowpass', 'Highpass', 'Bandpass', 'Trough'};
    figure('Color', 'w', 'Name', 'Tuning shape stability (ANOVA+gaussR2)', 'Position', [200, 200, 450, 400]);
    imagesc(shapeConfusion');
    axis xy; colorbar;
    set(gca, 'XTick', 1:4, 'XTickLabel', categoryNames, 'YTick', 1:4, 'YTickLabel', categoryNames);
    xlabel('Tuning shape (stationary)'); ylabel('Tuning shape (locomotion)');
    title(sprintf('Tuning shape stability (n=%d)\n%.1f%% same shape in both states', ...
        numel(gcs), pctShapeAgree), 'FontSize', 11);
    for r = 1:4
        for c = 1:4
            text(r, c, num2str(shapeConfusion(r,c)), 'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
        end
    end
end

%% PANEL: Tuning category proportions, stationary vs locomotion (grouped bar)
categoryNames = {'Lowpass', 'Highpass', 'Bandpass', 'Trough'};

gaussChar_stat_valid = cat(1, allDotUnits(validIdx_stat_alt).gaussChar_stat);
gaussChar_run_valid  = cat(1, allDotUnits(validIdx_run_alt).gaussChar_run);

catCountsStat = histcounts(gaussChar_stat_valid, 0.5:1:4.5);
catCountsRun  = histcounts(gaussChar_run_valid, 0.5:1:4.5);
catPropsStat  = catCountsStat / sum(catCountsStat);
catPropsRun   = catCountsRun  / sum(catCountsRun);

figure('Color', 'w', 'Name', 'Tuning category proportions (ANOVA+gaussR2)', 'Position', [200, 200, 500, 400]);
b = bar([catPropsStat; catPropsRun]', 'grouped');
b(1).FaceColor = [0.15 0.15 0.15];
b(2).FaceColor = [0.80 0.20 0.20];
set(gca, 'XTick', 1:4, 'XTickLabel', categoryNames);
ylabel('Proportion of responsive+well-fit boutons');
legend({sprintf('Stationary (n=%d)', numel(validIdx_stat_alt)), sprintf('Locomotion (n=%d)', numel(validIdx_run_alt))}, ...
    'Location', 'best');
title('Tuning shape category, by behavioural state');
box off;

%% Select representative dual-tuned boutons per category, for example overlay plots
% Ranked by the WORSE of the two states' gaussR2 (goodness-of-fit, not
% cross-validated reliability), same-shape-in-both-states requirement,
% diversity check across sessions -- same pattern as the R^2 version.
nPerCategory = 3;

exampleBoutons = cell(4,1);
if ~isempty(validIdx_both_alt)
    gaussChar_stat_both = cat(1, allDotUnits(validIdx_both_alt).gaussChar_stat);
    gaussChar_run_both  = cat(1, allDotUnits(validIdx_both_alt).gaussChar_run);
    gaussR2_stat_both = cat(1, allDotUnits(validIdx_both_alt).gaussR2_stat);
    gaussR2_run_both  = cat(1, allDotUnits(validIdx_both_alt).gaussR2_run);
    sessionKey_both = {allDotUnits(validIdx_both_alt).sessionLabel}';
    combinedGaussR2 = min(gaussR2_stat_both, gaussR2_run_both);

    for c = 1:4
        catIdx = find(gaussChar_stat_both == c & gaussChar_run_both == c);
        [~, order] = sort(combinedGaussR2(catIdx), 'descend');
        catIdx = catIdx(order);

        used = {}; picked = [];
        for k = 1:numel(catIdx)
            key = sessionKey_both{catIdx(k)};
            if ismember(key, used), continue; end
            picked(end+1) = validIdx_both_alt(catIdx(k));
            used{end+1} = key;
            if numel(picked) >= nPerCategory, break; end
        end
        exampleBoutons{c} = picked;
    end

    fprintf('\n--- Example boutons selected per category (responsive+well-fit, same shape both states) ---\n');
    for c = 1:4
        fprintf('%-10s: %d examples found\n', categoryNames{c}, numel(exampleBoutons{c}));
    end
end

%% PANEL: Superimposed stationary vs locomotion tuning curves, representative examples
speeds = [0 16 32 64 128 256];
xDense = linspace(1, 6, 100);

figure('Color', 'w', 'Name', 'Example tuning curves, stat vs run (ANOVA+gaussR2)', 'Position', [100, 50, 1200, 900]);
for c = 1:4
    picked = exampleBoutons{c};
    for k = 1:numel(picked)
        subplot(4, nPerCategory, (c-1)*nPerCategory + k); hold on;

        thisBouton = picked(k);
        y_data_stat = allDotUnits(thisBouton).tuning(:,1);
        y_data_run  = allDotUnits(thisBouton).tuning(:,2);
        p_stat = allDotUnits(thisBouton).gaussParams_stat;
        p_run  = allDotUnits(thisBouton).gaussParams_run;
        y_fit_stat = p_stat(1) + p_stat(2)*exp(-((xDense-p_stat(3)).^2)/(2*p_stat(4)^2));
        y_fit_run  = p_run(1)  + p_run(2)*exp(-((xDense-p_run(3)).^2)/(2*p_run(4)^2));

        plot(1:6, y_data_stat, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
        plot(xDense, y_fit_stat, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Stationary');
        plot(1:6, y_data_run, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
        plot(xDense, y_fit_run, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Locomotion');

        xlim([0.5 6.5]); xticks(1:6);
        xticklabels(arrayfun(@(v) sprintf('%g', v), speeds, 'UniformOutput', false));
        ylabel('\DeltaF/F'); xlabel('Visual speed (\circ/s)');
        title(sprintf('%s (bouton %d)\ngaussR^2: stat=%.2f, run=%.2f', categoryNames{c}, thisBouton, ...
            allDotUnits(thisBouton).gaussR2_stat, allDotUnits(thisBouton).gaussR2_run), 'FontSize', 9);
        if c == 1 && k == 1, legend('Location', 'best'); end
        box off;
    end
end
sgtitle('Representative responsive+well-fit boutons: Stationary (black) vs Locomotion (red)', ...
    'FontSize', 13, 'FontWeight', 'bold');
