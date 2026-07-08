% RFMapping_PlotGaussianFitSummary.m
%
% Summary figures built from the 2D Gaussian fit estimates (gaussFit_Az0,
% gaussFit_El0, gaussFit_sigmaAz, gaussFit_sigmaEl, gaussFit_isTrusted)
% written onto allRFMapping by RFMapping_Gaussian2DFit.m. Run that script
% first and keep allRFMapping in the workspace.
%
% Only TRUSTED fits are plotted (isStable AND ~isBeyondRange AND
% ~isDegenerateSigma). Boundary-adjacent-but-in-range fits are no longer
% distinguished by color -- all trusted fits are plotted uniformly.
% (isBeyondRange exclusion is unchanged: centers genuinely outside the
% tested/screen-capped visual field are still excluded from the trusted
% set entirely -- only the INTERIOR vs BOUNDARY-ADJACENT display
% distinction among the fits that ARE trusted has been removed.)
%
% Panels:
%   A. Preferred azimuth histogram (continuous, from Gaussian fit)
%   B. Preferred elevation histogram
%   C. 2D scatter of preferred position (azimuth vs elevation)
%   D. RF width asymmetry: histogram of paired difference
%      (sigmaAz - sigmaEl) per bouton, since the actual finding is a
%      paired, within-bouton difference, not a scatter-style correlation

%% gather trusted fits
hasFit = arrayfun(@(b) isfield(b, 'gaussFit_isTrusted') && ~isempty(b.gaussFit_isTrusted), allRFMapping);
if ~any(hasFit)
    error('No gaussFit_* fields found on allRFMapping -- run RFMapping_Gaussian2DFit.m first.');
end

isTrustedVec = false(numel(allRFMapping), 1);
az0Vec       = nan(numel(allRFMapping), 1);
el0Vec       = nan(numel(allRFMapping), 1);
sigAzVec     = nan(numel(allRFMapping), 1);
sigElVec     = nan(numel(allRFMapping), 1);

for i = find(hasFit)'   
    isTrustedVec(i) = allRFMapping(i).gaussFit_isTrusted;
    az0Vec(i)       = allRFMapping(i).gaussFit_Az0;
    el0Vec(i)       = allRFMapping(i).gaussFit_El0;
    sigAzVec(i)     = allRFMapping(i).gaussFit_sigmaAz;
    sigElVec(i)     = allRFMapping(i).gaussFit_sigmaEl;
end

plotIdx = find(isTrustedVec);

fprintf('Plotting %d trusted boutons.\n', numel(plotIdx));

if isempty(plotIdx)
    fitIdx = find(hasFit);
    nFit = numel(fitIdx);
    nStableOnly  = sum(arrayfun(@(i) isequal(allRFMapping(i).gaussFit_isStable, true), fitIdx));
    nBeyondRange = sum(arrayfun(@(i) isequal(allRFMapping(i).gaussFit_isBeyondRange, true), fitIdx));
    nDegenerate  = sum(arrayfun(@(i) isequal(allRFMapping(i).gaussFit_isDegenerateSigma, true), fitIdx));
    fprintf('\n*** No boutons passed all trusted criteria -- diagnostic breakdown (n = %d fitted boutons) ***\n', nFit);
    fprintf('  isStable = true for:            %d / %d\n', nStableOnly, nFit);
    fprintf('  isBeyondRange = true for:       %d / %d (these are EXCLUDED by design)\n', nBeyondRange, nFit);
    fprintf('  isDegenerateSigma = true for:   %d / %d (these are EXCLUDED by design)\n', nDegenerate, nFit);
    error('Zero boutons passed all trusted criteria -- check the breakdown above.');
end

%% compute paired sigma difference and stats
sigDiffVec = sigAzVec - sigElVec; % positive = azimuth wider than elevation, for THIS bouton
[pWilcoxSig, ~] = signrank(sigAzVec(plotIdx), sigElVec(plotIdx));
medianDiffSig = median(sigDiffVec(plotIdx), 'omitnan');

fitColor = [0.20 0.45 0.75];

%% ===================== figure =====================
figSummary = figure('Color', 'w', 'Position', [50 50 1200 950], 'Name', 'Gaussian fit summary');

%% Panel A: preferred azimuth histogram
subplot(2,2,1);
histogram(az0Vec(plotIdx), 12, 'FaceColor', fitColor, 'EdgeColor', 'k');
xlabel('Preferred Azimuth (\circ), Gaussian fit');
ylabel('Number of boutons');
title('Preferred azimuth (fitted x0)');
set(gca, 'Box', 'off', 'TickDir', 'out');

%% Panel B: preferred elevation histogram
subplot(2,2,2);
histogram(el0Vec(plotIdx), 10, 'FaceColor', fitColor, 'EdgeColor', 'k');
xlabel('Preferred Elevation (\circ), Gaussian fit');
ylabel('Number of boutons');
title('Preferred elevation (fitted y0)');
set(gca, 'Box', 'off', 'TickDir', 'out');

%% Panel C: 2D scatter of preferred position
subplot(2,2,3);
scatter(az0Vec(plotIdx), el0Vec(plotIdx), 30, fitColor, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Preferred Azimuth (\circ)');
ylabel('Preferred Elevation (\circ)');
title(sprintf('Preferred position (n = %d trusted boutons)', numel(plotIdx)));
axis equal;
set(gca, 'Box', 'off', 'TickDir', 'out');

%% Panel D: RF width asymmetry -- paired difference (sigmaAz - sigmaEl), not a scatter
subplot(2,2,4);
histogram(sigDiffVec(plotIdx), 12, 'FaceColor', fitColor, 'EdgeColor', 'k');
hold on;
xline(0, 'k--', 'LineWidth', 1);
xlabel('\sigma_{azimuth} - \sigma_{elevation} (\circ)');
ylabel('Number of boutons');
title(sprintf('RF width asymmetry (median = %.1f\\circ, p = %.4f)', medianDiffSig, pWilcoxSig), ...
    'FontWeight', 'normal');
set(gca, 'Box', 'off', 'TickDir', 'out');

sgtitle(sprintf('Gaussian RF fit summary (trusted fits only, n = %d)', numel(plotIdx)));

%% ===================== save: combined overview figure =====================
outputDirBase = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\gaussianfitsummary';

outputDir = fullfile(outputDirBase, 'overview');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figSummary, 'Visible', 'off');
saveFigureFormats(figSummary, fullfile(outputDir, 'gaussianfitsummary'));

%% ===================== save: preferred azimuth (standalone) =====================
figAz = figure('Color', 'w', 'Position', [50 50 500 400], 'Name', 'Preferred azimuth');
histogram(az0Vec(plotIdx), 12, 'FaceColor', fitColor, 'EdgeColor', 'k');
xlabel('Preferred Azimuth (\circ), Gaussian fit');
ylabel('Number of boutons');
title('Preferred azimuth (fitted x0)');
set(gca, 'Box', 'off', 'TickDir', 'out');

outputDir = fullfile(outputDirBase, 'preferred_azimuth');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figAz, 'Visible', 'off');
saveFigureFormats(figAz, fullfile(outputDir, 'gaussianfitsummary'));

%% ===================== save: preferred elevation (standalone) =====================
figEl = figure('Color', 'w', 'Position', [50 50 500 400], 'Name', 'Preferred elevation');
histogram(el0Vec(plotIdx), 10, 'FaceColor', fitColor, 'EdgeColor', 'k');
xlabel('Preferred Elevation (\circ), Gaussian fit');
ylabel('Number of boutons');
title('Preferred elevation (fitted y0)');
set(gca, 'Box', 'off', 'TickDir', 'out');

outputDir = fullfile(outputDirBase, 'preferred_elevation');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figEl, 'Visible', 'off');
saveFigureFormats(figEl, fullfile(outputDir, 'gaussianfitsummary'));

%% ===================== save: preferred position scatter (standalone) =====================
figPos = figure('Color', 'w', 'Position', [50 50 500 450], 'Name', 'Preferred position scatter');
scatter(az0Vec(plotIdx), el0Vec(plotIdx), 30, fitColor, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Preferred Azimuth (\circ)');
ylabel('Preferred Elevation (\circ)');
title(sprintf('Preferred position (n = %d trusted boutons)', numel(plotIdx)));
axis equal;
set(gca, 'Box', 'off', 'TickDir', 'out');

outputDir = fullfile(outputDirBase, 'preferred_position');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figPos, 'Visible', 'off');
saveFigureFormats(figPos, fullfile(outputDir, 'gaussianfitsummary'));

%% ===================== save: RF width asymmetry (standalone) =====================
figWidth = figure('Color', 'w', 'Position', [50 50 480 400], 'Name', 'RF width asymmetry sigmaAz-sigmaEl');
histogram(sigDiffVec(plotIdx), 12, 'FaceColor', fitColor, 'EdgeColor', 'k');
hold on;
xline(0, 'k--', 'LineWidth', 1);
xlabel('\sigma_{azimuth} - \sigma_{elevation} (\circ)');
ylabel('Number of boutons');
title(sprintf('RF width asymmetry (median = %.1f\\circ, p = %.4f)', medianDiffSig, pWilcoxSig), ...
    'FontWeight', 'normal');
set(gca, 'Box', 'off', 'TickDir', 'out');

outputDir = fullfile(outputDirBase, 'rf_width');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figWidth, 'Visible', 'off');
saveFigureFormats(figWidth, fullfile(outputDir, 'gaussianfitsummary'));

%% console summary stats
fprintf('\n=== Summary statistics (trusted fits, n = %d) ===\n', numel(plotIdx));
fprintf('Median preferred azimuth:   %.1f deg\n', median(az0Vec(plotIdx), 'omitnan'));
fprintf('Median preferred elevation: %.1f deg\n', median(el0Vec(plotIdx), 'omitnan'));
fprintf('Median sigma azimuth:       %.1f deg\n', median(sigAzVec(plotIdx), 'omitnan'));
fprintf('Median sigma elevation:     %.1f deg\n', median(sigElVec(plotIdx), 'omitnan'));
fprintf('Median sigma difference (az - el): %.1f deg\n', medianDiffSig);
fprintf('Wilcoxon signed-rank, sigmaAz vs sigmaEl: p = %.4f\n', pWilcoxSig);
fprintf(['NOTE: if p is significant, RF width differs systematically between axes on average --\n' ...
    'worth reporting which axis tends to be broader, rather than assuming RFs are round.\n']);
