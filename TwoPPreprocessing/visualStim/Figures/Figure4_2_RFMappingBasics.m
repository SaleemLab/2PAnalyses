
%% panel a
% exampleBoutons = [1968,3850,2408];   %1991,3189 % this is good 3137    
% narro 196/ 875 this is goof 
% figure to show basic rf mapping responses (boutons)
% panel A - example grid traces (3 responsive and 2 unresponsive);
% superimposed on the gaussian 2d fit
% panel B -
% panel C - preferred azimuth/elevation histograms
% panel D - within-FOV deviation (az/el separate, and joint 2D density)
% panel E - population RF coverage map

%% run to load all processed data; include the 2d gaussian fits
analyse_RFBoutons_current
%%

trustedMask = [gaussFitResults.isTrusted];
trustedGI   = find(trustedMask);                    % indices into the FULL gaussFitResults
trustedIROI = [gaussFitResults(trustedGI).iROI];     % actual bouton iROI numbers, trusted only

fprintf('Using %d / %d trusted boutons for all figure panels.\n', numel(trustedGI), numel(gaussFitResults));

% Filtered copy containing ONLY the trusted set, for panels that want to index sequentially
% (1:numel) rather than by boolean mask. ALWAYS index into gaussFitResultsTrusted (never the
% original gaussFitResults) when using an index derived from THIS variable.
gaussFitResultsTrusted = gaussFitResults(trustedGI);


%%
sigAzTrusted = [gaussFitResults(trustedMask).sigmaX];
sigElTrusted = [gaussFitResults(trustedMask).sigmaY];

prefAz = [gaussFitResults(trustedMask).x0];
prefEl = [gaussFitResults(trustedMask).y0];

fprintf('sigAz  = %.2f +/- %.2f (SD), n = %d\n', mean(sigAzTrusted), std(sigAzTrusted), numel(sigAzTrusted));
fprintf('sigEl  = %.2f +/- %.2f (SD), n = %d\n', mean(sigElTrusted), std(sigElTrusted), numel(sigElTrusted));
fprintf('prefAz = %.2f +/- %.2f (SD), n = %d\n', mean(prefAz), std(prefAz), numel(prefAz));
fprintf('prefEl = %.2f +/- %.2f (SD), n = %d\n', mean(prefEl), std(prefEl), numel(prefEl));


% sigAz  = 14.74 +/- 8.52 (SD), n = 62
% sigEl  = 8.40 +/- 3.63 (SD), n = 62
% prefAz = -27.53 +/- 26.04 (SD), n = 62
% prefEl = 11.33 +/- 12.63 (SD), n = 62


%%

% Among the CURRENT trusted set (n=82)
nBeyondTrusted = sum([gaussFitResults(trustedGI).isBeyondRange]);
fprintf('%d / %d trusted boutons have a center beyond the stimulus-extended range.\n', ...
    nBeyondTrusted, numel(trustedGI));

% Broken down by axis
nBeyondAzTrusted = sum([gaussFitResults(trustedGI).isBeyondRangeAz]);
nBeyondElTrusted = sum([gaussFitResults(trustedGI).isBeyondRangeEl]);
fprintf('  Azimuth beyond range: %d\n', nBeyondAzTrusted);
fprintf('  Elevation beyond range: %d\n', nBeyondElTrusted);

% For reference, among ALL above-floor fits (not just trusted)
nBeyondAll = sum([gaussFitResults.isBeyondRange]);
fprintf('%d / %d ALL above-floor fits have a center beyond the stimulus-extended range.\n', ...
    nBeyondAll, numel(gaussFitResults));

% Trusted boutons with Az center beyond the stimulus-extended range
azBeyondIdx = trustedGI([gaussFitResults(trustedGI).isBeyondRangeAz]);
fprintf('\n--- Trusted boutons with Az beyond range (n=%d) ---\n', numel(azBeyondIdx));
fprintf('%-8s %-10s %-10s %-8s\n', 'iROI', 'prefAz', 'prefEl', 'sigmaX');
for gi = azBeyondIdx
    fprintf('%-8d %-10.2f %-10.2f %-8.2f\n', gaussFitResults(gi).iROI, gaussFitResults(gi).x0, ...
        gaussFitResults(gi).y0, gaussFitResults(gi).sigmaX);
end

% Trusted boutons with El center beyond the stimulus-extended range
elBeyondIdx = trustedGI([gaussFitResults(trustedGI).isBeyondRangeEl]);
fprintf('\n--- Trusted boutons with El beyond range (n=%d) ---\n', numel(elBeyondIdx));
fprintf('%-8s %-10s %-10s %-8s\n', 'iROI', 'prefAz', 'prefEl', 'sigmaY');
for gi = elBeyondIdx
    fprintf('%-8d %-10.2f %-10.2f %-8.2f\n', gaussFitResults(gi).iROI, gaussFitResults(gi).x0, ...
        gaussFitResults(gi).y0, gaussFitResults(gi).sigmaY);
end

fprintf('\nFor reference, azRangeExt = [%.1f, %.1f], elRangeExt = [%.1f, %.1f]\n', ...
    azRangeExt(1), azRangeExt(2), elRangeExt(1), elRangeExt(2));
%%  Panel A: example bouton traces + Gaussian fits 
exampleBoutons = [1408 5518 1642 3406]; %3406
% exampleBoutons = [301 303 326 335];
exampleLabels  = {'Narrow RF', 'Medium RF', 'Broad RF', 'Broad RF 2'};

% Sanity check: confirm all example boutons are actually in the CURRENT trusted set before
% proceeding -- criteria have changed multiple times (r2Floor, degenerate-sigma fraction,
% beyond-range gating), so a previously-valid example bouton may no longer qualify.
missingExamples = exampleBoutons(~ismember(exampleBoutons, trustedIROI));
if ~isempty(missingExamples)
    error(['The following example bouton(s) are NOT in the current trusted set (n=%d): %s.\n' ...
           'Pick replacement(s) from trustedIROI, e.g. by re-running the "top N broadest/' ...
           'narrowest trusted boutons" check.'], numel(trustedGI), mat2str(missingExamples));
end

% allIROI here is built from gaussFitResultsTrusted, so `gi` found below indexes into
% gaussFitResultsTrusted specifically -- consistent with res = gaussFitResultsTrusted(gi) below.
allIROI = [gaussFitResultsTrusted.iROI];
exampleSet = nan(1, numel(exampleBoutons));
for i = 1:numel(exampleBoutons)
    gi = find(allIROI == exampleBoutons(i), 1);
    exampleSet(i) = gi;
end
nRows = numel(exampleSet);

figABC = figure('Color', 'w', 'Position', [100 50 900 1000], 'Name', 'Panel A: RF Profiles and Fits');
set(figABC, 'Visible', 'off');

leftMargin   = 0.08;
rightMargin  = 0.05;
topMargin    = 0.08;
bottomMargin = 0.08;
hGap         = 0.08;
vGap         = 0.06;

colWidth  = (1 - leftMargin - rightMargin - hGap) / 2;
rowHeight = (1 - topMargin - bottomMargin - (nRows-1)*vGap) / nRows;

[xq, yq] = meshgrid(linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100));

for r = 1:nRows
    gi    = exampleSet(r);
    res   = gaussFitResultsTrusted(gi);
    bIdx  = res.iROI;
    bData = allRFMapping(bIdx);
    yPos  = 1 - topMargin - r*rowHeight - (r-1)*vGap;

    trialMatrix = bData.baselineSubtracted;
    respVec = nan(numel(trialMatrix), 1);
    for posIdx = 1:numel(trialMatrix)
        trials = trialMatrix{posIdx};
        if ~isempty(trials)
            trials = double(trials);
            respVec(posIdx) = mean(mean(trials(:, respIdx), 2, 'omitnan'), 'omitnan');
        end
    end
    respGrid = reshape(respVec, size(AzGrid));

    traceShrinkFactor = 2;
    respWinForTraces  = [-1 4];

    boutonData = bData;
    boutonData.meanGridResponse     = respGrid;
    boutonData.meanTemporalResponse = bData.meanTemporalResponse;
    boutonData.meanBlankResponse    = bData.meanBlankResponse;
    boutonData.peakAmplitude        = max(respGrid(:)) * traceShrinkFactor;

    insetFrac = 0.85;
    padW = colWidth  * (1 - insetFrac) / 2;
    padH = rowHeight * (1 - insetFrac) / 2;

    axLine = axes('Position', [leftMargin+padW, yPos+padH, colWidth*insetFrac, rowHeight*insetFrac]);
    plotRFHeatmapWithTraces(axLine, boutonData, uAz, uEl_plot, timeVector, ...
        'Colormap', 'bone', 'Smooth', false, 'AxisMode', 'image', 'LabelUnits', 'deg', ...
        'PlotRespWin', true, 'RespWin', respWinForTraces);
    xlim(axLine, azLims);
    ylim(axLine, elLims);
    set(axLine, 'YDir', 'normal');
    set(axLine, 'XTickMode', 'auto', 'YTickMode', 'auto');
    title(axLine, sprintf('%s (Bouton %d)', exampleLabels{r}, bIdx), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    set(axLine, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 9);

    pFit = [res.A, res.x0, res.y0, res.sigmaX, res.sigmaY];
    fittedGrid = reshape(gaussFit2D(pFit, [xq(:), yq(:)]), size(xq));

    axFit = axes('Position', [leftMargin+colWidth+hGap+padW, yPos+padH, colWidth*insetFrac, rowHeight*insetFrac]);
    imagesc(axFit, linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100), fittedGrid);
    axis(axFit, 'image');
    colormap(axFit, 'bone');
    colorbar(axFit);
    xlim(axFit, azLims);
    ylim(axFit, elLims);
    set(axFit, 'YDir', 'normal');
    hold(axFit, 'on');
    plot(axFit, res.x0, res.y0, 'r+', 'MarkerSize', 6, 'LineWidth', 2);
    title(axFit, sprintf('2D Gaussian Fit  (R^2 = %.2f, boot CI_{low} = %.2f)', ...
        res.R2, res.bootR2LowerCI), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel(axFit, 'Azimuth (deg)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel(axFit, 'Elevation (deg)', 'FontName', 'Arial', 'FontSize', 9);
    set(axFit, 'XTickMode', 'auto', 'YTickMode', 'auto');
    set(axFit, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 9);
end

% outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\eg_rfs';
% if ~exist(outputDir, 'dir'), mkdir(outputDir); end
% saveFigureFormats(figABC, fullfile(outputDir, 'RF_profiles_and_2Dfits_grid'));
% close(figABC);



%%  Panel C: preferred azimuth/elevation histograms 


binWidth = 5;

% Azimuth — bin from actual data range, clipped to true screen limits
azBinLow  = -80;
azBinHigh = 20;
azEdges   = azBinLow:binWidth:azBinHigh;

% Elevation — same treatment
binWidthEL = 5;
elBinLow  = max(floor(min(prefEl)/binWidthEL)*binWidthEL, screenElLimits(1));
elBinHigh = min(ceil(max(prefEl)/binWidthEL)*binWidthEL, screenElLimits(2));
elEdges   = elBinLow:binWidthEL:elBinHigh;

figC = figure('Position', [100, 100, 600, 280]);

% --- Azimuth ---
axC1 = subplot(1,2,1);
azCounts  = histcounts(prefAz, azEdges);
azPct     = 100 * azCounts / numel(prefAz);
histogram(axC1, 'BinEdges', azEdges, 'BinCounts', azPct, 'FaceColor', [0.3 0.5 0.8]);
xlabel(axC1, 'Preferred Azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(axC1, '% of boutons', 'FontName', 'Arial', 'FontSize', 9);
title(axC1, sprintf('n = %d trusted', length(find(trustedMask))), 'FontName', 'Arial', 'FontSize', 9);
set(axC1, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8, 'XTick', azEdges);
drawnow;
defaultAxesProperties(axC1, true);

% --- Elevation ---
axC2 = subplot(1,2,2);
elCounts  = histcounts(prefEl, elEdges);
elPct     = 100 * elCounts / numel(prefEl);
histogram(axC2, 'BinEdges', elEdges, 'BinCounts', elPct, 'FaceColor', [0.8 0.4 0.3]);
xlabel(axC2, 'Preferred Elevation (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(axC2, '% of boutons', 'FontName', 'Arial', 'FontSize', 9);
title(axC2, sprintf('n = %d trusted', length(find(trustedMask))), 'FontName', 'Arial', 'FontSize', 9);
set(axC2, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8, 'XTick', elEdges);
drawnow;
defaultAxesProperties(axC2, true);

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\pref_az_el';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figC, 'Visible', 'off');
saveFigureFormats(figC, fullfile(outputDir, 'trusted_pref_az_el'));

%% Panel D: within-FOV deviation (az/el separate)
% Uses trustedGI/trustedIROI, defined once at the top -- NOT re-derived here.
trustedAz      = [gaussFitResults(trustedGI).x0];
trustedEl      = [gaussFitResults(trustedGI).y0];
trustedSession = sessionLabels(trustedIROI);

uniqueSessions = unique(trustedSession);

diffAz = nan(size(trustedAz));
diffEl = nan(size(trustedEl));

for si = 1:numel(uniqueSessions)
    sessMask = strcmp(trustedSession, uniqueSessions{si});
    if sum(sessMask) < 2
        continue; % need at least 2 boutons to define a meaningful FOV mean
    end
    diffAz(sessMask) = abs(trustedAz(sessMask) - mean(trustedAz(sessMask)));
    diffEl(sessMask) = abs(trustedEl(sessMask) - mean(trustedEl(sessMask)));
end

validIdx = ~isnan(diffAz);
nValid = sum(validIdx);

fprintf('Panel D: n = %d boutons (from FOVs with >=2 trusted boutons)\n', nValid);
fprintf('mean |dAz| = %.2f +/- %.2f\n', mean(diffAz(validIdx)), std(diffAz(validIdx)));
fprintf('mean |dEl| = %.2f +/- %.2f\n', mean(diffEl(validIdx)), std(diffEl(validIdx)));

figD = figure('Position', [100, 100, 600, 280]);

binWidthD = 5; % degrees

axD1 = subplot(1,2,1);
edgesAz = 0 : binWidthD : ceil(max(diffAz(validIdx))/binWidthD)*binWidthD;
countsAz = histcounts(diffAz(validIdx), edgesAz);
pctAz = 100 * countsAz / nValid;
centersAz = edgesAz(1:end-1) + diff(edgesAz)/2;
bar(centersAz, pctAz, 1, 'FaceColor', [0.3 0.5 0.8]);
xlabel('|\Delta Azimuth| from FOV mean (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel('% of boutons', 'FontName', 'Arial', 'FontSize', 9);
title(sprintf('mean = %.2f \\pm %.2f', mean(diffAz(validIdx)), std(diffAz(validIdx))), 'FontName', 'Arial', 'FontSize', 9);
set(axD1, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8, 'XTick', edgesAz);
drawnow;
defaultAxesProperties(axD1, true);

axD2 = subplot(1,2,2);
edgesEl = 0 : binWidthD : ceil(max(diffEl(validIdx))/binWidthD)*binWidthD;
countsEl = histcounts(diffEl(validIdx), edgesEl);
pctEl = 100 * countsEl / nValid;
centersEl = edgesEl(1:end-1) + diff(edgesEl)/2;
bar(centersEl, pctEl, 1, 'FaceColor', [0.8 0.4 0.3]);
xlabel('|\Delta Elevation| from FOV mean (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel('% of boutons', 'FontName', 'Arial', 'FontSize', 9);
title(sprintf('mean = %.2f \\pm %.2f', mean(diffEl(validIdx)), std(diffEl(validIdx))), 'FontName', 'Arial', 'FontSize', 9);
set(axD2, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8, 'XTick', edgesEl);
drawnow;
defaultAxesProperties(axD2, true);

%%  Panel D (alt): joint 2D density of within-FOV deviation 
% Reuses trustedAz/trustedEl/trustedSession/uniqueSessions from the section above -- no need to
% redefine trustedGI/trustedIROI again.
signedDiffAz = nan(size(trustedAz));
signedDiffEl = nan(size(trustedEl));

for si = 1:numel(uniqueSessions)
    sessMask = strcmp(trustedSession, uniqueSessions{si});
    if sum(sessMask) < 2
        continue;
    end
    signedDiffAz(sessMask) = trustedAz(sessMask) - mean(trustedAz(sessMask));
    signedDiffEl(sessMask) = trustedEl(sessMask) - mean(trustedEl(sessMask));
end

validIdx2D = ~isnan(signedDiffAz);
nValid2D = sum(validIdx2D);
fprintf('Panel D (joint density): n = %d boutons\n', nValid2D);

figD2 = figure('Position', [100, 100, 450, 400]);
axD2j = axes;

binWidth2D = 5; % degrees
maxAbsAz = ceil(max(abs(signedDiffAz(validIdx2D)))/binWidth2D)*binWidth2D;
maxAbsEl = ceil(max(abs(signedDiffEl(validIdx2D)))/binWidth2D)*binWidth2D;

edgesAz2D = -maxAbsAz : binWidth2D : maxAbsAz;
edgesEl2D = -maxAbsEl : binWidth2D : maxAbsEl;

histogram2(signedDiffAz(validIdx2D), signedDiffEl(validIdx2D), edgesAz2D, edgesEl2D, ...
    'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
colormap(axD2j, 'parula');
cb = colorbar;
cb.Label.String = 'Number of boutons';

xlabel('\Delta Azimuth from FOV mean (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel('\Delta Elevation from FOV mean (\circ)', 'FontName', 'Arial', 'FontSize', 9);
title(sprintf('Joint deviation from FOV mean (n=%d)', nValid2D), 'FontName', 'Arial', 'FontSize', 10);
axis(axD2j, 'equal');
set(axD2j, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
hold on;
xline(0, '--w', 'LineWidth', 1);
yline(0, '--w', 'LineWidth', 1);
drawnow;

%% 
% Reuses trustedGI, defined once at the top.
figCov = figure('Position', [100, 100, 650, 500], 'Color', 'w');
axCov = axes;
hold(axCov, 'on');

t_ellipse = linspace(0, 2*pi, 60);

for gi = trustedGI
    res = gaussFitResults(gi);
    xEllipse = res.x0 + res.sigmaX * cos(t_ellipse);
    yEllipse = res.y0 + res.sigmaY * sin(t_ellipse);
    patch(axCov, xEllipse, yEllipse, [0.2 0.3 0.6], ...
        'FaceAlpha', 0.06, 'EdgeColor', [0.2 0.3 0.6], 'EdgeAlpha', 0.15, 'LineWidth', 0.5);
end

xlabel(axCov, 'Azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 10);
ylabel(axCov, 'Elevation (\circ)', 'FontName', 'Arial', 'FontSize', 10);
title(axCov, sprintf('Visual Field Population Coverage Map (n = %d trusted)', numel(trustedGI)), ...
    'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');

axis(axCov, 'equal');
xlim(axCov, azRangeExt);
ylim(axCov, elRangeExt);
set(axCov, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 9);
drawnow;
%%
figSig = figure('Position', [100, 100, 900, 400], 'Color', 'w');

subplot(1,2,1);
histogram([gaussFitResults(trustedGI).sigmaX], 20, 'FaceColor', [0.2 0.3 0.6], 'FaceAlpha', 0.6);
hold on;
xline(degenerateSigmaFrac * range(azRange), 'r--', 'LineWidth', 1.5, 'Label', 'degenerate cutoff');
xlabel('sigmaX (Azimuth, deg)'); ylabel('Count');
title('Trusted sigmaX distribution');

subplot(1,2,2);
histogram([gaussFitResults(trustedGI).sigmaY], 20, 'FaceColor', [0.2 0.3 0.6], 'FaceAlpha', 0.6);
hold on;
xline(degenerateSigmaFrac * range(elRange), 'r--', 'LineWidth', 1.5, 'Label', 'degenerate cutoff');
xlabel('sigmaY (Elevation, deg)'); ylabel('Count');
title('Trusted sigmaY distribution');

%% within session 
% Mean absolute difference from each session's (FOV's) own mean position
madAzPerBouton = nan(size(prefAz));
madElPerBouton = nan(size(prefEl));

for s = 1:nSessions
    thisSessMask = strcmp(trustedSessLabels, uniqueSessions{s});
    sessAz = prefAz(thisSessMask);
    sessEl = prefEl(thisSessMask);
    
    madAzPerBouton(thisSessMask) = abs(sessAz - mean(sessAz));
    madElPerBouton(thisSessMask) = abs(sessEl - mean(sessEl));
end

% Combine az/el into a single 2D distance-from-FOV-mean, matching their single reported number
% (if they combined az+el into one distance rather than reporting separately)
madCombined = sqrt(madAzPerBouton.^2 + madElPerBouton.^2);

fprintf('Mean absolute difference from FOV mean:\n');
fprintf('  Azimuth:  %.2f +/- %.2f (SEM)\n', mean(madAzPerBouton,'omitnan'), std(madAzPerBouton,'omitnan')/sqrt(sum(~isnan(madAzPerBouton))));
fprintf('  Elevation: %.2f +/- %.2f (SEM)\n', mean(madElPerBouton,'omitnan'), std(madElPerBouton,'omitnan')/sqrt(sum(~isnan(madElPerBouton))));
fprintf('  Combined (2D distance): %.2f +/- %.2f (SEM)\n', mean(madCombined,'omitnan'), std(madCombined,'omitnan')/sqrt(sum(~isnan(madCombined))));