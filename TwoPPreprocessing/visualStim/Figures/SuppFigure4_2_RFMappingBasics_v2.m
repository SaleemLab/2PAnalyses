
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
analyse_RFBoutons_current % change session and run for bouton 
% then save the fit results
gaussFitResults_boutons = gaussFitResults; 

analyse_RFBoutons_current % change ids and run for somas 
gaussFitResults = gaussFitResults;



%% run this to save the .mat for each image type
save_RF_summary
%% Bootstrap robustness: is the (1-ciLevel) lower percentile of each fit's bootstrap R^2
% SOMAS
% distribution itself >= trustR2Threshold?
trustR2Threshold = 0.1;
ciLevel = 0.95;
lowerPct = 100 * (1 - ciLevel);

isRobust = false(1, numel(gaussFitResults));
for gi = 1:numel(gaussFitResults)
    validBoot = gaussFitResults(gi).bootR2(~isnan(gaussFitResults(gi).bootR2));
    if isempty(validBoot)
        continue;
    end
    isRobust(gi) = prctile(validBoot, lowerPct) >= trustR2Threshold;
end

% Trusted set: R^2 + bootstrap robustness + non-degenerate sigma + center not beyond range
% isDegenerateSigma already implements the >40% of tested range check (both axes).
% isBeyondRange already implements the stimulus-extended, screen-clipped center check --
% confirmed on this dataset that the raw-screen-limits check adds ZERO additional exclusions
% beyond isBeyondRange (182 vs 182, 0 beyond-screen-only), so it's redundant and omitted here.
trustedMask = ([gaussFitResults.R2] >= trustR2Threshold) & isRobust ...
    & ~[gaussFitResults.isDegenerateSigma] & ~[gaussFitResults.isBeyondRange];


trustedGI   = find(trustedMask);                    % indices into the FULL gaussFitResults
trustedIROI = [gaussFitResults(trustedGI).iROI];     % actual bouton iROI numbers, trusted only
nTrusted    = numel(trustedGI);

fprintf('Using %d / %d trusted boutons for all figure panels.\n', nTrusted, numel(gaussFitResults));
gaussFitResultsTrusted = gaussFitResults(trustedGI);
%% boutons

isRobustB = false(1, numel(gaussFitResults_boutons));
for gi = 1:numel(gaussFitResults_boutons)
    validBoot = gaussFitResults_boutons(gi).bootR2(~isnan(gaussFitResults_boutons(gi).bootR2));
    if isempty(validBoot)
        continue;
    end
    isRobustB(gi) = prctile(validBoot, lowerPct) >= trustR2Threshold;
end

% Trusted set: R^2 + bootstrap robustness + non-degenerate sigma + center not beyond range
% isDegenerateSigma already implements the >40% of tested range check (both axes).
% isBeyondRange already implements the stimulus-extended, screen-clipped center check --
% confirmed on this dataset that the raw-screen-limits check adds ZERO additional exclusions
% beyond isBeyondRange (182 vs 182, 0 beyond-screen-only), so it's redundant and omitted here.
trustedMaskB = ([gaussFitResults_boutons.R2] >= trustR2Threshold) & isRobustB ...
    & ~[gaussFitResults_boutons.isDegenerateSigma] & ~[gaussFitResults_boutons.isBeyondRange];


trustedGIB   = find(trustedMaskB);                    % indices into the FULL gaussFitResults
trustedIROIB = [gaussFitResults_boutons(trustedGIB).iROI];     % actual bouton iROI numbers, trusted only
nTrustedB    = numel(trustedGIB);

fprintf('Using %d / %d trusted boutons for all figure panels.\n', nTrustedB, numel(gaussFitResults_boutons));
gaussFitResultsTrustedB = gaussFitResults_boutons(trustedGIB);
%%
sigAzTrusted = [gaussFitResults(trustedMask).sigmaX];
sigElTrusted = [gaussFitResults(trustedMask).sigmaY];
prefAz       = [gaussFitResults(trustedMask).x0];
prefEl       = [gaussFitResults(trustedMask).y0];

fprintf('sigAz  = %.2f +/- %.2f (SD), n = %d\n', mean(sigAzTrusted), std(sigAzTrusted), numel(sigAzTrusted));
fprintf('sigEl  = %.2f +/- %.2f (SD), n = %d\n', mean(sigElTrusted), std(sigElTrusted), numel(sigElTrusted));
fprintf('prefAz = %.2f +/- %.2f (SD), n = %d\n', mean(prefAz), std(prefAz), numel(prefAz));
fprintf('prefEl = %.2f +/- %.2f (SD), n = %d\n', mean(prefEl), std(prefEl), numel(prefEl));



sigAzTrustedB = [gaussFitResults_boutons(trustedMaskB).sigmaX];
sigElTrustedB = [gaussFitResults_boutons(trustedMaskB).sigmaY];
prefAzB       = [gaussFitResults_boutons(trustedMaskB).x0];
prefElB       = [gaussFitResults_boutons(trustedMaskB).y0];

fprintf('Bouton sigAz  = %.2f +/- %.2f (SD), n = %d\n', mean(sigAzTrustedB), std(sigAzTrustedB), numel(sigAzTrustedB));
fprintf('Bouton sigEl  = %.2f +/- %.2f (SD), n = %d\n', mean(sigElTrustedB), std(sigElTrustedB), numel(sigElTrustedB));
fprintf('Bouton prefAz = %.2f +/- %.2f (SD), n = %d\n', mean(prefAzB), std(prefAzB), numel(prefAzB));
fprintf('Bouton prefEl = %.2f +/- %.2f (SD), n = %d\n', mean(prefAzB), std(prefAzB), numel(prefAzB));

%%  Panel A: example bouton traces + Gaussian fits  % somas 
exampleBoutons = [130 1 78]; 
exampleLabels  = {'Narrow RF', 'Medium RF', 'Broad RF', 'Broad RF 2'};

% Local axis limits, defined HERE rather than assumed to exist in the workspace from the fitting
% script. ONE SHARED range for BOTH panels -- matching the reference figure (Timplalexi et al.),
% which uses identical axis extents on the raw-data and fitted-model panels, not an extended
% range for the fit. Tight to the actual tested grid (+/- half a grid step), NOT the
% stimulus-footprint/screen-extended range -- deliberately reverting that earlier choice to match
% the published convention. Consequence (known and accepted): a fitted peak sitting right at the
% edge of the tested grid will sit flush against the plot wall, same as in the reference figure.
azStepLocal = mean(diff(sort(uAz)));
elStepLocal = mean(diff(sort(uEl_plot)));
sharedLimsAz = [min(uAz) - azStepLocal/2, max(uAz) + azStepLocal/2];
sharedLimsEl = [min(uEl_plot) - elStepLocal/2, max(uEl_plot) + elStepLocal/2];

% ALLOW PLOTTING ANY BOUTON, not just the current trusted set. Boutons not in trustedIROI are
% still plotted (indexed from the FULL gaussFitResults, since gaussFitResultsTrusted doesn't
% contain them), but flagged with a warning here and an explicit "(NOT TRUSTED)" tag in each
% panel's title -- so an untrusted example is still visibly distinguishable, never silently
% presented as if it passed the trust criteria.
notTrustedExamples = exampleBoutons(~ismember(exampleBoutons, trustedIROI));
if ~isempty(notTrustedExamples)
    warning(['The following example bouton(s) are NOT in the current trusted set (n=%d): %s.\n' ...
             'They will still be plotted (from the full gaussFitResults), but are tagged ' ...
             '"(NOT TRUSTED)" in their panel titles below.'], numel(trustedGI), mat2str(notTrustedExamples));
end

% Index into the FULL gaussFitResults (not gaussFitResultsTrusted) so any requested bouton --
% trusted or not -- can be found, as long as it was fit at all.
allIROIFull = [gaussFitResults.iROI];
exampleSet = nan(1, numel(exampleBoutons));
for i = 1:numel(exampleBoutons)
    gi = find(allIROIFull == exampleBoutons(i), 1);
    if isempty(gi)
        error('Bouton %d was not found in gaussFitResults at all (never fit, or excluded before the screening floor). Cannot plot.', exampleBoutons(i));
    end
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

[xq, yq] = meshgrid(linspace(sharedLimsAz(1), sharedLimsAz(2), 100), linspace(sharedLimsEl(1), sharedLimsEl(2), 100));

for r = 1:nRows
    gi    = exampleSet(r);
    res   = gaussFitResults(gi);
    bIdx  = res.iROI;
    bData = allRFMapping(bIdx);
    yPos  = 1 - topMargin - r*rowHeight - (r-1)*vGap;
    isThisTrusted = ismember(bIdx, trustedIROI);

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
    xlim(axLine, sharedLimsAz);
    ylim(axLine, sharedLimsEl);
    set(axLine, 'YDir', 'normal');
    % Ticks at the ACTUAL tested grid centers (uAz/uEl_plot), not auto-generated values --
    % labels what was really stimulated, not an arbitrary/extended axis range.
    set(axLine, 'XTick', sort(uAz), 'YTick', sort(uEl_plot));
    trustTag = '';
    if ~isThisTrusted
        trustTag = ' (NOT TRUSTED)';
    end
    title(axLine, sprintf('%s (Bouton %d)%s', exampleLabels{r}, bIdx, trustTag), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    set(axLine, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 9);

    pFit = [res.A, res.x0, res.y0, res.sigmaX, res.sigmaY];
    fittedGrid = reshape(gaussFit2D(pFit, [xq(:), yq(:)]), size(xq));

    axFit = axes('Position', [leftMargin+colWidth+hGap+padW, yPos+padH, colWidth*insetFrac, rowHeight*insetFrac]);
    imagesc(axFit, linspace(sharedLimsAz(1), sharedLimsAz(2), 100), linspace(sharedLimsEl(1), sharedLimsEl(2), 100), fittedGrid);
    axis(axFit, 'image');
    colormap(axFit, 'bone');
    colorbar(axFit);
    xlim(axFit, sharedLimsAz);
    ylim(axFit, sharedLimsEl);
    set(axFit, 'YDir', 'normal');
    hold(axFit, 'on');
    plot(axFit, res.x0, res.y0, 'r+', 'MarkerSize', 6, 'LineWidth', 2);
    title(axFit, sprintf('2D Gaussian Fit  (R^2 = %.2f, boot CI_{low} = %.2f)', ...
        res.R2, res.bootR2LowerCI), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel(axFit, 'Azimuth (deg)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel(axFit, 'Elevation (deg)', 'FontName', 'Arial', 'FontSize', 9);
    % Same tick convention as axLine: label the actual tested grid centers.
    set(axFit, 'XTick', sort(uAz), 'YTick', sort(uEl_plot));
    set(axFit, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 9);
end

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Supp_Section1_Fig4_1_VISp\eg_rfs';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figABC, fullfile(outputDir, 'RF_profiles_and_2Dfits_grid_boutons'));
close(figABC);



%% somas
exampleBoutons = [39 105 6]; %3406
% exampleBoutons = [301 303 326 335];
exampleLabels  = {'Narrow RF', 'Medium RF', 'Broad RF', 'Broad RF 2'};

% Local axis limits, defined HERE rather than assumed to exist in the workspace from the fitting
% script. ONE SHARED range for BOTH panels -- matching the reference figure (Timplalexi et al.),
% which uses identical axis extents on the raw-data and fitted-model panels, not an extended
% range for the fit. Tight to the actual tested grid (+/- half a grid step), NOT the
% stimulus-footprint/screen-extended range -- deliberately reverting that earlier choice to match
% the published convention. Consequence (known and accepted): a fitted peak sitting right at the
% edge of the tested grid will sit flush against the plot wall, same as in the reference figure.
azStepLocal = mean(diff(sort(uAz)));
elStepLocal = mean(diff(sort(uEl_plot)));
sharedLimsAz = [min(uAz) - azStepLocal/2, max(uAz) + azStepLocal/2];
sharedLimsEl = [min(uEl_plot) - elStepLocal/2, max(uEl_plot) + elStepLocal/2];

% ALLOW PLOTTING ANY BOUTON, not just the current trusted set. Boutons not in trustedIROI are
% still plotted (indexed from the FULL gaussFitResults, since gaussFitResultsTrusted doesn't
% contain them), but flagged with a warning here and an explicit "(NOT TRUSTED)" tag in each
% panel's title -- so an untrusted example is still visibly distinguishable, never silently
% presented as if it passed the trust criteria.
notTrustedExamples = exampleBoutons(~ismember(exampleBoutons, trustedIROI));
if ~isempty(notTrustedExamples)
    warning(['The following example bouton(s) are NOT in the current trusted set (n=%d): %s.\n' ...
             'They will still be plotted (from the full gaussFitResults), but are tagged ' ...
             '"(NOT TRUSTED)" in their panel titles below.'], numel(trustedGI), mat2str(notTrustedExamples));
end

% Index into the FULL gaussFitResults (not gaussFitResultsTrusted) so any requested bouton --
% trusted or not -- can be found, as long as it was fit at all.
allIROIFull = [gaussFitResults.iROI];
exampleSet = nan(1, numel(exampleBoutons));
for i = 1:numel(exampleBoutons)
    gi = find(allIROIFull == exampleBoutons(i), 1);
    if isempty(gi)
        error('Bouton %d was not found in gaussFitResults at all (never fit, or excluded before the screening floor). Cannot plot.', exampleBoutons(i));
    end
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

[xq, yq] = meshgrid(linspace(sharedLimsAz(1), sharedLimsAz(2), 100), linspace(sharedLimsEl(1), sharedLimsEl(2), 100));

for r = 1:nRows
    gi    = exampleSet(r);
    res   = gaussFitResults(gi);
    bIdx  = res.iROI;
    bData = allRFMapping(bIdx);
    yPos  = 1 - topMargin - r*rowHeight - (r-1)*vGap;
    isThisTrusted = ismember(bIdx, trustedIROI);

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
    xlim(axLine, sharedLimsAz);
    ylim(axLine, sharedLimsEl);
    set(axLine, 'YDir', 'normal');
    % Ticks at the ACTUAL tested grid centers (uAz/uEl_plot), not auto-generated values --
    % labels what was really stimulated, not an arbitrary/extended axis range.
    set(axLine, 'XTick', sort(uAz), 'YTick', sort(uEl_plot));
    trustTag = '';
    if ~isThisTrusted
        trustTag = ' (NOT TRUSTED)';
    end
    title(axLine, sprintf('%s (Bouton %d)%s', exampleLabels{r}, bIdx, trustTag), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    set(axLine, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 9);

    pFit = [res.A, res.x0, res.y0, res.sigmaX, res.sigmaY];
    fittedGrid = reshape(gaussFit2D(pFit, [xq(:), yq(:)]), size(xq));

    axFit = axes('Position', [leftMargin+colWidth+hGap+padW, yPos+padH, colWidth*insetFrac, rowHeight*insetFrac]);
    imagesc(axFit, linspace(sharedLimsAz(1), sharedLimsAz(2), 100), linspace(sharedLimsEl(1), sharedLimsEl(2), 100), fittedGrid);
    axis(axFit, 'image');
    colormap(axFit, 'bone');
    colorbar(axFit);
    xlim(axFit, sharedLimsAz);
    ylim(axFit, sharedLimsEl);
    set(axFit, 'YDir', 'normal');
    hold(axFit, 'on');
    plot(axFit, res.x0, res.y0, 'r+', 'MarkerSize', 6, 'LineWidth', 2);
    title(axFit, sprintf('2D Gaussian Fit  (R^2 = %.2f, boot CI_{low} = %.2f)', ...
        res.R2, res.bootR2LowerCI), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel(axFit, 'Azimuth (deg)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel(axFit, 'Elevation (deg)', 'FontName', 'Arial', 'FontSize', 9);
    % Same tick convention as axLine: label the actual tested grid centers.
    set(axFit, 'XTick', sort(uAz), 'YTick', sort(uEl_plot));
    set(axFit, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 9);
end

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Supp_Section1_Fig4_1_VISp\eg_rfs';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figABC, fullfile(outputDir, 'RF_profiles_and_2Dfits_grid_somas'));
close(figABC);



%% boutons vs somas 
binWidth = 5;
% Azimuth — bin from actual data range, clipped to true screen limits
azBinLow  = -80;
azBinHigh = 20;
azEdges   = azBinLow:binWidth:azBinHigh;
% Elevation — use combined soma+bouton range so both fit the same bins
binWidthEL = 5;
elBinLow  = max(floor(min([prefEl, prefElB])/binWidthEL)*binWidthEL, screenElLimits(1));
elBinHigh = min(ceil(max([prefEl, prefElB])/binWidthEL)*binWidthEL, screenElLimits(2));
elEdges   = elBinLow:binWidthEL:elBinHigh;
somaColor   = [0.3 0.5 0.8];
boutonColor = [0.9 0.5 0.1];
figC = figure('Position', [100, 100, 600, 280]);
% --- Azimuth ---
axC1 = subplot(1,2,1);
hold(axC1, 'on');
azCounts  = histcounts(prefAz, azEdges);
azPct     = 100 * azCounts / numel(prefAz);
histogram(axC1, 'BinEdges', azEdges, 'BinCounts', azPct, ...
'FaceColor', somaColor, 'FaceAlpha', 0.6, 'EdgeColor', somaColor*0.6, 'LineWidth', 1, 'DisplayName', 'Somas');
azCountsB = histcounts(prefAzB, azEdges);
azPctB    = 100 * azCountsB / numel(prefAzB);
histogram(axC1, 'BinEdges', azEdges, 'BinCounts', azPctB, ...
'FaceColor', boutonColor, 'FaceAlpha', 0.6, 'EdgeColor', boutonColor*0.6, 'LineWidth', 1, 'DisplayName', 'Boutons');
xlabel(axC1, 'Preferred Azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(axC1, '% of cells/boutons', 'FontName', 'Arial', 'FontSize', 9);
title(axC1, sprintf('Somas n=%d, Boutons n=%d', nTrusted, nTrustedB), 'FontName', 'Arial', 'FontSize', 9);
legend(axC1, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
set(axC1, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8, 'XTick', azEdges);
drawnow;
defaultAxesProperties(axC1, true);
% --- Elevation ---
axC2 = subplot(1,2,2);
hold(axC2, 'on');
elCounts  = histcounts(prefEl, elEdges);
elPct     = 100 * elCounts / numel(prefEl);
histogram(axC2, 'BinEdges', elEdges, 'BinCounts', elPct, ...
'FaceColor', somaColor, 'FaceAlpha', 0.6, 'EdgeColor', somaColor*0.6, 'LineWidth', 1, 'DisplayName', 'Somas');
elCountsB = histcounts(prefElB, elEdges);
elPctB    = 100 * elCountsB / numel(prefElB);
histogram(axC2, 'BinEdges', elEdges, 'BinCounts', elPctB, ...
'FaceColor', boutonColor, 'FaceAlpha', 0.6, 'EdgeColor', boutonColor*0.6, 'LineWidth', 1, 'DisplayName', 'Boutons');
xlabel(axC2, 'Preferred Elevation (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(axC2, '% of cells/boutons', 'FontName', 'Arial', 'FontSize', 9);
title(axC2, sprintf('Somas n=%d, Boutons n=%d', nTrusted, nTrustedB), 'FontName', 'Arial', 'FontSize', 9);
legend(axC2, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
set(axC2, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8, 'XTick', elEdges);
drawnow;
defaultAxesProperties(axC2, true);
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Supp_Section1_Fig4_1_VISp\pref_az_el';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figC, 'Visible', 'off');
saveFigureFormats(figC, fullfile(outputDir, 'trusted_pref_az_el_somas_vs_boutons'));

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
title(axC1, sprintf('n = %d trusted', nTrusted), 'FontName', 'Arial', 'FontSize', 9);
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
title(axC2, sprintf('n = %d trusted', nTrusted), 'FontName', 'Arial', 'FontSize', 9);
set(axC2, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8, 'XTick', elEdges);
drawnow;
defaultAxesProperties(axC2, true);

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Supp_Section1_Fig4_1_VISp\pref_az_el';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figC, 'Visible', 'off');
saveFigureFormats(figC, fullfile(outputDir, 'trusted_pref_az_el_somas'));

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
    patch(axCov, xEllipse, yEllipse, 'k', ...
        'FaceAlpha', 0.06, 'EdgeColor', 'k', 'EdgeAlpha', 0.15, 'LineWidth', 0.5);
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