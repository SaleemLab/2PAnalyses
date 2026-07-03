%% Panel A (redesigned): per-bouton Gaussian fit -> principled example selection
% Requires: fitGaussian2D.m, plotRFHeatmapWithTraces.m, plotRFHeatmapOnly.m on path
% Requires: allRFMapping, respIdxList, unrespIdxList, uAz, uEl_plot, timeVector,
%           RFMappingMetadata.respWin already built (as in your pooling script)

%% Step 1: fit a 2D Gaussian to every responsive bouton, store R2 + FWHM
nResp = numel(respIdxList);
fitR2            = nan(nResp, 1);
fitFWHM_az       = nan(nResp, 1);
fitFWHM_el       = nan(nResp, 1);
fitFWHM_combined = nan(nResp, 1);  % geometric mean of az/el FWHM, used for sorting

for k = 1:nResp
    b = allRFMapping(respIdxList(k));
    Z = b.meanGridResponse;

    try
        [fitFun, fwhmX, fwhmY] = fitGaussian2D(uAz, uEl_plot, Z);

        [Xg, Yg] = meshgrid(uAz, uEl_plot);
        Zhat = fitFun(Xg, Yg);
        validMask = ~isnan(Z);

        SSres = sum((Z(validMask) - Zhat(validMask)).^2);
        SStot = sum((Z(validMask) - mean(Z(validMask))).^2);
        r2 = 1 - SSres/SStot;

        fitR2(k)            = r2;
        fitFWHM_az(k)        = fwhmX;
        fitFWHM_el(k)        = fwhmY;
        fitFWHM_combined(k)  = sqrt(fwhmX * fwhmY);
    catch ME
        warning('Gaussian fit failed for bouton %d: %s', respIdxList(k), ME.message);
    end
end

% sanity check before choosing a threshold -- look for a bimodal split
figure('Name', 'Gaussian fit R2 distribution');
histogram(fitR2, 30);
xlabel('R^2 of 2D Gaussian fit'); ylabel('# boutons');
title('Check for a natural well-fit / poorly-fit split before setting r2Thresh');

%% Step 2: filter to well-fit boutons, pick narrow / medium / broad by percentile
r2Thresh = 0.5;  % <-- adjust based on the histogram above

wellFitIdxInResp = find(fitR2 > r2Thresh);
wellFitBoutonIdx = respIdxList(wellFitIdxInResp);
widths           = fitFWHM_combined(wellFitIdxInResp);

[sortedWidths, sortOrder] = sort(widths);
nWellFit = numel(sortedWidths);

pctTargets = [10 50 90];  % narrow, medium, broad
exampleRespIdxList = nan(1,3);
for i = 1:3
    targetRank = round(pctTargets(i)/100 * nWellFit);
    targetRank = max(1, min(nWellFit, targetRank));
    exampleRespIdxList(i) = wellFitBoutonIdx(sortOrder(targetRank));
end

fprintf('Selected responsive examples (bouton IDs): narrow=%d, medium=%d, broad=%d\n', ...
    exampleRespIdxList);
fprintf('Corresponding combined FWHM (deg): %.1f, %.1f, %.1f\n', ...
    sortedWidths(max(1,min(nWellFit,round(pctTargets/100*nWellFit)))));

%% Step 3: pick 2 non-responsive examples (unchanged from your existing approach)
unrespAmps       = [allRFMapping(unrespIdxList).peakAmplitude];
[~, medRankIdx]  = sort(unrespAmps);
exampleUnrespIdx1 = unrespIdxList(medRankIdx(round(numel(medRankIdx) * 0.25)));
exampleUnrespIdx2 = unrespIdxList(medRankIdx(round(numel(medRankIdx) * 0.75)));

%% Step 4: assemble the figure -- one smoothed (mean) heatmap + traces per example
% Single panel per bouton: smoothed mean map with traces superimposed at
% their correct Az/El position, axes extending to the true grid edges.

exampleSet    = [exampleRespIdxList, exampleUnrespIdx1, exampleUnrespIdx2];
exampleLabels = {'Narrow RF', 'Medium RF', 'Broad RF', 'Non-responsive', 'Non-responsive'};

nEx  = numel(exampleSet);
colW = 1 / nEx;

figA = figure('Color', 'w', 'Position', [50 50 nEx*550 750], 'Name', 'Panel A: RF examples');

for i = 1:nEx
    b = allRFMapping(exampleSet(i));
    xBase = (i-1)*colW;

    axPanel = axes('Position', [xBase+0.03 0.12 colW-0.06 0.75]);
    plotRFHeatmapWithTraces(axPanel, b, uAz, uEl_plot, timeVector, 'Smooth', true);
    title(axPanel, sprintf('%s (bouton %d)', exampleLabels{i}, exampleSet(i)), ...
        'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
end