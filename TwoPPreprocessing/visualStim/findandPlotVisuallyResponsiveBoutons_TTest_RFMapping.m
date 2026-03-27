function findandPlotVisuallyResponsiveBoutons_TTest_RFMapping(sessionFileInfo, stimName)
% 1. Setup and Pathing
iStim = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
load(sessionFileInfo.stimFiles(iStim).Response, 'response');
psth = response.psthData;
stimVs = vertcat(psth.stimValue);
nROI = size(psth(1).alignedResponses, 1);

% Identify Grid vs Blank ([200, 0])
blankIdx = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
gridMask = stimVs(:,1) ~= 200;
gridPSTH = psth(gridMask);
gridStim = stimVs(gridMask, :);
uAz = sort(unique(gridStim(:,1)), 'ascend');  
uEl_plot = sort(unique(gridStim(:,2)), 'descend'); 
uEl_bins = sort(unique(gridStim(:,2)), 'ascend'); 
nAz = length(uAz); nEl = length(uEl_plot);

% Half-pixel offsets for centering 
dAz = 20; dEl = 20; 
if nAz > 1, dAz = abs(uAz(2) - uAz(1)); end
if nEl > 1, dEl = abs(uEl_plot(1) - uEl_plot(2)); end
azLimFull = [min(uAz) - dAz/2, max(uAz) + dAz/2];
elLimFull = [min(uEl_plot) - dEl/2, max(uEl_plot) + dEl/2];

[~, fileName] = fileparts(sessionFileInfo.stimFiles(iStim).Response);
outDir  = fullfile(sessionFileInfo.Directories.save_folder,'Figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
pdfPath = fullfile(outDir, [fileName 'ResponsiveBoutonReport_UsingTTest_RFMapping.pdf']);
if exist(pdfPath, 'file'), delete(pdfPath); end

respWin = [0.5 2]; baseWin = [-2 0]; sigma_smooth = 0.8; alphaThresh = 0.05;

%% 2. Process Responsiveness (T-Test)
responsiveROIs = [];
tVecBlank = psth(blankIdx).timeVector(:);
bMask = tVecBlank >= respWin(1) & tVecBlank <= respWin(2);

for i = 1:nROI
    rawBlank = psth(blankIdx).alignedResponses(i, :, :);
    blankDist = squeeze(mean(rawBlank(:, bMask, :), 2, 'omitnan'));
    pVals = ones(numel(gridPSTH), 1);
    for k = 1:numel(gridPSTH)
        tVecK = gridPSTH(k).timeVector;
        kMask = tVecK >= respWin(1) & tVecK <= respWin(2);
        posDist = squeeze(mean(gridPSTH(k).alignedResponses(i, kMask, :), 2, 'omitnan'));
        if numel(posDist) > 1 && numel(blankDist) > 1
            [~, pVals(k)] = ttest2(posDist, blankDist, 'Tail', 'right');
        end
    end
    if any(pVals < alphaThresh), responsiveROIs = [responsiveROIs; i]; end
end

%% 3. Summary Cover Page
summaryFig = figure('Color', 'w', 'Position', [100 100 600 800], 'Visible', 'off');
annotation('textbox', [0.1 0.5 0.8 0.4], 'String', ...
    {['Session: ' stimName], ['Visually Responsive: ' num2str(length(responsiveROIs))], ...
     ['% Responsive: ' num2str(round(length(responsiveROIs)/nROI*100,1)) '%']}, ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
exportgraphics(summaryFig, pdfPath); close(summaryFig);

%% 4. ROI Pages & Center Collection
allCenters = [];
for p = 1:length(responsiveROIs)
    roiIdx = responsiveROIs(p);
    rfMatrix = nan(nEl, nAz); traceMax = 1e-6; 
    for k = 1:numel(gridPSTH)
        mu = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
        tV = gridPSTH(k).timeVector(:);
        mu = mu(:) - mean(mu(tV >= baseWin(1) & tV < baseWin(2)), 'omitnan');
        traceMax = max(traceMax, max(mu));
        r = find(uEl_plot == gridStim(k,2), 1); c = find(uAz == gridStim(k,1), 1);
        rfMatrix(r, c) = mean(mu(tV >= respWin(1) & tV <= respWin(2)), 'omitnan');
    end
    
    % Store the peak location for the population map
    [~, mI] = max(rfMatrix(:));
    [rPeak, cPeak] = ind2sub(size(rfMatrix), mI);
    allCenters = [allCenters; uAz(cPeak), uEl_plot(rPeak)];

    fig = figure('Color', 'w', 'Position', [100 100 700 600], 'Visible', 'off');
    ax1 = axes('Position', [0.15 0.15 0.7 0.7]);
    imagesc(uAz, uEl_plot, imgaussfilt(rfMatrix, sigma_smooth)); hold on; colormap(gray);
    set(ax1, 'YDir', 'normal', 'CLim', [0, max(rfMatrix(:)) + 1e-6]);
    xlim(azLimFull); ylim(elLimFull); colorbar;
    
    % Scaling to center traces [cite: 5, 21, 33]
    vS = dEl * 0.4; hS = dAz * 0.85;
    for k = 1:numel(gridPSTH)
        tr = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
        tr = tr(:) - mean(tr(tV >= baseWin(1) & tV < baseWin(2)), 'omitnan');
        sY = (tr / traceMax * vS) + gridStim(k,2);
        tNorm = (tV - min(tV)) / (max(tV) - min(tV));
        plot((tNorm-0.5)*hS + gridStim(k,1), sY, 'r', 'LineWidth', 1.2);
    end
    title(['ROI ' num2str(roiIdx) ' Heatmap + PSTHs']);
    exportgraphics(fig, pdfPath, 'Append', true); close(fig);
end

%% 5. POPULATION DENSITY SUMMARY
if ~isempty(allCenters)
    popFig = figure('Color', 'w', 'Position', [100 100 700 600], 'Visible', 'off');
    % Bin centers using monotonic edges
    [counts, ~, ~] = histcounts2(allCenters(:,1), allCenters(:,2), ...
        [uAz; max(uAz)+dAz], [uEl_bins; max(uEl_bins)+dEl]);
    
    imagesc(uAz, uEl_plot, flipud(counts')); colormap(flipud(bone)); 
    set(gca, 'YDir', 'normal'); xlim(azLimFull); ylim(elLimFull); colorbar;
    
    % Label counts in the middle of each grid cell 
    for r = 1:nEl
        for c = 1:nAz
            val = counts(c, nEl-r+1);
            if val > 0
                text(uAz(c), uEl_plot(r), num2str(val), 'Color', 'r', ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'center'); 
            end
        end
    end
    title('Population RF Density (Bouton Count)');
    xlabel('Azimuth (°)'); ylabel('Elevation (°)');
    exportgraphics(popFig, pdfPath, 'Append', true); close(popFig);
end
end