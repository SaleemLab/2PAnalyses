function findandPlotVisuallyResponsiveBoutons_RFMapping(sessionFileInfo, stimName)
% This plots all visually response rois using t-test 

iStim = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
load(sessionFileInfo.stimFiles(iStim).Response, 'response');
psth = response.psthData;
stimVs = vertcat(psth.stimValue);
nROI = size(psth(1).alignedResponses, 1);


blankIdx = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
gridMask = stimVs(:,1) ~= 200;
gridPSTH = psth(gridMask);
gridStim = stimVs(gridMask, :);
uAz = sort(unique(gridStim(:,1)), 'ascend');  
uEl_plot = sort(unique(gridStim(:,2)), 'descend'); % For visual display
uEl_bins = sort(unique(gridStim(:,2)), 'ascend');  % For histcounts2 monotonicity
nAz = length(uAz); nEl = length(uEl_plot);

% Calculate half-pixel offsets to center coordinates in the grid squares 
dAz = 20; dEl = 20; 
if nAz > 1, dAz = abs(uAz(2) - uAz(1)); end
if nEl > 1, dEl = abs(uEl_plot(1) - uEl_plot(2)); end
azLimFull = [min(uAz) - dAz/2, max(uAz) + dAz/2];
elLimFull = [min(uEl_plot) - dEl/2, max(uEl_plot) + dEl/2];

[~, fileName] = fileparts(sessionFileInfo.stimFiles(iStim).Response);
outDir  = fullfile(sessionFileInfo.Directories.save_folder,'Figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
pdfPath = fullfile(outDir, [fileName 'ResponsiveBoutonReport_UsingANOVA_RFMapping.pdf']);
if exist(pdfPath, 'file'), delete(pdfPath); end

% Parameters from Timplalexi et al. 2025 
respWin = [0.5 2];  
baseWin = [-2 0];   
sigma_smooth = 0.8;

%% 2. Run Responsiveness Analysis (ANOVA + 2SD Rule) [cite: 2]
responsiveROIs = [];
for i = 1:nROI
    dataForAnova = []; 
    groupLabels = [];
    prefResp = -inf;
    
    for g = 1:numel(psth)
        tVec = psth(g).timeVector;
        wMask = tVec >= respWin(1) & tVec <= respWin(2);
        trialMeans = squeeze(mean(psth(g).alignedResponses(i, wMask, :), 2, 'omitnan'));
        dataForAnova = [dataForAnova; trialMeans(:)]; 
        groupLabels = [groupLabels; repmat(g, numel(trialMeans), 1)];
        prefResp = max(prefResp, mean(trialMeans, 'omitnan'));
    end
    
    pValANOVA = anova1(dataForAnova, groupLabels, 'off');
    
    tVecB = psth(blankIdx).timeVector;
    bMask = tVecB >= respWin(1) & tVecB <= respWin(2);
    blankTrialMeans = squeeze(mean(psth(blankIdx).alignedResponses(i, bMask, :), 2, 'omitnan'));
    
    % Paper criteria: ANOVA p<0.05 and peak > blank + 2SD 
    if pValANOVA < 0.05 && prefResp > (mean(blankTrialMeans, 'omitnan') + 2 * std(blankTrialMeans, 'omitnan'))
        responsiveROIs = [responsiveROIs; i];
    end
end

%% 3. GENERATE SUMMARY COVER PAGE
coverFig = figure('Color', 'w', 'Position', [100 100 600 800], 'Visible', 'off');
summaryText = { ...
    ['Session: ' stimName], ...
    ['Total Boutons: ' num2str(nROI)], ...
    ['Responsive Boutons: ' num2str(length(responsiveROIs))], ...
    ['% Responsive: ' num2str(round(length(responsiveROIs)/nROI*100, 1)) '%'], ...
    '', ...
    'Criteria: ANOVA (p<0.05) & Peak > Blank+2SD'};

annotation('textbox', [0.1 0.1 0.8 0.8], 'String', summaryText, ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

exportgraphics(coverFig, pdfPath); 
close(coverFig);

%% 4. Generate ROI Pages (Heatmap + PSTHs Only) 
allCenters = [];
for p = 1:length(responsiveROIs)
    roiIdx = responsiveROIs(p);
    
    rfMatrix = nan(nEl, nAz); 
    traceMax = 1e-6; 
    
    for k = 1:numel(gridPSTH)
        mu = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
        tV = gridPSTH(k).timeVector(:);
        mu = mu(:) - mean(mu(tV >= baseWin(1) & tV < baseWin(2)), 'omitnan');
        traceMax = max(traceMax, max(mu));
        r = find(uEl_plot == gridStim(k,2), 1); c = find(uAz == gridStim(k,1), 1);
        rfMatrix(r, c) = mean(mu(tV >= respWin(1) & tV <= respWin(2)), 'omitnan');
    end
    
    % Store peak location for population summary
    [~, mI] = max(rfMatrix(:));
    [rPeak, cPeak] = ind2sub(size(rfMatrix), mI);
    allCenters = [allCenters; uAz(cPeak), uEl_plot(rPeak)];

    fig = figure('Color', 'w', 'Position', [50 50 700 600], 'Visible', 'off');
    ax1 = axes('Position', [0.15 0.15 0.7 0.7]);
    
    imagesc(uAz, uEl_plot, imgaussfilt(rfMatrix, sigma_smooth)); hold on; 
    colormap(ax1, gray);
    set(ax1, 'YDir', 'normal', 'CLim', [0, max(rfMatrix(:)) + 1e-6]);
    xlim(azLimFull); ylim(elLimFull); colorbar;
    
    vS_base = dEl * 0.4; % 40% cell height
    hS = dAz * 0.85;     % 85% cell width
    
    for k = 1:numel(gridPSTH)
        tr = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
        tr = tr(:) - mean(tr(tV >= baseWin(1) & tV < baseWin(2)), 'omitnan');
        sY = (tr / traceMax * vS_base) + gridStim(k,2);
        tNorm = (tV - min(tV)) / (max(tV) - min(tV));
        sX = (tNorm - 0.5) * hS + gridStim(k,1);
        plot(ax1, sX, sY, 'r', 'LineWidth', 1.2);
    end
    
    title(['ROI ' num2str(roiIdx) ' Heatmap + PSTHs']);
    xlabel('azimuth (°)'); ylabel('elevation (°)');
    
    exportgraphics(fig, pdfPath, 'Append', true);
    close(fig);
end

%% 5. GENERATE POPULATION DENSITY SUMMARY
if ~isempty(allCenters)
    popFig = figure('Color', 'w', 'Position', [100 100 700 600], 'Visible', 'off');
    [counts, ~, ~] = histcounts2(allCenters(:,1), allCenters(:,2), ...
        [uAz; max(uAz)+dAz], [uEl_bins; max(uEl_bins)+dEl]);
    
    imagesc(uAz, uEl_plot, flipud(counts')); colormap(flipud(bone)); 
    set(gca, 'YDir', 'normal'); xlim(azLimFull); ylim(elLimFull); colorbar;
    
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
    exportgraphics(popFig, pdfPath, 'Append', true); 
    close(popFig);
end
end


% function findVisuallyResponsiveBoutons_RFMapping(sessionFileInfo, stimName)
% % 1. Setup and Pathing
% iStim = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
% load(sessionFileInfo.stimFiles(iStim).Response, 'response');
% psth = response.psthData;
% stimVs = vertcat(psth.stimValue);
% nROI = size(psth(1).alignedResponses, 1);
% 
% % Identify Grid vs Blank
% blankIdx = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
% gridMask = stimVs(:,1) ~= 200;
% gridPSTH = psth(gridMask);
% gridStim = stimVs(gridMask, :);
% uAz = sort(unique(gridStim(:,1)), 'ascend');  
% uEl = sort(unique(gridStim(:,2)), 'descend');
% nAz = length(uAz); nEl = length(uEl);
% 
% % Calculate half-pixel offsets to center coordinates in the grid squares [cite: 4]
% dAz = 20; dEl = 20; 
% if nAz > 1, dAz = abs(uAz(2) - uAz(1)); end
% if nEl > 1, dEl = abs(uEl(1) - uEl(2)); end
% azLimFull = [min(uAz) - dAz/2, max(uAz) + dAz/2];
% elLimFull = [min(uEl) - dEl/2, max(uEl) + dEl/2];
% 
% [~, fileName] = fileparts(sessionFileInfo.stimFiles(iStim).Response);
% outDir  = fullfile(sessionFileInfo.Directories.save_folder,'Figures');
% if ~exist(outDir,'dir'), mkdir(outDir); end
% pdfPath = fullfile(outDir, [fileName 'ResponsiveBoutonReport_UsingANOVA_RFMapping.pdf']);
% if exist(pdfPath, 'file'), delete(pdfPath); end
% 
% % Parameters from Timplalexi et al. 2025 [cite: 2]
% respWin = [0.5 2];  
% baseWin = [-2 0];   
% sigma_smooth = 0.8;
% 
% %% 2. Run Responsiveness Analysis (ANOVA + 2SD Rule) [cite: 2]
% responsiveROIs = [];
% for i = 1:nROI
%     % Standardized variable name: dataForAnova [cite: 5]
%     dataForAnova = []; 
%     groupLabels = [];
%     prefResp = -inf;
% 
%     for g = 1:numel(psth)
%         tVec = psth(g).timeVector;
%         wMask = tVec >= respWin(1) & tVec <= respWin(2);
%         trialMeans = squeeze(mean(psth(g).alignedResponses(i, wMask, :), 2, 'omitnan'));
%         dataForAnova = [dataForAnova; trialMeans(:)]; % Fixed reference [cite: 5]
%         groupLabels = [groupLabels; repmat(g, numel(trialMeans), 1)];
%         prefResp = max(prefResp, mean(trialMeans, 'omitnan'));
%     end
% 
%     pValANOVA = anova1(dataForAnova, groupLabels, 'off');
% 
%     tVecB = psth(blankIdx).timeVector;
%     bMask = tVecB >= respWin(1) & tVecB <= respWin(2);
%     blankTrialMeans = squeeze(mean(psth(blankIdx).alignedResponses(i, bMask, :), 2, 'omitnan'));
% 
%     % Paper criteria: ANOVA p<0.05 and peak > blank + 2SD [cite: 2]
%     if pValANOVA < 0.05 && prefResp > (mean(blankTrialMeans, 'omitnan') + 2 * std(blankTrialMeans, 'omitnan'))
%         responsiveROIs = [responsiveROIs; i];
%     end
% end
% 
% %% 3. Generate ROI Pages (Heatmap + PSTHs Only) [cite: 4]
% for p = 1:length(responsiveROIs)
%     roiIdx = responsiveROIs(p);
% 
%     rfMatrix = nan(nEl, nAz); 
%     traceMax = 1e-6; 
% 
%     for k = 1:numel(gridPSTH)
%         mu = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
%         tV = gridPSTH(k).timeVector(:);
%         mu = mu(:) - mean(mu(tV >= baseWin(1) & tV < baseWin(2)), 'omitnan');
%         traceMax = max(traceMax, max(mu));
%         r = find(uEl == gridStim(k,2), 1); c = find(uAz == gridStim(k,1), 1);
%         rfMatrix(r, c) = mean(mu(tV >= respWin(1) & tV <= respWin(2)), 'omitnan');
%     end
% 
%     fig = figure('Color', 'w', 'Position', [50 50 700 600], 'Visible', 'off');
%     ax1 = axes('Position', [0.15 0.15 0.7 0.7]);
% 
%     % Plot Heatmap [cite: 4]
%     imagesc(uAz, uEl, imgaussfilt(rfMatrix, sigma_smooth)); hold on; 
%     colormap(ax1, gray);
% 
%     % Set axes to centered coordinates [cite: 4]
%     set(ax1, 'YDir', 'normal', 'CLim', [0, max(rfMatrix(:)) + 1e-6]);
%     xlim(azLimFull); ylim(elLimFull);
%     colorbar;
% 
%     % Trace Scaling and Centering 
%     vS_base = dEl * 0.4; % 40% 
%     hS = dAz * 0.85;     % 85% 
% 
%     for k = 1:numel(gridPSTH)
%         tr = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
%         tr = tr(:) - mean(tr(tV >= baseWin(1) & tV < baseWin(2)), 'omitnan');
% 
%         % Vertical normalization relative to ROI max [cite: 4]
%         sY = (tr / traceMax * vS_base) + gridStim(k,2);
% 
%         % Horizontal centering [cite: 4]
%         tNorm = (tV - min(tV)) / (max(tV) - min(tV));
%         sX = (tNorm - 0.5) * hS + gridStim(k,1);
% 
%         plot(ax1, sX, sY, 'r', 'LineWidth', 1.2);
%     end
% 
%     title(['ROI ' num2str(roiIdx) ' Heatmap + PSTHs']);
%     xlabel('azimuth (°)'); ylabel('elevation (°)');
% 
%     exportgraphics(fig, pdfPath, 'Append', true);
%     close(fig);
% end
% end


% function findVisuallyResponsiveBoutons_RFMapping(sessionFileInfo, stimName)
% % 1. Setup and Pathing
% iStim = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
% load(sessionFileInfo.stimFiles(iStim).Response, 'response');
% load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
% 
% psth = response.psthData;
% stimVs = vertcat(psth.stimValue);
% nROI = size(psth(1).alignedResponses, 1);
% 
% % Identify Grid vs Blank
% blankIdx = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
% gridMask = stimVs(:,1) ~= 200;
% gridPSTH = psth(gridMask);
% gridStim = stimVs(gridMask, :);
% uAz = sort(unique(gridStim(:,1)), 'ascend');  
% uEl = sort(unique(gridStim(:,2)), 'descend');
% nAz = length(uAz); nEl = length(uEl);
% 
% [dirPath, fileName] = fileparts(sessionFileInfo.stimFiles(iStim).Response);
% pdfPath = fullfile(dirPath, [fileName '_RFMapping_Publication.pdf']);
% if exist(pdfPath, 'file'), delete(pdfPath); end
% 
% % Parameters from Timplalexi et al. 2025 [cite: 460, 483]
% respWin = [0.5 2];  
% baseWin = [-2 0];   
% sigma_smooth = 0.8;
% 
% %% 2. Run Responsiveness Analysis (ANOVA + 2SD Rule)
% responsiveROIs = [];
% for i = 1:nROI
%     % ANOVA across all stimuli + blank [cite: 459, 482]
%     dataForAnova = [];
%     groupLabels = [];
%     prefResp = -inf;
% 
%     for g = 1:numel(psth)
%         tVec = psth(g).timeVector;
%         wMask = tVec >= respWin(1) & tVec <= respWin(2);
%         trialMeans = squeeze(mean(psth(g).alignedResponses(i, wMask, :), 2, 'omitnan'));
%         dataForAnova = [dataForAnova; trialMeans(:)];
%         groupLabels = [groupLabels; repmat(g, numel(trialMeans), 1)];
%         prefResp = max(prefResp, mean(trialMeans, 'omitnan'));
%     end
% 
%     pValANOVA = anova1(dataForAnova, groupLabels, 'off');
% 
%     tVecB = psth(blankIdx).timeVector;
%     bMask = tVecB >= respWin(1) & tVecB <= respWin(2);
%     blankTrialMeans = squeeze(mean(psth(blankIdx).alignedResponses(i, bMask, :), 2, 'omitnan'));
% 
%     % Paper criteria: ANOVA p<0.05 and peak > blank + 2SD [cite: 459, 468, 482]
%     if pValANOVA < 0.05 && prefResp > (mean(blankTrialMeans, 'omitnan') + 2 * std(blankTrialMeans, 'omitnan'))
%         responsiveROIs = [responsiveROIs; i];
%     end
% end
% 
% %% 3. Generate ROI Pages (Figure 1C Style)
% for p = 1:length(responsiveROIs)
%     roiIdx = responsiveROIs(p);
% 
%     rfMatrix = nan(nEl, nAz); 
%     traceMax = 1e-6; 
% 
%     for k = 1:numel(gridPSTH)
%         mu = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
%         tV = gridPSTH(k).timeVector(:);
%         mu = mu(:) - mean(mu(tV >= baseWin(1) & tV < baseWin(2)), 'omitnan');
%         traceMax = max(traceMax, max(mu));
%         r = find(uEl == gridStim(k,2), 1); c = find(uAz == gridStim(k,1), 1);
%         rfMatrix(r, c) = mean(mu(tV >= respWin(1) & tV <= respWin(2)), 'omitnan');
%     end
% 
%     % 2D Gaussian Fit Coordinates [cite: 462, 463, 485, 486]
%     [X, Y] = meshgrid(uAz, uEl);
%     xdata_fit = cat(3, X, Y); 
%     gauss2D = @(x, xdata) x(1) * exp( -( (xdata(:,:,1)-x(2)).^2/(2*x(4)^2) + ...
%                                          (xdata(:,:,2)-x(3)).^2/(2*x(5)^2) ) );
% 
%     [mA, mI] = max(rfMatrix(:));
%     fitP = lsqcurvefit(gauss2D, [mA, X(mI), Y(mI), 15, 15], xdata_fit, rfMatrix, ...
%                        [0 min(uAz) min(uEl) 2 2], [inf max(uAz) max(uEl) 60 60], optimset('Display','off'));
% 
%     fig = figure('Color', 'w', 'Position', [50 50 1100 500], 'Visible', 'off');
%     tlo = tiledlayout(1, 2, 'TileSpacing', 'compact');
% 
%     % LEFT PANEL: Heatmap + Traces 
%     ax1 = nexttile(tlo);
%     imagesc(uAz, uEl, imgaussfilt(rfMatrix, sigma_smooth)); hold on; colormap(ax1, gray);
%     vS = abs(uEl(1)-uEl(2))*0.7; hS = abs(uAz(1)-uAz(2))*0.85;
%     for k = 1:numel(gridPSTH)
%         tr = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
%         tr = tr(:) - mean(tr(tV >= baseWin(1) & tV < baseWin(2)), 'omitnan');
%         % Baseline-anchored scaling
%         sY = (tr/traceMax * vS) + gridStim(k,2);
%         % Horizontal time centering
%         tNorm = (tV - min(tV)) / (max(tV) - min(tV));
%         sX = (tNorm - 0.5) * hS + gridStim(k,1);
%         plot(ax1, sX, sY, 'r', 'LineWidth', 1.5);
%     end
%     set(ax1, 'YDir', 'normal'); title(['ROI ' num2str(roiIdx) ' Heatmap + PSTHs']);
%     xlabel('azimuth (°)'); ylabel('elevation (°)');
% 
%     % RIGHT PANEL: Smoothed Gaussian Fit 
%     ax2 = nexttile(tlo);
%     [Xf, Yf] = meshgrid(linspace(min(uAz), max(uAz), 100), linspace(min(uEl), max(uEl), 100));
%     xfdata_smooth = cat(3, Xf, Yf);
%     Zfit = gauss2D(fitP, xfdata_smooth);
% 
%     imagesc(linspace(min(uAz), max(uAz), 100), linspace(min(uEl), max(uEl), 100), Zfit); hold on;
%     colormap(ax2, parula); shading interp; % Publication yellow-blue style
%     % 1-sigma contour (RF extent) [cite: 86, 89]
%     contour(Xf, Yf, Zfit, [fitP(1)*exp(-0.5) fitP(1)*exp(-0.5)], 'w--', 'LineWidth', 1.5);
%     set(ax2, 'YDir', 'normal'); colorbar;
%     title(sprintf('\\sigma_{az}=%.1f, \\sigma_{el}=%.1f', fitP(4), fitP(5)));
%     xlabel('azimuth (°)');
% 
%     exportgraphics(fig, pdfPath, 'Append', true);
%     close(fig);
% end
% end