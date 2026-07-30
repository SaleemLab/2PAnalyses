function figHandle = plotSingleROI_ConditionStack(sessionFileInfo, response, neuronIdx, signalToUse, saveIndividualPNG, applySmoothing)
% plotSingleROI_ConditionStack: single row layout.
% Column 1 = ONE combined raster (all laps stacked, baseline first, then each
%            condition in turn) with TRUE WHITE GAPS between condition blocks
%            (not separator lines) and a colored gutter marking condition
%            membership. Per-lap normalized to [0,1] purely for display
%            contrast (single shared colorbar) — does not affect the tuning
%            curves.
% Columns 2..N+1 = one tuning curve per condition (mean +/- SEM), computed
%            from the SAME per-lap normalized data as the raster, so both
%            panels share a consistent scale. Each panel shows the baseline
%            as a dashed gray reference for comparison.
%
% saveIndividualPNG (optional, default true): if false, skips saving the
% per-ROI PNG (useful when this function is called in a loop to build a
% single multi-page PDF via plotAllROIs_ConditionStack_PDF).
%
% applySmoothing (optional, default true): if true, each lap's trace is
% smoothed along the position-bin axis with a zero-phase Gaussian filter
% (gausswin + filtfilt) before plotting. NaNs are preserved (masked out
% before filtering, restored after).

if nargin < 4 || isempty(signalToUse), signalToUse = ''; end
if nargin < 5 || isempty(saveIndividualPNG), saveIndividualPNG = true; end
if nargin < 6 || isempty(applySmoothing), applySmoothing = true; end

%% 1. Data Prep
data = response.lapPositionActivity.(signalToUse);
roiActivity = squeeze(data(neuronIdx, :, :));
conds = fieldnames(response.trialIndicesByCondition);

% Filter for conditions that actually have trials to avoid empty columns
validConds = {};
for i = 1:length(conds)
    if ~isempty(response.trialIndicesByCondition.(conds{i}))
        validConds{end+1} = conds{i};
    end
end
numValid = length(validConds);
nBins = size(roiActivity, 2);
xPos = 1:nBins;
tickLocs = [40, 80, 120, 160];
tickLabels = {'40', '80', '120', '160'};
colors = {[0 0 0], [0 0.45 0.74], [0.85 0.33 0.1], [0.2 0.6 0.2], [0.5 0.2 0.5], [0.3 0.7 0.9]};

%% 1b. Optional Smoothing (zero-phase Gaussian filter per lap, NaN-safe)
if applySmoothing
    w = gausswin(10); w = w / sum(w);
    nLapsTotal = size(roiActivity, 1);
    for iL = 1:nLapsTotal
        trace = roiActivity(iL, :);
        if ~all(isnan(trace))
            nanMask = isnan(trace);
            trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace);
            smoothed(nanMask) = NaN;
            roiActivity(iL, :) = smoothed;
        end
    end
end

%% 2. Per-Lap [0,1] Normalization (used for BOTH the raster and tuning curves)
% Each lap normalized to its own min/max along position bins. Using the
% SAME normalized matrix for both panels keeps them on a consistent scale.
roiActivityNorm = normalize(roiActivity, 2, 'range');

%% 2b. Y-Axis Limit for Tuning Curves (Based on NORMALIZED Means)
maxMean = 0;
for iC = 1:numValid
    mu = mean(roiActivityNorm(response.trialIndicesByCondition.(validConds{iC}), :), 1, 'omitnan');
    maxMean = max(maxMean, max(mu));
end
yLimit = maxMean * 1;
if yLimit <= 0 || isnan(yLimit), yLimit = 1; end

% Reference Baseline
baseIdx = find(contains(lower(validConds), 'baseline'), 1);
if isempty(baseIdx), baseIdx = 1; end
baseMu = mean(roiActivityNorm(response.trialIndicesByCondition.(validConds{baseIdx}), :), 1, 'omitnan');

%% 3. Build the combined stacked raster: baseline first, then remaining
% conditions in their original order. A block of NaN rows is inserted
% between each condition so the separation is a genuine WHITE GAP (via
% AlphaData transparency) rather than a drawn line.
gapRows = max(5, round(0.12 * sum(cellfun(@(c) numel(response.trialIndicesByCondition.(c)), validConds))));
stackOrder = [baseIdx, setdiff(1:numValid, baseIdx, 'stable')];

plotMatrix = [];
blockStart = zeros(1, numValid);
blockEnd = zeros(1, numValid);
rowPtr = 1;
for k = 1:numValid
    iC = stackOrder(k);
    theseIdx = response.trialIndicesByCondition.(validConds{iC});
    theseIdx = theseIdx(:)';
    n = numel(theseIdx);

    blockData = roiActivityNorm(theseIdx, :);

    blockStart(iC) = rowPtr;
    blockEnd(iC) = rowPtr + n - 1;
    rowPtr = rowPtr + n;

    plotMatrix = [plotMatrix; blockData]; %#ok<AGROW>

    if k < numValid
        plotMatrix = [plotMatrix; nan(gapRows, nBins)];
        rowPtr = rowPtr + gapRows;
    end
end
totalRows = rowPtr - 1;

%% 4. Figure Setup (single row: raster column + N tuning-curve columns)
nCols = numValid + 1;
figWidth = max(1100, 280*nCols + 200);
fig = figure('Color', 'w', 'Position', [50 50 figWidth 250]);
t = tiledlayout(1, nCols, 'TileSpacing', 'compact', 'Padding', 'loose');
title(t, sprintf('%s | ROI %d | %s', sessionFileInfo.session_name, neuronIdx, signalToUse), ...
    'FontSize', 14, 'FontWeight', 'bold');

%% 5. COLUMN 1: combined raster with white gaps + colored gutter
axTop = nexttile(t, 1); hold(axTop, 'on');

himg = imagesc(axTop, xPos, 1:totalRows, plotMatrix);
set(himg, 'AlphaData', ~isnan(plotMatrix)); % NaN gap rows render as transparent -> white
set(axTop, 'Color', 'w'); % axes background shows through the transparent gaps
colormap(axTop, flipud(gray));
set(axTop, 'CLim', [0 1]);

% Colored gutter: one dot per lap, just left of the raster, colored by
% which condition that lap belongs to (same colors used for the tuning
% curve columns), visually linking the raster to its matching column.
gutterX = -0.06 * nBins;
for iC = 1:numValid
    currColor = colors{mod(iC-1, length(colors))+1};
    rows = blockStart(iC):blockEnd(iC);
    s = scatter(axTop, repmat(gutterX, size(rows)), rows, 2, currColor, 'filled');
    s.Clipping = 'off';
end

for p = tickLocs
    xline(axTop, p, 'k--', 'Alpha', 0.3);
end

set(axTop, 'YDir', 'reverse', 'Box', 'off', 'XTick', tickLocs, 'XTickLabel', tickLabels, ...
    'XLim', [gutterX - 0.03*nBins, nBins], 'YLim', [0.5, totalRows + 0.5]);
ylabel(axTop, 'Laps', 'FontSize', 11, 'FontWeight', 'bold');
xlabel(axTop, 'Position (cm)', 'FontWeight', 'bold');
title(axTop, 'Lap Activity', 'FontSize', 11, 'FontWeight', 'bold');

cb = colorbar(axTop);
cb.Label.String = 'Norm. Activity';

defaultAxesProperties(axTop, true);
offsetAxes(axTop)

%% 6. COLUMNS 2..N+1: one raw-units tuning curve per condition
for iC = 1:numValid
    currCond = validConds{iC};
    trialIDs = response.trialIndicesByCondition.(currCond);
    trialIDs = trialIDs(:)';
    nLapsInBlock = length(trialIDs);
    currColor = colors{mod(iC-1, length(colors))+1};

    axR = nexttile(t, 1 + iC); hold(axR, 'on');
    if iC ~= baseIdx
        % Gray dashed baseline reference on OTHER conditions' panels only —
        % on baseline's own panel it would just overlap its own solid line.
        plot(axR, xPos, baseMu, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
    end

    mu = mean(roiActivityNorm(trialIDs, :), 1, 'omitnan');
%     sem = std(roiActivityNorm(trialIDs, :), 0, 1, 'omitnan') ./ sqrt(nLapsInBlock);
%     fill(axR, [xPos, fliplr(xPos)], [mu+sem, fliplr(mu - sem)], currColor, ...
%         'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(axR, xPos, mu, 'Color', currColor, 'LineWidth', 2);

    for p = tickLocs
        xline(axR, p, 'k--', 'Alpha', 0.3);
    end

    title(axR, sprintf('%s (n=%d)', strrep(currCond, '_', ' '), nLapsInBlock), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', currColor);

    set(axR, 'Box', 'off', 'XTick', tickLocs, 'XTickLabel', tickLabels, 'ylim', [0 yLimit]);
    xlabel(axR, 'Position (cm)', 'FontWeight', 'bold');

    if iC == 1
        ylabel(axR, 'Norm. Activity \pm SEM', 'FontSize', 10, 'FontWeight', 'bold');
    end

    defaultAxesProperties(axR, true);
    offsetAxes(axR)
end

%% 7. Save (optional per-ROI PNG)
if saveIndividualPNG
    saveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'ROI_Summaries');
    if ~exist(saveDir, 'dir'), mkdir(saveDir); end
    saveas(fig, fullfile(saveDir, sprintf('ROI_%d_Summary_Stacked.png', neuronIdx)));
end

if nargout > 0, figHandle = fig; end
end
% function figHandle = plotSingleROI_ConditionStack(sessionFileInfo, response, neuronIdx, signalToUse)
% % plotSingleROI_ConditionStack: Robust vertical stack with sequential lap counts.
% % FIXED: Handles explicitly passed to all graphics functions to ensure 
% % defaultAxesProperties and offsetAxes are accurately targeted and executed.
% 
% if nargin < 4, signalToUse = 'dFFNeuropilCorrected'; end
% 
% %% 1. Data Prep
% data = response.lapPositionActivity.(signalToUse);
% roiActivity = squeeze(data(neuronIdx, :, :));
% conds = fieldnames(response.trialIndicesByCondition);
% 
% % Filter for conditions that actually have trials to avoid empty rows
% validConds = {};
% for i = 1:length(conds)
%     if ~isempty(response.trialIndicesByCondition.(conds{i}))
%         validConds{end+1} = conds{i};
%     end
% end
% numValid = length(validConds);
% nBins = size(roiActivity, 2);
% xPos = 1:nBins;
% tickLocs = [40, 80, 120, 160]; 
% tickLabels = {'40', '80', '120', '160'};
% colors = {[0 0 0], [0 0.45 0.74], [0.85 0.33 0.1], [0.2 0.6 0.2], [0.5 0.2 0.5], [0.3 0.7 0.9]};
% 
% %% 2. Optimized Y-Axis Limit (Based on Means)
% maxMean = 0;
% for iC = 1:numValid
%     mu = mean(roiActivity(response.trialIndicesByCondition.(validConds{iC}), :), 1, 'omitnan');
%     maxMean = max(maxMean, max(mu));
% end
% yLimit = maxMean * 1; 
% if yLimit <= 0 || isnan(yLimit), yLimit = 1; end
% 
% %% 3. Figure Setup
% fig = figure('Color', 'w', 'Position', [100 50 1000 950]);
% t = tiledlayout(numValid, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
% title(t, sprintf('%s | ROI %d | %s', sessionFileInfo.session_name, neuronIdx, signalToUse), ...
%     'FontSize', 14, 'FontWeight', 'bold');
% 
% % Reference Baseline
% baseIdx = find(contains(lower(validConds), 'baseline'), 1);
% if isempty(baseIdx), baseIdx = 1; end
% baseMu = mean(roiActivity(response.trialIndicesByCondition.(validConds{baseIdx}), :), 1, 'omitnan');
% 
% %% 4. Vertical Stack Loop
% for iC = 1:numValid
%     currCond = validConds{iC};
%     trialIDs = response.trialIndicesByCondition.(currCond);
%     nLapsInBlock = length(trialIDs);
%     currColor = colors{mod(iC-1, length(colors))+1};
%     
%     % --- LEFT COLUMN: RASTER ---
%     axL = nexttile(t); hold(axL, 'on'); % Explicit tile destination assignment
%     
%     imagesc(axL, xPos, 1:nLapsInBlock, roiActivity(trialIDs, :));
%     colormap(axL, flipud(gray));
%     
%     % FIXED: Switched from gca to explicit axL handle targeting
%     set(axL, 'YDir', 'reverse', 'Box', 'off', 'XTick', tickLocs, ...
%              'YTick', [1, nLapsInBlock], 'YTickLabel', {'1', num2str(nLapsInBlock)});
%     
%     ylabel(axL, sprintf('%s\n(Laps)', strrep(currCond, '_', ' ')), ...
%         'FontSize', 11, 'FontWeight', 'bold', 'Color', currColor);
%     
%     if iC == numValid
%         set(axL, 'XTickLabel', tickLabels);
%         xlabel(axL, 'Position (cm)', 'FontWeight', 'bold');
%     else
%         set(axL, 'XTickLabel', []);
%     end
%     
%     % Call properties on the handle directly
%     defaultAxesProperties(axL, true);
%     offsetAxes(axL)
%     
%     % --- RIGHT COLUMN: TUNING CURVE ---
%     axR = nexttile(t); hold(axR, 'on'); % Explicit tile destination assignment
%     plot(axR, xPos, baseMu, 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'LineStyle', '--');
%     
%     mu = mean(roiActivity(trialIDs, :), 1, 'omitnan');
%     sem = std(roiActivity(trialIDs, :), 0, 1, 'omitnan') ./ sqrt(nLapsInBlock);
%     fill(axR, [xPos, fliplr(xPos)], [mu+sem, fliplr(mu - sem)], currColor, ...
%         'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     plot(axR, xPos, mu, 'Color', currColor, 'LineWidth', 3);
%     
%     for p = tickLocs
%         xline(axR, p, 'k--', 'Alpha', 0.3); 
%     end
%     
%     % FIXED: Switched from gca to explicit axR handle targeting
%     set(axR, 'Box', 'off', 'XTick', tickLocs, 'ylim', [-1 yLimit]);
%     
%     if iC == numValid
%         set(axR, 'XTickLabel', tickLabels);
%         xlabel(axR, 'Position (cm)', 'FontWeight', 'bold');
%     else
%         set(axR, 'XTickLabel', []);
%     end
%     
%     if iC == 1
%         %title(axR, 'Mean \Delta F/F \pm SEM', 'FontSize', 10);
%         title(axR, 'Mean Activity \pm SEM', 'FontSize', 10);
%     end
%     
%     % Call properties on the handle directly
%     defaultAxesProperties(axR, true);
%     offsetAxes(axR)
% end
% 
% %% 5. Save
% saveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'ROI_Summaries');
% if ~exist(saveDir, 'dir'), mkdir(saveDir); end
% saveas(fig, fullfile(saveDir, sprintf('ROI_%d_Summary_Stacked.png', neuronIdx)));
% 
% if nargout > 0, figHandle = fig; end
% end