function figHandle = plotSpeedPositionActivity_ForROI(response, targetROI, applySmoothing, smoothSigma)
% Plots panel b (speed-stratified position tuning curves, low vs high only)
% and panel c (continuous 2D speed x position heatmap) for a single ROI.
%
% Requires BOTH to have been run first:
%   - getPositionSpeedMatrix               -> panel c heatmap (speedPositionActivity.continuous)
%   - getLowHighSpeedPositionMatrix        -> panel b lines   (speedPositionActivity.lowHigh)
%
% Inputs:
%   response       : loaded response struct
%   targetROI      : ROI index to plot
%   applySmoothing : true (default)
%   smoothSigma    : [2.0, 2.5] (default) [Speed, Position]

if nargin < 3, applySmoothing = true; end
if nargin < 4, smoothSigma    = [2.0, 2.5]; end

%% Check both pipelines have been run
if ~isfield(response, 'speedPositionActivity')
    error('speedPositionActivity not found. Run getPositionSpeedMatrix and getLowHighSpeedPositionMatrix first.');
end
if ~isfield(response.speedPositionActivity, 'continuous')
    error('continuous matrix not found. Run getPositionSpeedMatrix first.');
end
if ~isfield(response.speedPositionActivity, 'lowHigh')
    error('lowHigh not found. Run getLowHighSpeedPositionMatrix first.');
end

%% Panel c: continuous heatmap matrix
tuningSurface = response.speedPositionActivity.continuous.matrix(:, :, targetROI);
speedEdges    = response.speedPositionActivity.continuous.speedEdges;
[~, numPosBins] = size(tuningSurface);

%% Panel b: low / high lines from getLowHighSpeedPositionMatrix
lowLine  = response.speedPositionActivity.lowHigh.matrixLow(targetROI,  :)'; % [numPosBins x 1]
highLine = response.speedPositionActivity.lowHigh.matrixHigh(targetROI, :)';

% Per-position median threshold line
medianThreshLine    = response.speedPositionActivity.lowHigh.medianThreshLine;
medianThreshSummary = median(medianThreshLine, 'omitnan');

%% Smoothing
if applySmoothing
    % Panel c — 2D gaussian smooth
    mask       = ~isnan(tuningSurface);
    dataZeroed = tuningSurface;
    dataZeroed(~mask) = 0;
    blurredData   = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
    blurredMask   = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
    tuningSurface = blurredData ./ blurredMask;

    % Panel b — 1D position smoothing
    lowLine  = smoothdata(lowLine,  'gaussian', smoothSigma(2) * 3);
    highLine = smoothdata(highLine, 'gaussian', smoothSigma(2) * 3);
end

tuningSurface(isnan(tuningSurface)) = 0;

%% Colours and layout
colorLow  = [0.00, 0.60, 0.80];  % blue
colorHigh = [0.80, 0.00, 0.60];  % magenta

landmarkCentres = [40, 80, 120, 160];
x_pos      = 1:numPosBins;
leftMargin = 0.12;
plotWidth  = 0.62;

minY = speedEdges(1);
maxY = speedEdges(end);

%% Figure
figHandle = figure('Name', sprintf('ROI %d - Speed Position Activity', targetROI), ...
    'Position', [100 100 720 600], 'Color', 'w', 'Visible', 'off');

%% Panel b — low vs high tuning curves
ax1 = subplot('Position', [leftMargin, 0.68, plotWidth, 0.24]);
hold on;

for c = landmarkCentres
    if c <= numPosBins
        xline(c, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
    end
end

hLow  = plot(ax1, x_pos, lowLine,  'Color', colorLow,  'LineWidth', 2);
hHigh = plot(ax1, x_pos, highLine, 'Color', colorHigh, 'LineWidth', 2);

legend([hLow, hHigh], ...
    {sprintf('Low  (<%.0f cm/s)',  medianThreshSummary), ...
     sprintf('High (>%.0f cm/s)', medianThreshSummary)}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 8);

ylabel('\DeltaF/F (NeuC)');
set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
title(sprintf('ROI %d — %s', targetROI, response.stimName), 'FontWeight', 'bold');

if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax1, 0); end
if exist('offsetAxes',            'file') == 2, offsetAxes(ax1); end

%% Panel c — continuous 2D heatmap
ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.50]);

displayMatrix = padarray(tuningSurface, [1 1], 0, 'post');
xEdges = 0.5:(numPosBins + 0.5);

p = pcolor(xEdges, speedEdges, displayMatrix);
set(p, 'EdgeColor', 'none');
set(ax2, 'YDir', 'normal');
colormap(ax2, parula);

activeData = tuningSurface(tuningSurface > 0);
maxVal = 1;
if ~isempty(activeData), maxVal = prctile(activeData, 98.5); end
set(ax2, 'CLim', [0, maxVal]);

cb = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.025, 0.50]);
cb.Label.String = '\DeltaF/F [NeuC]';

hold on;

for c = landmarkCentres
    if c <= numPosBins
        xline(c, '--w', 'Alpha', 0.3, 'LineWidth', 1.2);
    end
end

% Per-position median threshold line
plot(ax2, 1:numPosBins, medianThreshLine, '--', 'Color', [0.9 0.9 0.9], 'LineWidth', 1.5);

xticks([1 40 80 120 160 200]);
xticklabels({'1', '40', '80', '120', '160', '200'});
xlim([0.5, numPosBins + 0.5]);
ylim([minY, maxY]);

ytickVals = unique([round(minY), round(medianThreshSummary), round(maxY)]);
yticks(ytickVals);
yticklabels(arrayfun(@num2str, ytickVals, 'UniformOutput', false));

xlabel('Position (cm)');
ylabel('Running speed (cm/s)');
set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax2, 0); end

%% Stratification colour bar — two zones only (low / high)
axStrat = axes('Position', [leftMargin + plotWidth + 0.11, 0.12, 0.02, 0.50]);
hold on;

patch([0 1 1 0], [minY                minY                medianThreshSummary medianThreshSummary], colorLow,  'EdgeColor', 'none', 'FaceAlpha', 0.5);
patch([0 1 1 0], [medianThreshSummary medianThreshSummary maxY                maxY],                colorHigh, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

set(axStrat, 'YLim', [minY, maxY], 'XLim', [0, 1]);
yTickPlacements = [mean([minY, medianThreshSummary]), mean([medianThreshSummary, maxY])];
set(axStrat, 'YTick', yTickPlacements, ...
             'YTickLabel', {'Low', 'High'}, ...
             'YAxisLocation', 'right', ...
             'XTick', [], 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
ylabel(axStrat, 'Speed zones', 'FontSize', 10);

%% Link x axes
linkaxes([ax1, ax2], 'x');
set(figHandle, 'Visible', 'on');
end

% function figHandle = plotSpeedPositionActivity_ForROI(response, targetROI, applySmoothing, smoothSigma)
% % Plots panel b (speed-stratified position tuning curves, low vs high only)
% % and panel c (continuous 2D speed x position heatmap) for a single ROI.
% %
% % Changes from previous version:
% %   - Baseline trials filter applied to speedData before threshold computation
% %   - Occupancy-based per-position tertile thresholds (not fixed 15/27 cm/s)
% %   - Low vs High only (median split) — medium dropped
% %   - Tuning curves in panel b recomputed using same occupancy thresholds
% %   - Stratification colorbar updated to two zones (low/high)
% %   - Heatmap yticks reflect actual occupancy threshold values
% %
% % Requires getPositionSpeedMatrix to have been run first (panel c heatmap).
% %
% % Inputs:
% %   response       : loaded response struct
% %   targetROI      : ROI index to plot
% %   applySmoothing : true (default)
% %   smoothSigma    : [1.0, 2.0] (default) [Speed, Position]
% 
% if nargin < 3, applySmoothing = true; end
% if nargin < 4, smoothSigma    = [2.0, 2.5]; end
% 
% %% Check pipeline has been run
% if ~isfield(response, 'speedPositionActivity') % || ~isfield(response.speedPositionActivity, 'continuous')
%     error('speedPositionActivity.continuous not found. Run getPositionSpeedMatrix first.');
% end
% 
% %% Load raw data — baseline trials only
% rawROIData     = response.lapPositionActivity.dFFNeuropilCorrected; % [numROIs x numTrials x numPosBins]
% allSpeedData   = response.lapPositionRunningSpeed;                  % [numTrials x numPosBins]
% baselineTrials = response.trialIndicesByCondition.Baseline;
% 
% rawROIData = rawROIData(:, baselineTrials, :);
% speedData  = allSpeedData(baselineTrials, :);
% 
% [~, ~, numPosBins] = size(rawROIData);
% 
% %% Panel c: continuous heatmap matrix
% tuningSurface = response.speedPositionActivity.matrix(:, :, targetROI);
% speedCenters  = response.speedPositionActivity.speedBinCenters;
% speedEdges    = response.speedPositionActivity.speedEdges;
% 
% %% Compute per-position occupancy-based thresholds (median split → low / high)
% medianThreshLine = nan(1, numPosBins);
% 
% for b = 1:numPosBins
%     binSpeeds   = speedData(:, b);
%     validSpeeds = binSpeeds(binSpeeds > 1 & ~isnan(binSpeeds));
%     if length(validSpeeds) >= 6
%         medianThreshLine(b) = median(validSpeeds);
%     end
% end
% 
% %% Recompute low / high tuning curves using same occupancy thresholds
% lowLine  = nan(numPosBins, 1);
% highLine = nan(numPosBins, 1);
% 
% for b = 1:numPosBins
%     if isnan(medianThreshLine(b)), continue; end
%     binSpeeds = speedData(:, b);
%     roiData   = squeeze(rawROIData(targetROI, :, b));
% 
%     lowIdx  = binSpeeds > 1 & binSpeeds <= medianThreshLine(b) & ~isnan(binSpeeds);
%     highIdx = binSpeeds > medianThreshLine(b) & ~isnan(binSpeeds);
% 
%     if sum(lowIdx)  >= 2, lowLine(b)  = mean(roiData(lowIdx),  'omitnan'); end
%     if sum(highIdx) >= 2, highLine(b) = mean(roiData(highIdx), 'omitnan'); end
% end
% 
% %% Smoothing
% if applySmoothing
%     % Panel c — 2D gaussian smooth
%     mask       = ~isnan(tuningSurface);
%     dataZeroed = tuningSurface;
%     dataZeroed(~mask) = 0;
%     blurredData   = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
%     blurredMask   = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
%     tuningSurface = blurredData ./ blurredMask;
% 
%     % Panel b — 1D position smoothing
%     lowLine  = smoothdata(lowLine,  'gaussian', smoothSigma(2) * 3);
%     highLine = smoothdata(highLine, 'gaussian', smoothSigma(2) * 3);
% end
% 
% tuningSurface(isnan(tuningSurface)) = 0;
% 
% %% Compute representative threshold values for legend/yticks
% % Use median of the per-position threshold line as a single summary value
% medianThreshSummary = median(medianThreshLine, 'omitnan');
% 
% %% Colours and layout
% colorLow  = [0.00, 0.60, 0.80];  % blue
% colorHigh = [0.80, 0.00, 0.60];  % magenta
% 
% landmarkCentres = [40, 80, 120, 160];
% x_pos      = 1:numPosBins;
% leftMargin = 0.12;
% plotWidth  = 0.62;
% 
% minY = speedEdges(1);
% maxY = speedEdges(end);
% 
% %% Figure
% figHandle = figure('Name', sprintf('ROI %d - Speed Position Activity', targetROI), ...
%     'Position', [100 100 720 600], 'Color', 'w', 'Visible', 'off');
% 
% %% Panel b — low vs high tuning curves
% ax1 = subplot('Position', [leftMargin, 0.68, plotWidth, 0.24]);
% hold on;
% 
% for c = landmarkCentres
%     if c <= numPosBins
%         xline(c, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
%     end
% end
% 
% hLow  = plot(ax1, x_pos, lowLine,  'Color', colorLow,  'LineWidth', 2);
% hHigh = plot(ax1, x_pos, highLine, 'Color', colorHigh, 'LineWidth', 2);
% 
% legend([hLow, hHigh], ...
%     {sprintf('Low  (<%.0f cm/s)',  medianThreshSummary), ...
%      sprintf('High (>%.0f cm/s)', medianThreshSummary)}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
% 
% ylabel('\DeltaF/F (NeuC)');
% set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
% title(sprintf('ROI %d — %s', targetROI, response.stimName), 'FontWeight', 'bold');
% 
% % if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax1, 0); end
% % if exist('offsetAxes',            'file') == 2, offsetAxes(ax1); end
% 
% %% Panel c — continuous 2D heatmap
% ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.50]);
% 
% displayMatrix = padarray(tuningSurface, [1 1], 0, 'post');
% xEdges = 0.5:(numPosBins + 0.5);
% 
% p = pcolor(xEdges, speedEdges, displayMatrix);
% set(p, 'EdgeColor', 'none');
% set(ax2, 'YDir', 'normal');
% colormap(ax2, parula);
% 
% activeData = tuningSurface(tuningSurface > 0);
% maxVal = 1;
% if ~isempty(activeData), maxVal = prctile(activeData, 98.5); end
% set(ax2, 'CLim', [0, maxVal]);
% 
% cb = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.025, 0.50]);
% cb.Label.String = '\DeltaF/F [NeuC]';
% 
% hold on;
% 
% % Landmark lines
% for c = landmarkCentres
%     if c <= numPosBins
%         xline(c, '--w', 'Alpha', 0.3, 'LineWidth', 1.2);
%     end
% end
% 
% % Occupancy-based threshold line (position-varying median)
% plot(ax2, x_pos, medianThreshLine, '--', 'Color', [0.9 0.9 0.9], 'LineWidth', 1.5);
% 
% xticks([1 40 80 120 160 200]);
% xticklabels({'1', '40', '80', '120', '160', '200'});
% xlim([0.5, numPosBins + 0.5]);
% ylim([minY, maxY]);
% 
% % yticks: just show min, summary threshold, max
% ytickVals = unique([round(minY), round(medianThreshSummary), round(maxY)]);
% yticks(ytickVals);
% yticklabels(arrayfun(@num2str, ytickVals, 'UniformOutput', false));
% 
% xlabel('Position (cm)');
% ylabel('Running speed (cm/s)');
% set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
% % if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax2, 0); end
% 
% %% Stratification colour bar — two zones only (low / high)
% axStrat = axes('Position', [leftMargin + plotWidth + 0.11, 0.12, 0.02, 0.50]);
% hold on;
% 
% patch([0 1 1 0], [minY                 minY                 medianThreshSummary  medianThreshSummary], colorLow,  'EdgeColor', 'none', 'FaceAlpha', 0.5);
% patch([0 1 1 0], [medianThreshSummary  medianThreshSummary  maxY                 maxY],                colorHigh, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% 
% set(axStrat, 'YLim', [minY, maxY], 'XLim', [0, 1]);
% yTickPlacements = [mean([minY, medianThreshSummary]), mean([medianThreshSummary, maxY])];
% set(axStrat, 'YTick', yTickPlacements, ...
%              'YTickLabel', {'Low', 'High'}, ...
%              'YAxisLocation', 'right', ...
%              'XTick', [], 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
% ylabel(axStrat, 'Speed zones', 'FontSize', 10);
% 
% %% Link x axes
% linkaxes([ax1, ax2], 'x');
% set(figHandle, 'Visible', 'on');
% end

% function figHandle = plotSpeedPositionActivity_ForROI(response, targetROI, applySmoothing, smoothSigma)
% % Plots panel b (speed-stratified position tuning curves) and
% % panel c (continuous 2D speed x position heatmap) for a single ROI.
% %
% % Requires BOTH functions to have been run first:
% %   - getPositionSpeedMatrix               -> panel c heatmap
% %   - getLowMedHigh_SpeedPosActivityMatrix -> panel b lines
% %
% % Inputs:
% %   response       : loaded response struct containing both outputs
% %   targetROI      : ROI index to plot
% %   applySmoothing : true (default) - smooths both panel b and panel c
% %   smoothSigma    : [1.0, 5.0] (default) [Speed, Position]
% 
% if nargin < 3, applySmoothing = true; end
% if nargin < 4, smoothSigma    = [1.0, 2.0]; end
% 
% %% Check both pipelines have been run
% if ~isfield(response, 'speedPositionActivity')
%     error('speedPositionActivity not found. Run getPositionSpeedMatrix and getLowMedHigh_SpeedPosActivityMatrix first.');
% end
% if ~isfield(response.speedPositionActivity, 'continuous')
%     error('Continuous matrix not found. Run getPositionSpeedMatrix first.');
% end
% if ~isfield(response.speedPositionActivity, 'stratified')
%     error('Stratified ranges not found. Run getLowMedHigh_SpeedPosActivityMatrix first.');
% end
% 
% %% Panel c: continuous matrix
% tuningSurface   = response.speedPositionActivity.continuous.matrix(:, :, targetROI);
% speedCenters    = response.speedPositionActivity.continuous.speedBinCenters;
% speedEdges      = response.speedPositionActivity.continuous.speedEdges;
% [~, numPosBins] = size(tuningSurface);
% 
% 
% 
% %% Compute position-varying thresholds from occupancy
% speedData = response.lapPositionRunningSpeed; % 
% 
% 
% lowThreshLine  = nan(1, numPosBins);
% highThreshLine = nan(1, numPosBins);
% 
% for b = 1:numPosBins
%     binSpeeds = speedData(:, b);
%     validSpeeds = binSpeeds(binSpeeds > 1 & ~isnan(binSpeeds));
%     if length(validSpeeds) >= 3
%         thresholds = quantile(validSpeeds, [1/3 2/3]);
%         lowThreshLine(b)  = thresholds(1);
%         highThreshLine(b) = thresholds(2);
%     end
% end
% %% Panel b: low/med/high lines from getLowMedHigh output
% lowLine  = squeeze(response.speedPositionActivity.stratified.speedRanges.low(:,  targetROI));
% medLine  = squeeze(response.speedPositionActivity.stratified.speedRanges.med(:,  targetROI));
% highLine = squeeze(response.speedPositionActivity.stratified.speedRanges.high(:, targetROI));
% 
% fixedEdges = response.speedPositionActivity.stratified.speedEdges;
% lowThresh  = fixedEdges(2);  % 15 cm/s
% highThresh = fixedEdges(3);  % 27 cm/s
% 
% %% Smoothing — applied to both panel b and panel c consistently
% if applySmoothing
%     % panel c — 2D gaussian smooth [Speed x Position]
%     mask       = ~isnan(tuningSurface);
%     dataZeroed = tuningSurface;
%     dataZeroed(~mask) = 0;
%     blurredData   = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
%     blurredMask   = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
%     tuningSurface = blurredData ./ blurredMask;
% 
%     % panel b — same position sigma as panel c
%     lowLine  = smoothdata(lowLine,  'gaussian', smoothSigma(2) * 3);
%     medLine  = smoothdata(medLine,  'gaussian', smoothSigma(2) * 3);
%     highLine = smoothdata(highLine, 'gaussian', smoothSigma(2) * 3);
% end
% 
% % set any remaining nans to 0 for display
% tuningSurface(isnan(tuningSurface)) = 0;
% 
% %% Colours and layout
% colorLow  = [0.00, 0.60, 0.80];  % blue
% colorMed  = [0.20, 0.20, 0.20];  % dark grey
% colorHigh = [0.80, 0.00, 0.60];  % magenta
% 
% landmarkCentres = [40, 80, 120, 160];
% x_pos      = 1:numPosBins;
% leftMargin = 0.12;
% plotWidth  = 0.62;
% 
% minY = speedEdges(1);
% maxY = speedEdges(end);
% 
% %% Figure
% figHandle = figure('Name', sprintf('ROI %d - Speed Position Activity', targetROI), ...
%     'Position', [100 100 720 600], 'Color', 'w', 'Visible', 'off');
% 
% %% Panel b — speed-stratified tuning curves
% ax1 = subplot('Position', [leftMargin, 0.68, plotWidth, 0.24]);
% hold on;
% 
% for c = landmarkCentres
%     if c <= numPosBins
%         xline(c, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
%     end
% end
% 
% hLow  = plot(ax1, x_pos, lowLine,  'Color', colorLow,  'LineWidth', 2);
% hMed  = plot(ax1, x_pos, medLine,  'Color', colorMed,  'LineWidth', 2);
% hHigh = plot(ax1, x_pos, highLine, 'Color', colorHigh, 'LineWidth', 2);
% 
% legend([hLow, hMed, hHigh], ...
%     {sprintf('Low  (<%.0f cm/s)',      lowThresh), ...
%      sprintf('Med  (%.0f-%.0f cm/s)',  lowThresh, highThresh), ...
%      sprintf('High (>%.0f cm/s)',      highThresh)}, ...
%     'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
% 
% ylabel('\DeltaF/F (NeuC)');
% set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
% title(sprintf('ROI %d — %s', targetROI, response.stimName), 'FontWeight', 'bold');
% 
% if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax1, 0); end
% if exist('offsetAxes',            'file') == 2, offsetAxes(ax1); end
% 
% %% Panel c — continuous 2D heatmap
% ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.50]);
% 
% displayMatrix = padarray(tuningSurface, [1 1], 0, 'post');
% xEdges = 0.5:(numPosBins + 0.5);
% 
% p = pcolor(xEdges, speedEdges, displayMatrix);
% set(p, 'EdgeColor', 'none');
% set(ax2, 'YDir', 'normal');
% colormap(ax2, parula);
% 
% activeData = tuningSurface(tuningSurface > 0);
% maxVal = 1;
% if ~isempty(activeData), maxVal = prctile(activeData, 98.5); end
% set(ax2, 'CLim', [0, maxVal]);
% 
% cb = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.025, 0.50]);
% cb.Label.String = '\DeltaF/F [NeuC]';
% 
% hold on;
% for c = landmarkCentres
%     if c <= numPosBins
%         xline(c, '--w', 'Alpha', 0.3, 'LineWidth', 1.2);
%     end
% end
% 
% % yline(ax2, lowThresh,  '--', 'Color', colorLow,  'LineWidth', 1.5, 'Alpha', 0.8);
% % yline(ax2, highThresh, '--', 'Color', colorHigh, 'LineWidth', 1.5, 'Alpha', 0.8);
% 
% % replace yline(ax2, lowThresh ...) and yline(ax2, highThresh ...) with:
% plot(ax2, x_pos, lowThreshLine,  '--', 'Color', colorLow,  'LineWidth', 1.5);
% plot(ax2, x_pos, highThreshLine, '--', 'Color', colorHigh, 'LineWidth', 1.5);
% 
% xticks([1 40 80 120 160 200]);
% xticklabels({'1', '40', '80', '120', '160', '200'});
% xlim([0.5, numPosBins + 0.5]);
% ylim([minY, maxY]);
% yticks(unique([1, lowThresh, highThresh, round(maxY)]));
% yticklabels(arrayfun(@num2str, unique([1, lowThresh, highThresh, round(maxY)]), 'UniformOutput', false));
% 
% xlabel('Position (cm)');
% ylabel('Running speed (cm/s)');
% set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
% if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax2, 0); end
% 
% %% Stratification colour bar (right side)
% axStrat = axes('Position', [leftMargin + plotWidth + 0.11, 0.12, 0.02, 0.50]);
% hold on;
% 
% patch([0 1 1 0], [minY       minY       lowThresh  lowThresh],  colorLow,  'EdgeColor', 'none', 'FaceAlpha', 0.5);
% patch([0 1 1 0], [lowThresh  lowThresh  highThresh highThresh], colorMed,  'EdgeColor', 'none', 'FaceAlpha', 0.3);
% patch([0 1 1 0], [highThresh highThresh maxY       maxY],       colorHigh, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% 
% set(axStrat, 'YLim', [minY, maxY], 'XLim', [0, 1]);
% yTickPlacements = [mean([minY, lowThresh]), mean([lowThresh, highThresh]), mean([highThresh, maxY])];
% set(axStrat, 'YTick', yTickPlacements, ...
%              'YTickLabel', {'Low', 'Med', 'High'}, ...
%              'YAxisLocation', 'right', ...
%              'XTick', [], 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
% ylabel(axStrat, 'Speed zones', 'FontSize', 10);
% 
% %% Link x axes
% linkaxes([ax1, ax2], 'x');
% set(figHandle, 'Visible', 'on');
% end

% function figHandle = plotSpeedPositionActivity_ForROI(response, targetROI, applySmoothing, smoothSigma)
%     
%     if nargin < 4, applySmoothing = true; end
%     if nargin < 5, smoothSigma = [2.1, 2.5]; end
%     if nargin < 6, sessionName = 'session'; end
%     
%     tuningSurface = response.speedPositionActivity.matrix(:, :, targetROI);
%     speedCenters  = response.speedPositionActivity.speedBinCenters;
%     speedEdges    = response.speedPositionActivity.speedEdges;
%     [numSpeedBins, numPosBins] = size(tuningSurface);
%     
%     if applySmoothing
%         mask = ~isnan(tuningSurface);
%         dataZeroed = tuningSurface;
%         dataZeroed(~mask) = 0; 
%         blurredData = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
%         blurredMask = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
%         tuningSurface = blurredData ./ blurredMask;
%         tuningSurface(isnan(tuningSurface)) = 0;
%     else
%         tuningSurface(isnan(tuningSurface)) = 0;
%     end
%     
%     binsPerTier = 3; 
%     
%     lowRows  = 1:binsPerTier;
%     medRows  = (binsPerTier + 1):(2 * binsPerTier);
%     highRows = (2 * binsPerTier + 1):numSpeedBins;
%     
%     x_pos = 1:numPosBins;
%     
%     [lowLine,  lowSEM]  = getProfileStats(tuningSurface, lowRows);
%     [medLine,  medSEM]  = getProfileStats(tuningSurface, medRows);
%     [highLine, highSEM] = getProfileStats(tuningSurface, highRows);
%     
%     lowThreshold  = speedEdges(max(lowRows) + 1);
%     highThreshold = speedEdges(max(medRows) + 1);
%     
%     figHandle = figure('Name', sprintf('%s - Bouton %d', sessionName, targetROI), ...
%         'Position', [100 100 720 600], 'Color', 'w', 'Visible', 'off');
%     
%     leftMargin = 0.12;
%     plotWidth  = 0.62; 
%     
%     ax1 = subplot('Position', [leftMargin, 0.68, plotWidth, 0.24]); 
%     hold on;
%     
%     landmarkCentres = [40, 80, 120, 160];
%     for c = landmarkCentres
%         if c <= numPosBins
%             xline(c, ':', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2);
%         end
%     end
%     
%     colorLow  = [0, 0.6, 0.8];   
%     colorMed  = [0.2, 0.2, 0.2]; 
%     colorHigh = [0.8, 0, 0.6];   
%     
%     renderShadedError(ax1, x_pos, lowLine,  lowSEM,  colorLow);
%     renderShadedError(ax1, x_pos, medLine,  medSEM,  colorMed);
%     renderShadedError(ax1, x_pos, highLine, highSEM, colorHigh);
%     
%     legendTraces = []; legendLabels = {};
%     
%     if ~isempty(lowRows)
%         hLow = plot(ax1, x_pos, lowLine, 'Color', colorLow, 'LineWidth', 2);
%         legendTraces(end+1) = hLow; 
%         legendLabels{end+1} = sprintf('Low (<%.1f cm/s), n=%d rows', lowThreshold, length(lowRows));
%     end
%     if ~isempty(medRows)
%         hMed = plot(ax1, x_pos, medLine, 'Color', colorMed, 'LineWidth', 2);
%         legendTraces(end+1) = hMed; 
%         legendLabels{end+1} = sprintf('Med (%.1f-%.1f cm/s), n=%d rows', speedEdges(min(medRows)), highThreshold, length(medRows));
%     end
%     if ~isempty(highRows)
%         hHigh = plot(ax1, x_pos, highLine, 'Color', colorHigh, 'LineWidth', 2);
%         legendTraces(end+1) = hHigh; 
%         legendLabels{end+1} = sprintf('High (>%.1f cm/s), n=%d rows', speedEdges(min(highRows)), length(highRows));
%     end
%     
%     ylabel('\DeltaF/F (Neu)');
%     set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
%     title(sprintf('Bouton %d - Speed Position Activity (with Shaded SEM)\n %s', targetROI, response.stimName), 'FontWeight', 'bold');
%     
%     if ~isempty(legendTraces)
%         legend(legendTraces, legendLabels, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
%     end
%     if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax1,0); end
%     if exist('offsetAxes', 'file') == 2, offsetAxes(ax1); end
%     
%     ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.50]); 
%     
%     displayMatrix = padarray(tuningSurface, [1 1], 0, 'post');
%     xEdges = 0.5 : (numPosBins + 0.5);
%     
%     p = pcolor(xEdges, speedEdges, displayMatrix);
%     set(p, 'EdgeColor', 'none');
%     set(ax2, 'YDir', 'normal'); 
%     colormap(parula);
%     
%     activeData = tuningSurface(tuningSurface > 0);
%     maxVal = 1; if ~isempty(activeData), maxVal = prctile(activeData, 98.5); end
%     set(ax2, 'CLim', [0, maxVal]);
%     
%     c = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.025, 0.50]);
%     c.Label.String = '\DeltaF/F [NeuC]';
%     
%     hold on;
%     for c_vert = landmarkCentres
%         if c_vert <= numPosBins
%             xline(c_vert, '--w', 'Alpha', 0.3, 'LineWidth', 1.2);
%         end
%     end
%     
%     xticks([1 40 80 120 160 200]);
%     xticklabels({'1', '40', '80', '120', '160', '200'});
%     xlim([0.5, numPosBins + 0.5]);
%     
%     minY = speedEdges(1);
%     maxY = speedEdges(end);
%     ylim([minY, maxY]);
%     
%     set(ax2, 'TickDir', 'out');
%     yticks([1, 15, 30, 45]); 
%     yticklabels({'1', '15', '30', '45'});
%     
%     xlabel('Position (cm)'); ylabel('Running speed (cm/s)');
%     set(ax2, 'Box', 'off', 'FontSize', 11);
%     if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax2,0); end
%     
%     axStrat = axes('Position', [leftMargin + plotWidth + 0.11, 0.12, 0.02, 0.50]);
%     hold on;
%     
%     patch([0 1 1 0], [minY minY lowThreshold lowThreshold], [0.85, 0.95, 1], 'EdgeColor', 'none');
%     patch([0 1 1 0], [lowThreshold lowThreshold highThreshold highThreshold], [0.92, 0.92, 0.92], 'EdgeColor', 'none');
%     patch([0 1 1 0], [highThreshold highThreshold maxY maxY], [0.98, 0.85, 0.90], 'EdgeColor', 'none');
%     
%     set(axStrat, 'YLim', [minY, maxY], 'XLim', [0, 1]);
%     
%     yTickPlacements = [ ...
%         mean([minY, lowThreshold]), ...
%         mean([lowThreshold, highThreshold]), ...
%         mean([highThreshold, maxY]) ...
%     ];
%     
%     set(axStrat, 'YTick', yTickPlacements, ...
%                  'YTickLabel', {'Low', 'Med', 'High'}, 'YAxisLocation', 'right', ...
%                  'XTick', [], 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
%     ylabel(axStrat, 'Stratification Zones', 'FontSize', 10);
%     
%     linkaxes([ax1, ax2], 'x');
%     set(ax1, 'XTickLabel', []);
%     set(figHandle, 'Visible', 'on');
% end
% 
% function [mLine, semLine] = getProfileStats(surface, rowIndices)
%     if isempty(rowIndices)
%         mLine = nan(1, size(surface, 2));
%         semLine = nan(1, size(surface, 2));
%         return;
%     end
%     mLine = mean(surface(rowIndices, :), 1, 'omitnan');
%     nRows = length(rowIndices);
%     if nRows > 1
%         semLine = std(surface(rowIndices, :), 0, 1, 'omitnan') ./ sqrt(nRows);
%     else
%         semLine = zeros(1, size(surface, 2)); 
%     end
% end
% 
% function renderShadedError(ax, x, m, err, col)
%     if all(isnan(m)) || all(err == 0), return; end
%     x_patch = [x, fliplr(x)];
%     y_patch = [(m + err), fliplr(m - err)];
%     invalidData = isnan(y_patch);
%     x_patch(invalidData) = [];
%     y_patch(invalidData) = [];
%     patch(ax, x_patch, y_patch, col, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
% end

% function figHandle = plotSpeedPositionActivity_ForROI(response, targetROI, applySmoothing, smoothSigma)
%     % log scale 
%     if nargin < 4, applySmoothing = true; end
%     if nargin < 5, smoothSigma = [2.1, 2.5]; end
%     if nargin < 6, sessionName = 'session'; end
%     
%     tuningSurface = response.speedPositionActivity.matrix(:, :, targetROI);
%     speedCenters  = response.speedPositionActivity.speedBinCenters;
%     speedEdges    = response.speedPositionActivity.speedEdges;
%     [numSpeedBins, numPosBins] = size(tuningSurface);
%     
%     if applySmoothing
%         mask = ~isnan(tuningSurface);
%         dataZeroed = tuningSurface;
%         dataZeroed(~mask) = 0; 
%         blurredData = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
%         blurredMask = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
%         tuningSurface = blurredData ./ blurredMask;
%         tuningSurface(isnan(tuningSurface)) = 0;
%     else
%         tuningSurface(isnan(tuningSurface)) = 0;
%     end
%     
%     binsPerTier = 3; 
%     
%     lowRows  = 1:binsPerTier;
%     medRows  = (binsPerTier + 1):(2 * binsPerTier);
%     highRows = (2 * binsPerTier + 1):numSpeedBins;
%     
%     x_pos = 1:numPosBins;
%     
%     [lowLine,  lowSEM]  = getProfileStats(tuningSurface, lowRows);
%     [medLine,  medSEM]  = getProfileStats(tuningSurface, medRows);
%     [highLine, highSEM] = getProfileStats(tuningSurface, highRows);
%     
%     lowThreshold  = speedEdges(max(lowRows) + 1);
%     highThreshold = speedEdges(max(medRows) + 1);
%     
%     figHandle = figure('Name', sprintf('%s - Bouton %d', sessionName, targetROI), ...
%         'Position', [100 100 720 600], 'Color', 'w', 'Visible', 'off');
%     
%     leftMargin = 0.12;
%     plotWidth  = 0.62; 
%     
%     ax1 = subplot('Position', [leftMargin, 0.68, plotWidth, 0.24]); 
%     hold on;
%     
%     landmarkCentres = [40, 80, 120, 160];
%     for c = landmarkCentres
%         if c <= numPosBins
%             xline(c, ':', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2);
%         end
%     end
%     
%     colorLow  = [0, 0.6, 0.8];   
%     colorMed  = [0.2, 0.2, 0.2]; 
%     colorHigh = [0.8, 0, 0.6];   
%     
%     renderShadedError(ax1, x_pos, lowLine,  lowSEM,  colorLow);
%     renderShadedError(ax1, x_pos, medLine,  medSEM,  colorMed);
%     renderShadedError(ax1, x_pos, highLine, highSEM, colorHigh);
%     
%     legendTraces = []; legendLabels = {};
%     
%     if ~isempty(lowRows)
%         hLow = plot(ax1, x_pos, lowLine, 'Color', colorLow, 'LineWidth', 2);
%         legendTraces(end+1) = hLow; 
%         legendLabels{end+1} = sprintf('Low (<%.1f cm/s), n=%d rows', lowThreshold, length(lowRows));
%     end
%     if ~isempty(medRows)
%         hMed = plot(ax1, x_pos, medLine, 'Color', colorMed, 'LineWidth', 2);
%         legendTraces(end+1) = hMed; 
%         legendLabels{end+1} = sprintf('Med (%.1f-%.1f cm/s), n=%d rows', speedEdges(min(medRows)), highThreshold, length(medRows));
%     end
%     if ~isempty(highRows)
%         hHigh = plot(ax1, x_pos, highLine, 'Color', colorHigh, 'LineWidth', 2);
%         legendTraces(end+1) = hHigh; 
%         legendLabels{end+1} = sprintf('High (>%.1f cm/s), n=%d rows', speedEdges(min(highRows)), length(highRows));
%     end
%     
%     ylabel('\DeltaF/F (Neu)');
%     set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
%     title(sprintf('Bouton %d - Speed Position Activity (with Shaded SEM)\n %s', targetROI, response.stimName), 'FontWeight', 'bold');
%     
%     if ~isempty(legendTraces)
%         legend(legendTraces, legendLabels, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
%     end
%     if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax1,0); end
%     if exist('offsetAxes', 'file') == 2, offsetAxes(ax1); end
%     
%     ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.50]); 
%     
%     displayMatrix = padarray(tuningSurface, [1 1], 0, 'post');
%     xEdges = 0.5 : (numPosBins + 0.5);
%     
%     p = pcolor(xEdges, speedEdges, displayMatrix);
%     set(p, 'EdgeColor', 'none');
%     set(ax2, 'YScale', 'log', 'YDir', 'normal');
%     colormap(parula);
%     
%     activeData = tuningSurface(tuningSurface > 0);
%     maxVal = 1; if ~isempty(activeData), maxVal = prctile(activeData, 98.5); end
%     set(ax2, 'CLim', [0, maxVal]);
%     
%     c = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.025, 0.50]);
%     c.Label.String = '\DeltaF/F [NeuC]';
%     
%     hold on;
%     for c_vert = landmarkCentres
%         if c_vert <= numPosBins
%             xline(c_vert, '--w', 'Alpha', 0.3, 'LineWidth', 1.2);
%         end
%     end
%     
%     xticks([1 40 80 120 160 200]);
%     xticklabels({'1', '40', '80', '120', '160', '200'});
%     xlim([0.5, numPosBins + 0.5]);
%     
%     minY = speedEdges(1);
%     maxY = speedEdges(end);
%     ylim([minY, maxY]);
%     
%     set(ax2, 'YMinorTick', 'on', 'TickDir', 'out');
%     yticks([2, 5, 10, 20, 30, 45]); 
%     yticklabels({'2', '5', '10', '20', '30', '45'});
%     
%     xlabel('Position (cm)'); ylabel('Running speed (cm/s)');
%     set(ax2, 'Box', 'off', 'FontSize', 11);
%     if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax2,0); end
%     
%     axStrat = axes('Position', [leftMargin + plotWidth + 0.11, 0.12, 0.02, 0.50]);
%     hold on;
%     
%     patch([0 1 1 0], [minY minY lowThreshold lowThreshold], [0.85, 0.95, 1], 'EdgeColor', 'none');
%     patch([0 1 1 0], [lowThreshold lowThreshold highThreshold highThreshold], [0.92, 0.92, 0.92], 'EdgeColor', 'none');
%     patch([0 1 1 0], [highThreshold highThreshold maxY maxY], [0.98, 0.85, 0.90], 'EdgeColor', 'none');
%     
%     set(axStrat, 'YScale', 'log', 'YLim', [minY, maxY], 'XLim', [0, 1]);
%     
%     yTickPlacements = [ ...
%         exp(mean([log(minY), log(lowThreshold)])), ...
%         exp(mean([log(lowThreshold), log(highThreshold)])), ...
%         exp(mean([log(highThreshold), log(maxY)])) ...
%     ];
%     
%     set(axStrat, 'YTick', yTickPlacements, ...
%                  'YTickLabel', {'Low', 'Med', 'High'}, 'YAxisLocation', 'right', ...
%                  'XTick', [], 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
%     ylabel(axStrat, 'Stratification Zones', 'FontSize', 10);
%     
%     linkaxes([ax1, ax2], 'x');
%     set(ax1, 'XTickLabel', []);
%     set(figHandle, 'Visible', 'on');
% end
% 
% function [mLine, semLine] = getProfileStats(surface, rowIndices)
%     if isempty(rowIndices)
%         mLine = nan(1, size(surface, 2));
%         semLine = nan(1, size(surface, 2));
%         return;
%     end
%     mLine = mean(surface(rowIndices, :), 1, 'omitnan');
%     nRows = length(rowIndices);
%     if nRows > 1
%         semLine = std(surface(rowIndices, :), 0, 1, 'omitnan') ./ sqrt(nRows);
%     else
%         semLine = zeros(1, size(surface, 2)); 
%     end
% end
% 
% function renderShadedError(ax, x, m, err, col)
%     if all(isnan(m)) || all(err == 0), return; end
%     x_patch = [x, fliplr(x)];
%     y_patch = [(m + err), fliplr(m - err)];
%     invalidData = isnan(y_patch);
%     x_patch(invalidData) = [];
%     y_patch(invalidData) = [];
%     patch(ax, x_patch, y_patch, col, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
% end
% 






% function plotSpeedPositionActivity_ForROI(sessionFileInfo, response, targetROI, applySmoothing, smoothSigma)
% % plotspeedpositionactivity_forroi - plots a 2d heatmap with a shared-axis 1d profile.
% % ensures perfect vertical alignment by accounting for the colorbar width.
% 
%     %% set default options
%     if nargin < 4, applySmoothing = true; end
%     if nargin < 5, smoothSigma = [1.1, 1.5]; end
%     
%     %% 2. extract data
%     tuningSurface = response.speedPositionActivity.matrix(:, :, targetROI);
%     speedCenters  = response.speedPositionActivity.speedBinCenters;
%     numPosBins    = size(tuningSurface, 2);
% 
%     %% smoothing logic 
%     if applySmoothing
%         mask = ~isnan(tuningSurface);
%         dataZeroed = tuningSurface;
%         dataZeroed(~mask) = 0; 
%         blurredData = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
%         blurredMask = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
%         tuningSurface = blurredData ./ blurredMask;
%         tuningSurface(isnan(tuningSurface)) = 0;
%     else
%         tuningSurface(isnan(tuningSurface)) = 0;
%     end
% 
%     %% create figure
%     figHandle = figure('Name', sprintf('Soma %d', targetROI), ...
%         'Position', [100 100 650 550], 'Color', 'w');
%     
%     % define a shared width and left margin to ensure alignment
%     leftMargin = 0.13;
%     plotWidth  = 0.70; 
% 
%     % collapsed spatial profile ---
%     ax1 = subplot('Position', [leftMargin, 0.72, plotWidth, 0.20]); 
%     collapsedSpatial = mean(tuningSurface, 1, 'omitnan');
%     plot(1:numPosBins, collapsedSpatial, 'k', 'LineWidth', 1.5);
%     hold on;
%     xline([40, 80, 120, 160, 200], '--', 'Color', [0.7 0.7 0.7], 'Alpha', 0.5);
%     
%     ylabel('\DeltaF/F');
%     set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
%     title(sprintf('Bouton %d: speed-position activity', targetROI));
% 
%     %speed-position heatmap
%     ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.55]); 
%     imagesc(1:numPosBins, speedCenters, tuningSurface);
%     set(ax2, 'YDir', 'normal'); 
%     
%     colormap(parula);
%     
%     % contrast (98.5th percentile)
%     activeData = tuningSurface(tuningSurface > 0);
%     if ~isempty(activeData)
%         maxVal = prctile(activeData, 98.5);
%     else
%         maxVal = 1;
%     end
%     set(ax2, 'CLim', [0, maxVal]);
% 
%     % add colorbar WITHOUT resizing the axis
%     c = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.03, 0.55]);
%     c.Label.String = '\DeltaF/F [NeuC]';
% 
%     %% reference lines and axes
%     hold on;
%     xline([40, 80, 120, 160, 200], '--w', 'Alpha', 0.4, 'LineWidth', 1.5);
%     
%     % x-axis formatting
%     xticks([1 40 80 120 160 200]);
%     xticklabels({'1', '40', '80', '120', '160', '200'});
%     xlim([1, numPosBins]);
% 
%     % y-axis formatting (log scale)
%     set(ax2, 'YScale', 'log'); 
%     set(ax2, 'YMinorTick', 'on', 'TickDir', 'out');
%     yticks([2, 5, 10, 20, 30]); 
%     yticklabels({'2', '5', '10', '20', '30'});
%     ylim([min(speedCenters), max(speedCenters)]);
% 
%     % labels and styling
%     xlabel('Position (cm)');
%     ylabel('Running speed (cm/s)');
%     set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
% 
%     %% link the x-axes for interactive use
%     linkaxes([ax1, ax2], 'x');
% 
%     %% 7. directory setup and saving
%     figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'SpeedPositionMaps');
%     if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
%     
%     smoothLabel = 'raw'; if applySmoothing, smoothLabel = 'smoothed'; end
%     saveName = sprintf('%s_%s_bouton%03d_%s_speedpos.png', ...
%         sessionFileInfo.animal_name, sessionFileInfo.session_name, targetROI, smoothLabel);
%     
%     exportgraphics(figHandle, fullfile(figSaveDir, saveName), 'Resolution', 300);
%     fprintf('figure saved for roi %d: %s\n', targetROI, saveName);
% end


% function figHandle = plotSpeedPositionActivity_ForROI(sessionFileInfo, response, targetROI, applySmoothing, smoothSigma)
%     if nargin < 4, applySmoothing = true; end
%     if nargin < 5, smoothSigma = [1.1, 1.5]; end
% 
%     tuningSurface = response.speedPositionActivity.matrix(:, :, targetROI);
%     speedCenters  = response.speedPositionActivity.speedBinCenters;
%     speedEdges    = response.speedPositionActivity.speedEdges;
%     [numSpeedBins, numPosBins] = size(tuningSurface);
% 
%     if numSpeedBins ~= 10
%         error('This script is optimized specifically to parse a 10-bin speed matrix.');
%     end
% 
%     if applySmoothing
%         mask = ~isnan(tuningSurface);
%         dataZeroed = tuningSurface;
%         dataZeroed(~mask) = 0; 
%         blurredData = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
%         blurredMask = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
%         tuningSurface = blurredData ./ blurredMask;
%         tuningSurface(isnan(tuningSurface)) = 0;
%     else
%         tuningSurface(isnan(tuningSurface)) = 0;
%     end
% 
%     figHandle = figure('Name', sprintf('Bouton %d', targetROI), ...
%         'Position', [100 100 720 600], 'Color', 'w');
% 
%     leftMargin = 0.12;
%     plotWidth  = 0.62; 
% 
%     ax1 = subplot('Position', [leftMargin, 0.68, plotWidth, 0.24]); 
%     hold on;
% 
%     landmarkCentres = [40, 80, 120, 160];
%     for c = landmarkCentres
%         if c <= numPosBins
%             xline(c, ':', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2);
%         end
%     end
% 
%     x_pos = 1:numPosBins;
% 
%     lowLine  = mean(tuningSurface(1:3, :), 1, 'omitnan');
%     medLine  = mean(tuningSurface(4:7, :), 1, 'omitnan');
%     highLine = mean(tuningSurface(8:10, :), 1, 'omitnan');
% 
%     lowThreshold  = speedEdges(4);
%     highThreshold = speedEdges(8);
% 
%     plot(x_pos, lowLine, 'Color', [0, 0.9, 1], 'LineWidth', 2);
%     plot(x_pos, medLine, 'Color', 'k', 'LineWidth', 2);
%     plot(x_pos, highLine, 'Color', [1, 0, 0.9], 'LineWidth', 2);
% 
%     ylabel('\DeltaF/F (Neu)');
%     set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
%     title(sprintf('Bouton %d - Stratified Speed Profile', targetROI), 'FontWeight', 'normal');
% 
%     labels = { ...
%         sprintf('Low (<%.1f cm/s)', lowThreshold), ...
%         sprintf('Med (%.1f-%.1f cm/s)', lowThreshold, highThreshold), ...
%         sprintf('High (>%.1f cm/s)', highThreshold) ...
%     };
%     legend(ax1, labels, 'Location', 'northeast', 'Box', 'off', 'FontSize', 9);
% 
%     ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.50]); 
%     imagesc(1:numPosBins, speedCenters, tuningSurface);
%     set(ax2, 'YDir', 'normal'); 
%     colormap(parula);
% 
%     activeData = tuningSurface(tuningSurface > 0);
%     if ~isempty(activeData)
%         maxVal = prctile(activeData, 98.5);
%     else
%         maxVal = 1;
%     end
%     set(ax2, 'CLim', [0, maxVal]);
% 
%     c = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.025, 0.50]);
%     c.Label.String = '\DeltaF/F [NeuC]';
% 
%     hold on;
%     for c_vert = landmarkCentres
%         if c_vert <= numPosBins
%             xline(c_vert, '--w', 'Alpha', 0.3, 'LineWidth', 1.2);
%         end
%     end
% 
%     xticks([1 40 80 120 160 200]);
%     xticklabels({'1', '40', '80', '120', '160', '200'});
%     xlim([1, numPosBins]);
% 
%     set(ax2, 'YScale', 'log'); 
%     set(ax2, 'YMinorTick', 'on', 'TickDir', 'out');
%     yticks([2, 5, 10, 20, 30]); 
%     yticklabels({'2', '5', '10', '20', '30'});
%     minY = min(speedCenters);
%     maxY = max(speedCenters);
%     ylim([minY, maxY]);
% 
%     xlabel('Position (cm)');
%     ylabel('Running speed (cm/s)');
%     set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
% 
%     axStrat = axes('Position', [leftMargin + plotWidth + 0.11, 0.12, 0.02, 0.50]);
%     hold on;
% 
%     patch([0 1 1 0], [minY minY lowThreshold lowThreshold], [0.68, 0.92, 0.98], 'EdgeColor', 'none');
%     patch([0 1 1 0], [lowThreshold lowThreshold highThreshold highThreshold], [0.88, 0.88, 0.88], 'EdgeColor', 'none');
%     patch([0 1 1 0], [highThreshold highThreshold maxY maxY], [0.96, 0.64, 0.76], 'EdgeColor', 'none');
% 
%     set(axStrat, 'YScale', 'log', 'YLim', [minY, maxY], 'XLim', [0, 1]);
%     set(axStrat, 'YTick', [mean([minY, lowThreshold]), mean([lowThreshold, highThreshold]), mean([highThreshold, maxY])], ...
%                  'YTickLabel', {'Low', 'Med', 'High'}, 'YAxisLocation', 'right', ...
%                  'XTick', [], 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
%     ylabel(axStrat, 'Stratification Zones', 'FontSize', 10);
% 
%     linkaxes([ax1, ax2], 'x');
% 
%     figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'SpeedPositionMaps');
%     if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
% 
%     smoothLabel = 'raw'; if applySmoothing, smoothLabel = 'smoothed'; end
%     saveName = sprintf('%s_%s_bouton%03d_%s_zone_labeled.png', ...
%         sessionFileInfo.animal_name, sessionFileInfo.session_name, targetROI, smoothLabel);
% 
%     exportgraphics(figHandle, fullfile(figSaveDir, saveName), 'Resolution', 300);
%     fprintf('Figure saved for ROI %d: %s\n', targetROI, saveName);
% end
% 
% function figHandle = plotSpeedPositionActivity_ForROI(response, targetROI, applySmoothing, smoothSigma)
%     % Plots a 2D Position-Speed heatmap alongside 1D stratified profiles 
%     % with shaded SEM bands based on fixed physical speed thresholds.
%       % this can be used with the log physical measures 
%     
%     if nargin < 4, applySmoothing = true; end
% 
%     if nargin < 5, smoothSigma = [1.1, 1.5]; end
%     if nargin < 6, sessionName = 'session'; end
%     
%     %  copy matrix variables from response struct
%     tuningSurface = response.speedPositionActivity.matrix(:, :, targetROI);
%     speedCenters  = response.speedPositionActivity.speedBinCenters;
%     speedEdges    = response.speedPositionActivity.speedEdges;
%     [~, numPosBins] = size(tuningSurface);
%     
%     %  Spatial Smoothing (dealing with NaNs natively)
%     if applySmoothing
%         mask = ~isnan(tuningSurface);
%         dataZeroed = tuningSurface;
%         dataZeroed(~mask) = 0; 
%         blurredData = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
%         blurredMask = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
%         tuningSurface = blurredData ./ blurredMask;
%         tuningSurface(isnan(tuningSurface)) = 0;
%     else
%         tuningSurface(isnan(tuningSurface)) = 0;
%     end
%     
%     FIXED_LOW_CUTOFF  = 15; % Upper boundary of low speed 
%     FIXED_HIGH_CUTOFF = 27; % Upper boundary of medium speed 
%     
%     % Locate which matrix rows naturally belong to each physical velocity tier
%     lowRows  = find(speedCenters < FIXED_LOW_CUTOFF);
%     medRows  = find(speedCenters >= FIXED_LOW_CUTOFF & speedCenters <= FIXED_HIGH_CUTOFF);
%     highRows = find(speedCenters > FIXED_HIGH_CUTOFF);
%     
%     x_pos = 1:numPosBins;
%     
%     % Compute Mean and SEM across pooled matrix rows for each spatial bin
%     [lowLine,  lowSEM]  = getProfileStats(tuningSurface, lowRows);
%     [medLine,  medSEM]  = getProfileStats(tuningSurface, medRows);
%     [highLine, highSEM] = getProfileStats(tuningSurface, highRows);
%     
%     % Track precise physical transitions based on the active row edges
%     if ~isempty(lowRows)
%         lowThreshold = speedEdges(max(lowRows) + 1);
%     else
%         lowThreshold = min(speedEdges);
%     end
%     
%     if ~isempty(medRows)
%         highThreshold = speedEdges(max(medRows) + 1);
%     elseif ~isempty(lowRows)
%         highThreshold = speedEdges(max(lowRows) + 1);
%     else
%         highThreshold = min(speedEdges);
%     end
%  
%     
%     % 
%     figHandle = figure('Name', sprintf('%s - Bouton %d', sessionName, targetROI), ...
%         'Position', [100 100 720 600], 'Color', 'w', 'Visible', 'off');
%     
%     leftMargin = 0.12;
%     plotWidth  = 0.62; 
%     
%     %% 1D line plots
%     ax1 = subplot('Position', [leftMargin, 0.68, plotWidth, 0.24]); 
%     hold on;
%     
%     % Spatial position landmark indicator lines
%     landmarkCentres = [40, 80, 120, 160];
%     for c = landmarkCentres
%         if c <= numPosBins
%             xline(c, ':', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2);
%         end
%     end
%     
%     colorLow  = [0, 0.9, 1];   % Cyan
%     colorMed  = [0, 0, 0];     % Black
%     colorHigh = [1, 0, 0.9];   % Magenta
%     
%     % Render Shaded SEM Clouds FIRST (so solid lines sit cleanly on top)
%     renderShadedError(ax1, x_pos, lowLine,  lowSEM,  colorLow);
%     renderShadedError(ax1, x_pos, medLine,  medSEM,  colorMed);
%     renderShadedError(ax1, x_pos, highLine, highSEM, colorHigh);
%     
%     legendTraces = []; legendLabels = {};
%     
%     % Plot solid traces ONLY if rows exist in that velocity segment
%     if ~isempty(lowRows)
%         hLow = plot(ax1, x_pos, lowLine, 'Color', colorLow, 'LineWidth', 2);
%         legendTraces(end+1) = hLow; 
%         legendLabels{end+1} = sprintf('Low (<%.1f cm/s), n=%d rows', lowThreshold, length(lowRows));
%     end
%     if ~isempty(medRows)
%         hMed = plot(ax1, x_pos, medLine, 'Color', colorMed, 'LineWidth', 2);
%         legendTraces(end+1) = hMed; 
%         legendLabels{end+1} = sprintf('Med (%.1f-%.1f cm/s), n=%d rows', speedEdges(min(medRows)), highThreshold, length(medRows));
%     end
%     if ~isempty(highRows)
%         hHigh = plot(ax1, x_pos, highLine, 'Color', colorHigh, 'LineWidth', 2);
%         legendTraces(end+1) = hHigh; 
%         legendLabels{end+1} = sprintf('High (>%.1f cm/s), n=%d rows', speedEdges(min(highRows)), length(highRows));
%     end
%     
%     ylabel('\DeltaF/F (Neu)');
%     set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
%     title(sprintf('Bouton %d - Speed Position Activity (with Shaded SEM)\n %s', targetROI, response.stimName), 'FontWeight', 'bold');
%     
%     if ~isempty(legendTraces)
%         legend(legendTraces, legendLabels, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
%     end
%     defaultAxesProperties(ax1,0)
%     offsetAxes(ax1)
%     
%     %% 2D Position-Speed Heatmap (bottom left) ---
%     ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.50]); 
%     imagesc(1:numPosBins, speedCenters, tuningSurface);
%     set(ax2, 'YDir', 'normal'); 
%     colormap(parula);
%     
%     activeData = tuningSurface(tuningSurface > 0);
%     maxVal = 1; if ~isempty(activeData), maxVal = prctile(activeData, 98.5); end
%     set(ax2, 'CLim', [0, maxVal]);
%     
%     c = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.025, 0.50]);
%     c.Label.String = '\DeltaF/F [NeuC]';
%     
%     hold on;
%     for c_vert = landmarkCentres
%         if c_vert <= numPosBins
%             xline(c_vert, '--w', 'Alpha', 0.3, 'LineWidth', 1.2);
%         end
%     end
%     
%     xticks([1 40 80 120 160 200]);
%     xticklabels({'1', '40', '80', '120', '160', '200'});
%     xlim([1, numPosBins]);
%     
%     set(ax2, 'YScale', 'log', 'YMinorTick', 'on', 'TickDir', 'out');
%     yticks([2, 5, 10, 20, 30]); yticklabels({'2', '5', '10', '20', '30'});
%     
%     % FIXED: Set limits using the exact outer speed edges rather than bin centers
%     minY = speedEdges(1); 
%     maxY = speedEdges(end);
%     ylim([minY, maxY]);
%     
%     xlabel('Position (cm)'); ylabel('Running speed (cm/s)');
%     set(ax2, 'Box', 'off', 'FontSize', 11);
%     defaultAxesProperties(ax2,0)
% %     offsetAxes(ax2)
%     
%     %% -Subplot 3: vertical color stratification strip (Far-Right Edge) ---
%     axStrat = axes('Position', [leftMargin + plotWidth + 0.11, 0.12, 0.02, 0.50]);
%     hold on;
%     
%     patch([0 1 1 0], [minY minY lowThreshold lowThreshold], [0.68, 0.92, 0.98], 'EdgeColor', 'none');
%     patch([0 1 1 0], [lowThreshold lowThreshold highThreshold highThreshold], [0.88, 0.88, 0.88], 'EdgeColor', 'none');
%     patch([0 1 1 0], [highThreshold highThreshold maxY maxY], [0.96, 0.64, 0.76], 'EdgeColor', 'none');
%     
%     set(axStrat, 'YScale', 'log', 'YLim', [minY, maxY], 'XLim', [0, 1]);
%     
%     yTickPlacements = [mean([minY, lowThreshold]), mean([lowThreshold, highThreshold]), mean([highThreshold, maxY])];
%     yTickPlacements(isnan(yTickPlacements)) = minY; 
%     
%     set(axStrat, 'YTick', yTickPlacements, ...
%                  'YTickLabel', {'Low', 'Med', 'High'}, 'YAxisLocation', 'right', ...
%                  'XTick', [], 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
%     ylabel(axStrat, 'Stratification Zones', 'FontSize', 10);
%     
%     linkaxes([ax1, ax2], 'x');
%     set(ax1, 'XTickLabel', []);
%     set(figHandle, 'Visible', 'on');
% end

% function [mLine, semLine] = getProfileStats(surface, rowIndices)
%     if isempty(rowIndices)
%         mLine = nan(1, size(surface, 2));
%         semLine = nan(1, size(surface, 2));
%         return;
%     end
%     mLine = mean(surface(rowIndices, :), 1, 'omitnan');
%     nRows = length(rowIndices);
%     if nRows > 1
%         semLine = std(surface(rowIndices, :), 0, 1, 'omitnan') ./ sqrt(nRows);
%     else
%         semLine = zeros(1, size(surface, 2)); 
%     end
% end
% 
% % i think this is an isse becaue it caps it at the bin cenre and does it mean i cant see the last bin? 
% function renderShadedError(ax, x, m, err, col)
%     if all(isnan(m)) || all(err == 0), return; end
%     x_patch = [x, fliplr(x)];
%     y_patch = [(m + err), fliplr(m - err)];
%     invalidData = isnan(y_patch);
%     x_patch(invalidData) = [];
%     y_patch(invalidData) = [];
%     patch(ax, x_patch, y_patch, col, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
% end


% function figHandle = plotSpeedPositionActivity_ForROI(response, targetROI, applySmoothing, smoothSigma)
%     % Plots a 2D Position-Speed heatmap alongside 1D stratified profiles 
%     % with shaded SEM bands using a true logarithmic speed axis.
%     
%     if nargin < 4, applySmoothing = true; end
%     if nargin < 5, smoothSigma = [2.1, 2.5]; end
%     if nargin < 6, sessionName = 'session'; end
%     
%     %  copy matrix variables from response struct
%     tuningSurface = response.speedPositionActivity.matrix(:, :, targetROI);
%     speedCenters  = response.speedPositionActivity.speedBinCenters;
%     speedEdges    = response.speedPositionActivity.speedEdges;
%     [numSpeedBins, numPosBins] = size(tuningSurface);
%     
%     %  Spatial Smoothing (dealing with NaNs natively)
%     if applySmoothing
%         mask = ~isnan(tuningSurface);
%         dataZeroed = tuningSurface;
%         dataZeroed(~mask) = 0; 
%         blurredData = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
%         blurredMask = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
%         tuningSurface = blurredData ./ blurredMask;
%         tuningSurface(isnan(tuningSurface)) = 0;
%     else
%         tuningSurface(isnan(tuningSurface)) = 0;
%     end
%     
%     % -------------------------------------------------------------------------
%     % FIXED TIER ROW SPLITTING (3 ROWS PER VELOCITY TIER)
%     % -------------------------------------------------------------------------
%     binsPerTier = 3; 
%     
%     lowRows  = 1:binsPerTier;
%     medRows  = (binsPerTier + 1):(2 * binsPerTier);
%     highRows = (2 * binsPerTier + 1):numSpeedBins;
%     
%     x_pos = 1:numPosBins;
%     
%     % Compute Mean and SEM across pooled matrix rows for each spatial bin
%     [lowLine,  lowSEM]  = getProfileStats(tuningSurface, lowRows);
%     [medLine,  medSEM]  = getProfileStats(tuningSurface, medRows);
%     [highLine, highSEM] = getProfileStats(tuningSurface, highRows);
%     
%     % Track precise physical transitions based on the active row edges
%     lowThreshold  = speedEdges(max(lowRows) + 1);
%     highThreshold = speedEdges(max(medRows) + 1);
%     % -------------------------------------------------------------------------
%     
%     figHandle = figure('Name', sprintf('%s - Bouton %d', sessionName, targetROI), ...
%         'Position', [100 100 720 600], 'Color', 'w', 'Visible', 'off');
%     
%     leftMargin = 0.12;
%     plotWidth  = 0.62; 
%     
%     %% 1D line plots
%     ax1 = subplot('Position', [leftMargin, 0.68, plotWidth, 0.24]); 
%     hold on;
%     
%     % Spatial position landmark indicator lines
%     landmarkCentres = [40, 80, 120, 160];
%     for c = landmarkCentres
%         if c <= numPosBins
%             xline(c, ':', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2);
%         end
%     end
%     
%     colorLow  = [0, 0.6, 0.8];   
%     colorMed  = [0.2, 0.2, 0.2]; 
%     colorHigh = [0.8, 0, 0.6];   
%     
%     renderShadedError(ax1, x_pos, lowLine,  lowSEM,  colorLow);
%     renderShadedError(ax1, x_pos, medLine,  medSEM,  colorMed);
%     renderShadedError(ax1, x_pos, highLine, highSEM, colorHigh);
%     
%     legendTraces = []; legendLabels = {};
%     
%     if ~isempty(lowRows)
%         hLow = plot(ax1, x_pos, lowLine, 'Color', colorLow, 'LineWidth', 2);
%         legendTraces(end+1) = hLow; 
%         legendLabels{end+1} = sprintf('Low (<%.1f cm/s), n=%d rows', lowThreshold, length(lowRows));
%     end
%     if ~isempty(medRows)
%         hMed = plot(ax1, x_pos, medLine, 'Color', colorMed, 'LineWidth', 2);
%         legendTraces(end+1) = hMed; 
%         legendLabels{end+1} = sprintf('Med (%.1f-%.1f cm/s), n=%d rows', speedEdges(min(medRows)), highThreshold, length(medRows));
%     end
%     if ~isempty(highRows)
%         hHigh = plot(ax1, x_pos, highLine, 'Color', colorHigh, 'LineWidth', 2);
%         legendTraces(end+1) = hHigh; 
%         legendLabels{end+1} = sprintf('High (>%.1f cm/s), n=%d rows', speedEdges(min(highRows)), length(highRows));
%     end
%     
%     ylabel('\DeltaF/F (Neu)');
%     set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
%     title(sprintf('Bouton %d - Speed Position Activity (with Shaded SEM)\n %s', targetROI, response.stimName), 'FontWeight', 'bold');
%     
%     if ~isempty(legendTraces)
%         legend(legendTraces, legendLabels, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
%     end
%     if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax1,0); end
%     if exist('offsetAxes', 'file') == 2, offsetAxes(ax1); end
%     
%     %% 2D Position-Speed Heatmap (True Log Scale) ---
%     ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.50]); 
%     
%     % Using pcolor with explicit edges guarantees perfect edge-to-edge flushness on log axes
%     % We append a dummy row/col because pcolor drops the last row/col of the matrix data
%     displayMatrix = padarray(tuningSurface, [1 1], 0, 'post');
%     xEdges = 0.5 : (numPosBins + 0.5);
%     
%     p = pcolor(xEdges, speedEdges, displayMatrix);
%     set(p, 'EdgeColor', 'none');
%     set(ax2, 'YScale', 'log', 'YDir', 'normal');
%     colormap(parula);
%     
%     activeData = tuningSurface(tuningSurface > 0);
%     maxVal = 1; if ~isempty(activeData), maxVal = prctile(activeData, 98.5); end
%     set(ax2, 'CLim', [0, maxVal]);
%     
%     c = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.025, 0.50]);
%     c.Label.String = '\DeltaF/F [NeuC]';
%     
%     hold on;
%     for c_vert = landmarkCentres
%         if c_vert <= numPosBins
%             xline(c_vert, '--w', 'Alpha', 0.3, 'LineWidth', 1.2);
%         end
%     end
%     
%     xticks([1 40 80 120 160 200]);
%     xticklabels({'1', '40', '80', '120', '160', '200'});
%     xlim([0.5, numPosBins + 0.5]);
%     
%     minY = speedEdges(1);
%     maxY = speedEdges(end);
%     ylim([minY, maxY]);
%     
%     % True log ticks clean layout
%     set(ax2, 'YMinorTick', 'on', 'TickDir', 'out');
%     yticks([2, 5, 10, 20, 30, 45]); 
%     yticklabels({'2', '5', '10', '20', '30', '45'});
%     
%     xlabel('Position (cm)'); ylabel('Running speed (cm/s)');
%     set(ax2, 'Box', 'off', 'FontSize', 11);
%     if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax2,0); end
%     
%     %% -Subplot 3: vertical color stratification strip (Far-Right Edge) ---
%     axStrat = axes('Position', [leftMargin + plotWidth + 0.11, 0.12, 0.02, 0.50]);
%     hold on;
%     
%     patch([0 1 1 0], [minY minY lowThreshold lowThreshold], [0.85, 0.95, 1], 'EdgeColor', 'none');
%     patch([0 1 1 0], [lowThreshold lowThreshold highThreshold highThreshold], [0.92, 0.92, 0.92], 'EdgeColor', 'none');
%     patch([0 1 1 0], [highThreshold highThreshold maxY maxY], [0.98, 0.85, 0.90], 'EdgeColor', 'none');
%     
%     set(axStrat, 'YScale', 'log', 'YLim', [minY, maxY], 'XLim', [0, 1]);
%     
%     % Ticks placed log-arithmetically centered within each patch zone
%     yTickPlacements = [
%         exp(mean([log(minY), log(lowThreshold)])), ...
%         exp(mean([log(lowThreshold), log(highThreshold)])), ...
%         exp(mean([log(highThreshold), log(maxY)]))
%     ];
%     
%     set(axStrat, 'YTick', yTickPlacements, ...
%                  'YTickLabel', {'Low', 'Med', 'High'}, 'YAxisLocation', 'right', ...
%                  'XTick', [], 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
%     ylabel(axStrat, 'Stratification Zones', 'FontSize', 10);
%     
%     linkaxes([ax1, ax2], 'x');
%     set(ax1, 'XTickLabel', []);
%     set(figHandle, 'Visible', 'on');
% end
% 
% function [mLine, semLine] = getProfileStats(surface, rowIndices)
%     if isempty(rowIndices)
%         mLine = nan(1, size(surface, 2));
%         semLine = nan(1, size(surface, 2));
%         return;
%     end
%     mLine = mean(surface(rowIndices, :), 1, 'omitnan');
%     nRows = length(rowIndices);
%     if nRows > 1
%         semLine = std(surface(rowIndices, :), 0, 1, 'omitnan') ./ sqrt(nRows);
%     else
%         semLine = zeros(1, size(surface, 2)); 
%     end
% end
% 
% function renderShadedError(ax, x, m, err, col)
%     if all(isnan(m)) || all(err == 0), return; end
%     x_patch = [x, fliplr(x)];
%     y_patch = [(m + err), fliplr(m - err)];
%     invalidData = isnan(y_patch);
%     x_patch(invalidData) = [];
%     y_patch(invalidData) = [];
%     patch(ax, x_patch, y_patch, col, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
%end