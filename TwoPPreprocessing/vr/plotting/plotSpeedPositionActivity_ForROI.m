function figHandle = plotSpeedPositionActivity_ForROI(response, targetROI, applySmoothing, smoothSigma)
% Plots panel b (speed-stratified position tuning curves, low vs high only)
% and panel c (continuous 2D speed x position heatmap) for a single ROI.
%
% Panel b's y-axis is pinned to the same range used for panel c's color
% scale (CLim), so the two panels are directly visually comparable — a
% line height in panel b sits at the same ΔF/F value as the equivalent
% color intensity in panel c. The shared range uses a 1.5th-98.5th
% percentile clip on the heatmap data (not the true min/max — a few
% noisy, likely low-lap-count outlier cells otherwise dominate the scale
% and wash out real structure), extended if needed to cover panel b's own
% line peaks/troughs. This does NOT make panel b's y-axis a substitute
% for the colorbar: the colorbar decodes COLOR into a value; panel b's
% axis only encodes LINE HEIGHT. Matching the ranges makes
% cross-referencing easier, it doesn't let you read a color off an axis
% with no color information on it. The colorbar is kept.
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

%% Shared colour scale — computed once, used for both panel b's y-limits
% and panel c's CLim, so the two panels are on the same ΔF/F range.
% Uses a percentile clip (1.5th-98.5th) rather than the true min/max —
% the true extremes are dominated by one or a few noisy, likely
% low-lap-count outlier cells (check occupancy/lap counts for the
% continuous matrix if you want to confirm and mask those directly rather
% than clip around them). Clipping keeps real landmark structure visible
% instead of being washed out by stretching the scale to fit a handful of
% unreliable bins. This still keeps the ORIGINAL goal — panel b and panel
% c share one consistent, comparable scale — just not the unclipped
% extremes, which weren't actually representative of the data.
% Computed BEFORE the NaN->0 replacement below, using omitnan directly,
% so genuine zero-valued bins aren't confused with the former-NaN bins.
validSurfaceVals = tuningSurface(~isnan(tuningSurface));
minVal = 0;
maxVal = 1;
if ~isempty(validSurfaceVals)
    minVal = prctile(validSurfaceVals, 1.5);
    maxVal = prctile(validSurfaceVals, 98.5);
end

% Ensure panel b's own line peaks/troughs are covered too — lowLine/highLine
% come from a different computation than tuningSurface, so their range
% isn't guaranteed to fall inside tuningSurface's clipped range. Using the
% true min/max of the lines themselves (not clipped) since they're only
% two smoothed curves, not hundreds of individual cells — far less prone
% to single-outlier domination than the heatmap.
minVal = min(minVal, min([lowLine(:); highLine(:)], [], 'omitnan'));
maxVal = max(maxVal, max([lowLine(:); highLine(:)], [], 'omitnan'));

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
set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins], 'YLim', [minVal, maxVal]);
title(sprintf('ROI %d — %s', targetROI, response.stimName), 'FontWeight', 'bold');

if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax1, 0); end
if exist('offsetAxes',            'file') == 2, offsetAxes(ax1); end

% Re-assert YLim — defaultAxesProperties/offsetAxes may auto-scale or
% round the axis limits internally, which can silently override the
% [minVal, maxVal] set above. Setting it again here guarantees it sticks.
ylim(ax1, [minVal, maxVal]);

% Explicit shared ticks — computed once from [minVal, maxVal] and applied
% identically to panel b's y-axis AND the colorbar below, so the two
% don't drift apart via two independent auto-tick generators. Rounded to
% multiples of tickStep (0.2) rather than raw linspace values, so ticks
% read as clean numbers (0, 0.2, 0.4, ...) instead of arbitrary decimals.
tickStep = 0.2; % change this to whatever increment you want
firstTick = ceil(minVal / tickStep) * tickStep;
lastTick  = floor(maxVal / tickStep) * tickStep;
sharedTicks = firstTick:tickStep:lastTick;
if isempty(sharedTicks) % range narrower than one tickStep — fall back to the two endpoints
    sharedTicks = [minVal, maxVal];
end
set(ax1, 'YTick', sharedTicks);

%% Panel c — continuous 2D heatmap
ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.50]);

% pcolor, not imagesc — each cell needs to be a separate, individually
% editable/movable/recolourable vector shape in Inkscape, which imagesc
% (one flat raster block) can't provide. EdgeColor='flat' (not 'none')
% draws each cell's border in its own fill colour rather than leaving it
% undrawn, closing the visible seam between cells while keeping them as
% genuinely separate shapes — the standard fix for pcolor's seam artifact
% that doesn't sacrifice per-cell editability.
displayMatrix = padarray(tuningSurface, [1 1], 0, 'post');
xEdges = 0.5:(numPosBins + 0.5);

p = pcolor(xEdges, speedEdges, displayMatrix);
set(p, 'EdgeColor', 'flat');
set(ax2, 'YDir', 'normal');
colormap(ax2, parula);
set(ax2, 'CLim', [minVal, maxVal]);

% Colorbar positioned to align with panel b (ax1), not panel c — panel
% b's native y-axis is ΔF/F, the same quantity the colorbar decodes,
% whereas panel c's native y-axis is running speed. Same shared ticks as
% ax1's YTick above, so the two read identically rather than relying on
% two independent auto-tick generators to happen to agree.
cb = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.68, 0.025, 0.24]);
cb.Label.String = '\DeltaF/F [NeuC]';
cb.Ticks = sharedTicks;

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