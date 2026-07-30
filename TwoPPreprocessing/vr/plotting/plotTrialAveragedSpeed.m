function figA = plotTrialAveragedSpeed(response, binSize, landmarkPositions, landmarkWidth)
% Plots trial-averaged running speed vs. position, binned into fixed-width
% position bins (binSize cm), with SEM error bars and landmark bands.
%
% binSize: width of position bins in cm (e.g. 4, 8, 10)
% landmarkPositions: vector of corridor positions (cm) to mark, e.g. [40 80 120 160]
% landmarkWidth: width in cm of shaded landmark bands (default 8)

if nargin < 2, binSize = 4; end
if nargin < 3, landmarkPositions = [40 80 120 160]; end
if nargin < 4, landmarkWidth = 8; end

%% Interpolate lap speeds onto 1cm bins first (no smoothing needed, binning will average)
numLaps    = length(response.lapRunningSpeed);
numPosBins = 200;
allInterpSpeeds = nan(numLaps, numPosBins);

for l = 1:numLaps
    rawV = response.lapRunningSpeed{l};
    if length(rawV) < 5, continue; end
    binnedV = interp1(linspace(1, numPosBins, length(rawV)), rawV, 1:numPosBins, 'linear', 'extrap');
    allInterpSpeeds(l, :) = binnedV;
end

%% Bin into fixed-width position bins
binEdges   = 1:binSize:numPosBins+1;
numBins    = length(binEdges) - 1;
binCentres = binEdges(1:end-1) + binSize/2;

meanSpeedBinned = nan(1, numBins);
semSpeedBinned  = nan(1, numBins);

for b = 1:numBins
    idx = binEdges(b):(binEdges(b+1)-1);
    idx = idx(idx <= numPosBins);
    binVals = allInterpSpeeds(:, idx);
    binVals = mean(binVals, 2, 'omitnan');  % average across positions within bin, per lap
    validLaps = ~isnan(binVals);
    meanSpeedBinned(b) = mean(binVals(validLaps), 'omitnan');
    semSpeedBinned(b)  = std(binVals(validLaps), 'omitnan') / sqrt(max(sum(validLaps), 1));
end

%% Axis limits
minY = 0;
dataMaxY = max(meanSpeedBinned + semSpeedBinned, [], 'omitnan');
if isnan(dataMaxY), dataMaxY = 50; end
maxY = ceil(dataMaxY / 10) * 10;
yLim = maxY + 10;

%% Plot
figA = figure('Name', 'Trial-Averaged Speed (Binned)', ...
    'Position', [100, 100, 650, 520], 'Color', 'w', 'Visible', 'off');
ax = axes('Position', [0.15, 0.15, 0.75, 0.75]);
hold on;

% Landmark shaded bands (drawn first, behind everything)
for lp = landmarkPositions
    xPatch = [lp - landmarkWidth/2, lp + landmarkWidth/2, ...
              lp + landmarkWidth/2, lp - landmarkWidth/2];
    yPatch = [minY, minY, yLim, yLim];
    fill(ax, xPatch, yPatch, [0.3 0.3 0.3], ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');
end

% Binned mean speed with error bars
errorbar(ax, binCentres, meanSpeedBinned, semSpeedBinned, ...
    'Color', [0.8, 0.0, 0.4], 'LineWidth', 2, 'CapSize', 0, ...
    'Marker', 'o', 'MarkerFaceColor', [0.8, 0.0, 0.4], 'MarkerSize', 5);

%% Formatting
ylabel(ax, 'Speed (cm/s)', 'FontSize', 14);
xlabel(ax, 'Position (cm)', 'FontSize', 14);
box off;
set(ax, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 12);
xlim([1, numPosBins]);
xticks([1, 40, 80, 120, 160, 200]);
ylim([minY, yLim]);
yticks(unique([0, 20, 40, maxY]));

if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax, 0); end
if exist('offsetAxes',            'file') == 2, offsetAxes(ax); end

set(figA, 'Visible', 'on');

end