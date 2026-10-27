function figA = plotSpeedTrajectories_positionMedian(response, smoothSigma)
% Loads lap running speed data and plots traces across position bins.
% Draws the per-position occupancy-based median speed line (the same
% medianThreshLine computed and saved by getLowHighSpeedPositionMatrix),
% and shades the region above it pink and below it green, directly
% visualizing the low/high speed split used throughout the analysis.
if nargin < 2, smoothSigma = 2; end

%% Interpolate and smooth lap speeds
numLaps    = length(response.lapRunningSpeed);
numPosBins = 200;
allInterpSpeeds = nan(numLaps, numPosBins);
for l = 1:numLaps
    rawV = response.lapRunningSpeed{l};
    if length(rawV) < 5, continue; end
    binnedV = interp1(linspace(1, numPosBins, length(rawV)), rawV, 1:numPosBins, 'linear', 'extrap');
    allInterpSpeeds(l, :) = smoothdata(binnedV, 'gaussian', smoothSigma * 3);
end

%% Get the per-position median speed line
% Prefer the saved medianThreshLine from getLowHighSpeedPositionMatrix so
% this plot matches the actual analysis exactly, not a recomputed
% approximation. Falls back to computing it here (same convention:
% exclude stationary frames <=1 cm/s, median per bin across laps) if it
% isn't present in response.
if isfield(response, 'speedPositionActivity') && isfield(response.speedPositionActivity, 'lowHigh') ...
        && isfield(response.speedPositionActivity.lowHigh, 'medianThreshLine')
    medianLine = response.speedPositionActivity.lowHigh.medianThreshLine;
    fprintf('Using saved medianThreshLine from getLowHighSpeedPositionMatrix\n');
else
    warning('No saved medianThreshLine found in response — computing directly from lap speeds.');
    medianLine = nan(1, numPosBins);
    for b = 1:numPosBins
        binSpeeds = allInterpSpeeds(:, b);
        validBin  = binSpeeds(binSpeeds > 1 & ~isnan(binSpeeds));
        if ~isempty(validBin)
            medianLine(b) = median(validBin);
        end
    end
end

% Fill any NaN bins (insufficient laps at that position) by linear
% interpolation across neighboring bins, purely for a continuous shaded
% region — doesn't affect the underlying analysis, just the visualization.
if any(isnan(medianLine))
    validIdx = ~isnan(medianLine);
    medianLine = interp1(find(validIdx), medianLine(validIdx), 1:numPosBins, 'linear', 'extrap');
end

%% Axis limits
minY = 0;
maxY = max(allInterpSpeeds(:)) + 1;
if isnan(maxY), maxY = 50; end

%% Plot
figA = figure('Name', 'Speed Trajectories', ...
    'Position', [100, 100, 650, 520], 'Color', 'w', 'Visible', 'off');
ax = axes('Position', [0.15, 0.15, 0.75, 0.75]);
hold on;

% Shading relative to the (position-varying) median line: pink above,
% green below. Polygon vertices trace along the median line, then back
% along the axis limit, to form a closed region for patch().
xVals = 1:numPosBins;
patch(ax, [xVals, fliplr(xVals)], [medianLine, maxY*ones(1, numPosBins)], ...
    [0.98 0.85 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.6); % pink, above median
patch(ax, [xVals, fliplr(xVals)], [minY*ones(1, numPosBins), fliplr(medianLine)], ...
    [0.85 0.95 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6); % green, below median

% Speed traces
for l = 1:numLaps
    if all(isnan(allInterpSpeeds(l,:))), continue; end
    plot(ax, 1:numPosBins, allInterpSpeeds(l,:), 'Color', [0.2, 0.2, 0.2, 0.25], 'LineWidth', 1.5);
end

% Median line itself, drawn on top of the traces and shading
plot(ax, xVals, medianLine, 'k-', 'LineWidth', 2);

%% Formatting
ylabel(ax, 'Speed (cm/s)', 'FontSize', 14);
xlabel(ax, 'Position (cm)',  'FontSize', 14);
box off;
set(ax, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 12);
xlim([1, numPosBins]);
xticks([1, 40, 80, 120, 160, 200]);
xticklabels({'1', '40', '80', '120', '160', '200'});
ylim([minY, maxY]);
yticks(unique([1, round(min(medianLine)), round(max(medianLine)), round(maxY)]));
if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax, 0); end
if exist('offsetAxes',            'file') == 2, offsetAxes(ax); end
set(figA, 'Visible', 'on');
end