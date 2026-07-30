function figA = plotSpeedTrajectoriesSimple(response, smoothSigma)
% Loads lap running speed data and plots traces across position bins.
% No colour patches or thresholds — just the raw lap traces with the
% across-lap mean superimposed on top. X-axis uses fixed ticks every 10 cm.

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

meanSpeed = mean(allInterpSpeeds, 1, 'omitnan');

%% Axis limits
minY = 0;
dataMaxY = max(allInterpSpeeds(:));
if isnan(dataMaxY), dataMaxY = 50; end
maxY = ceil(dataMaxY / 10) * 10;  % round up to nearest 10, e.g. 63 -> 70
yLim = maxY + 10;                 % extra headroom so traces aren't clipped

%% Plot
figA = figure('Name', 'Speed Trajectories', ...
    'Position', [100, 100, 650, 520], 'Color', 'w', 'Visible', 'off');
ax = axes('Position', [0.15, 0.15, 0.75, 0.75]);
hold on;

% Speed traces
for l = 1:numLaps
    if all(isnan(allInterpSpeeds(l,:))), continue; end
    plot(ax, 1:numPosBins, allInterpSpeeds(l,:), 'Color', [0.2, 0.2, 0.2, 0.25], 'LineWidth', 1.5);
end

% Superimposed mean speed trace (drawn on top of individual laps)
plot(ax, 1:numPosBins, meanSpeed, 'Color', [0.8, 0.0, 0.4], 'LineWidth', 2.5);

%% Formatting
ylabel(ax, 'Speed (cm/s)', 'FontSize', 14);
xlabel(ax, 'Position (cm)',  'FontSize', 14);
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