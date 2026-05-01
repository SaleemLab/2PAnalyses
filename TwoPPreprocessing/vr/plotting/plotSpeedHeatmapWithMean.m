function plotSpeedHeatmapWithMean(sessionFileInfo, VRStimName, response)
% plotSpeedHeatmapWithMean: Shorter plot format with thicker individual traces.

% 1. Data Prep
binnedSpeed = response.lapPositionRunningSpeed;
[nLaps, nBins] = size(binnedSpeed);
meanSpeed = mean(binnedSpeed, 1, 'omitnan');
posBins = 1:nBins;

% 2. Figure Setup: Reduced height (600 -> 400)
fig = figure('Color', 'w', 'Position', [100 100 700 300]);
t = tiledlayout(1, 1, 'Padding', 'compact');
ax = nexttile; hold on;

% 3. Plot Individual Laps (Light Gray, Thicker)
% Increased LineWidth to 1.2 for better visibility
for iL = 1:nLaps
    plot(posBins, binnedSpeed(iL, :), 'Color', [0.8 0.8 0.8 0.6], 'LineWidth', 1.2);
end

% 4. Plot Population Mean (Bold Black)
plot(posBins, meanSpeed, 'k', 'LineWidth', 3);

% 5. Formatting & Landmarks
tickLocs = [40 80 120 160];
set(gca, 'Box', 'off', 'XTick', tickLocs, 'TickDir', 'out');
xlabel('Position (cm)');
ylabel('Running Speed (cm/s)');
title(sprintf('%s | Individual Laps vs Mean', VRStimName));
ylim([0 max(meanSpeed, [], 'omitnan') * 1.5]);
xlim([1 nBins]);

% Add landmark lines
for p = tickLocs
    xline(p, 'k:', 'Alpha', 0.5);
end

% 6. Save
saveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
saveas(fig, fullfile(saveDir, sprintf('%s_SpeedBehavior_Overlay.png', VRStimName)));
end