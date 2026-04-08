function [stratStruct, speedMetrics] = stratifyLapSpeeds(response, thresholds)
% STRATIFYLAPSPEEDS Groups laps into speed bins using a robust trimmean
% thresholds: [low_min, med_min, high_min] e.g., [1, 10, 30]

if nargin < 2
    thresholds = [1, 10, 25]; % Default paper thresholds
end

numLaps = length(response.lapRunningSpeed);
Speeds = zeros(numLaps, 1);

for l = 1:numLaps
    currentLapRaw = response.lapRunningSpeed{l};

    % trimmean(..., 25) ignores the top/bottom 12.5% of the 'up and down'
    if ~isempty(currentLapRaw)
        Speeds(l) = trimmean(currentLapRaw, 25);
    else
        Speeds(l) = 0;
    end
end


stratStruct.lowIdx  = Speeds >= thresholds(1) & Speeds < thresholds(2);
% 
stratStruct.medIdx  = Speeds >= thresholds(2) & Speeds < thresholds(3);
% 
stratStruct.highIdx = Speeds >= thresholds(3);


[speedMetrics.sortedSpeeds, speedMetrics.sortIdx] = sort(Speeds);
speedMetrics.SpeedsPerLap = Speeds;

fprintf('Stratification Results:\n');
fprintf('Low Speed Laps:  %d\n', sum(stratStruct.lowIdx));
fprintf('Med Speed Laps:  %d\n', sum(stratStruct.medIdx));
fprintf('High Speed Laps: %d\n', sum(stratStruct.highIdx));


if iscell(response.lapPositionRunningSpeed)
 
    speedMatrix = vertcat(response.lapPositionRunningSpeed{:});
else
    speedMatrix = response.lapPositionRunningSpeed;
end


figure('Name', 'Stratification: Sorted Speed Heatmaps', 'Position', [100, 100, 800, 900]);
t = tiledlayout(3,1, 'TileSpacing', 'compact');

% Group metadata for the loop
masks = {stratStruct.highIdx, stratStruct.medIdx, stratStruct.lowIdx};
titles = {'High Speed Group (>25 cm/s)', 'Medium Speed Group (10-25 cm/s)', 'Low Speed Group (1-10 cm/s)'};

for i = 1:3
    nexttile;
    currentMask = masks{i};


    speedCats = speedMetrics.SpeedsPerLap(currentMask);
    % Find the sort order within this specific group (slowest at bottom, Fastest at top)
    [~, localSortIdx] = sort(speedCats, 'ascend');

    % Extract the raw 200-bin speed data for this group
    categoryData = speedMatrix(currentMask, :);

    % Apply the local sort to the rows of the category matrix
    sortedCategoryData = categoryData(localSortIdx, :);

    % Plot Speed vs. Position
    imagesc(1:200, 1:sum(currentMask), sortedCategoryData);

    % Formatting
    colormap(jet);
    c = colorbar;
    c.Label.String = 'cm/s';
    clim([0 60]); % scale to see difference between subplots 

    title([titles{i}, ' (n=', num2str(sum(currentMask)), ')']);
    ylabel('Laps (Sorted)');
end

xlabel(t, 'Position Bin (1-200)');
ylabel(t, 'Individual Laps (Within-Category Speed Sort)');
end


% target = 11; % 
% neuronData = squeeze(response.lapPositionActivity.dFF(target, :, :));
% 
% % Calculate the 3 "Paper-Style" captures
% lowCurve  = nanmean(neuronData(stratStruct.lowIdx, :), 1);
% medCurve  = nanmean(neuronData(stratStruct.medIdx, :), 1);
% highCurve = nanmean(neuronData(stratStruct.highIdx, :), 1);
% 
% nLow  = sum(stratStruct.lowIdx);
% nMed  = sum(stratStruct.medIdx);
% nHigh = sum(stratStruct.highIdx);
% 
% % Plot them together
% figure; hold on;
% plot(lowCurve, 'c', 'LineWidth', 2); 
% plot(medCurve, 'k', 'LineWidth', 2); 
% plot(highCurve, 'm', 'LineWidth', 2);
% title(sprintf('Bouton %d Speed-Stratified Tuning (n: Low=%d, Med=%d, High=%d)', ...
%       target, nLow, nMed, nHigh));
% legend('low', 'medium', 'high')
% xline(40, 'k--', 'LineWidth', 2.5);
% xline(80, 'k--', 'LineWidth', 2.5);
% xline(120, 'k--', 'LineWidth', 2.5);
% xline(160, 'k--', 'LineWidth', 2.5);
% 
