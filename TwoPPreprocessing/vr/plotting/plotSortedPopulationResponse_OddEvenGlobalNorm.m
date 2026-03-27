function plotSortedPopulationResponse_OddEvenGlobalNorm(sessionFileInfo, response, signalToUse, applySmoothing)
% plotSortedPopulationResponse_GlobalNorm: Population heatmaps with linked scaling.
% Sorts by Odd laps, but normalizes relative to the max of BOTH Odd and Even.

if nargin < 3; signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 4; applySmoothing = true; end

%% Output path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' signalToUse '_ROISandNonROIS_SortedbyOdd_GlobalNorm.png']);

%% Extract activity matrix
lapActivity = response.lapPositionActivity.(signalToUse);

% Optional spatial smoothing
if applySmoothing
    w = gausswin(5); w = w / sum(w);
    for iCell = 1:size(lapActivity, 1)
        for iLap = 1:size(lapActivity, 2)
            trace = squeeze(lapActivity(iCell, iLap, :));
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace);
            trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace);
            smoothed(nanMask) = NaN;
            lapActivity(iCell, iLap, :) = smoothed;
        end
    end
end

%% Split odd and even laps
oddLaps = lapActivity(:, 1:2:end, :);
evenLaps = lapActivity(:, 2:2:end, :);

% Average across laps
meanOdd = squeeze(mean(oddLaps, 2, 'omitnan'));
meanEven = squeeze(mean(evenLaps, 2, 'omitnan'));

%% Global Normalization across position bins
% Concatenate to find the global min/max for each ROI across both conditions
combinedActivity = [meanOdd, meanEven];
globalMin = min(combinedActivity, [], 2);
globalMax = max(combinedActivity, [], 2);
rangeVal = globalMax - globalMin;
rangeVal(rangeVal == 0) = 1; % Avoid division by zero

normOdd = (meanOdd - globalMin) ./ rangeVal;
normEven = (meanEven - globalMin) ./ rangeVal;

% Sort cells by peak location in odd lap average
[~, peakIdx] = max(normOdd, [], 2);
[~, sortIdx] = sort(peakIdx);

%% Plot
figure('Position', [100 100 1100 500]);

% --- Odd laps ---
subplot(1, 2, 1);
imagesc(normOdd(sortIdx, :));
clim([0 1]); colormap(flipud(gray));
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
xline(40, 'k--', 'LineWidth', 2.5);
xline(80, 'k--', 'LineWidth', 2.5);
xline(120, 'k--', 'LineWidth', 2.5);
xline(160, 'k--', 'LineWidth', 2.5);
xticks([1 40 80 120 160]);
xticklabels({'1', '40', '80', '120', '160', '200'});
xlabel('Position (cm)');
ylabel('ROIs');
title([sessionFileInfo.animal_name ' - Odd laps sorted (' signalToUse ')']);
colorbar; ylabel(colorbar, 'Activity (Global Norm)');

% --- Even laps ---
subplot(1, 2, 2);
imagesc(normEven(sortIdx, :));
clim([0 1]); colormap(flipud(gray));
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
xline(40, 'k--', 'LineWidth', 2.5);
xline(80, 'k--', 'LineWidth', 2.5);
xline(120, 'k--', 'LineWidth', 2.5);
xline(160, 'k--', 'LineWidth', 2.5);
xticks([1 40 80 120 160]);
xticklabels({'1', '40', '80', '120', '160', '200'});
xlabel('Position (cm)');
ylabel('ROI');
title([sessionFileInfo.animal_name ' - Even laps sorted (by odd)']);
colorbar; ylabel(colorbar, 'Activity (Global Norm)');

%% Save
set(gcf, 'PaperUnits', 'inches', ...
         'PaperPosition', [0 0 11 8.5], ...
         'PaperOrientation', 'landscape');
print(gcf, filename, '-dpng', '-r300');
end