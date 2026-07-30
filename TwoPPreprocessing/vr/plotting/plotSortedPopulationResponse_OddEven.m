function plotSortedPopulationResponse_OddEven(sessionFileInfo, response, signalToUse, applySmoothing)
%   Plots normalised population response heatmaps (odd vs even laps),
%   sorted by peak response location in odd laps, with population mean
%   traces below each heatmap for alignment checking.
%
% Inputs:
%   sessionFileInfo : struct
%       Metadata and file paths for the session
%
%   response : struct
%       Must include response.lapPositionActivity and response.signalUsed
%
%   applySmoothing : logical (optional)
%       Whether to smooth lapPositionActivity across position bins (default = false)
%
% Example:
%   plotSortedPopulationResponse(sessionFileInfo, response, true);
%
% Aman and Sonali - April 2025
if nargin < 3; signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 4; applySmoothing = true; end
%% Output path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end
filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' signalToUse '_ROISandNonROIS_SortedbyOdd.png']);
%% Extract activity matrix
lapActivity = response.lapPositionActivity.(signalToUse);
% Optional spatial smoothning
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
% Normalize across position bins
normOdd = normalize(meanOdd, 2, 'range');
normEven = normalize(meanEven, 2, 'range');
% Sort cells by peak location in odd lap average
[~, peakIdx] = max(normOdd, [], 2);
[~, sortIdx] = sort(peakIdx);

%% Population mean across ROIs (for alignment check at the bottom)
% Computed on the normalised per-ROI traces so scale matches the heatmap
popMeanOdd = mean(normOdd, 1, 'omitnan');
popMeanEven = mean(normEven, 1, 'omitnan');
popSemOdd  = std(normOdd, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(normOdd), 1));
popSemEven = std(normEven, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(normEven), 1));
numPosBins = size(normOdd, 2);
posVec = 1:numPosBins;

%% Determine smoothing label for figure
% if applySmoothing
%     smoothingLabel = 'smoothed (plot-time)';
% elseif isfield(response, 'smoothingApplied') && response.smoothingApplied
%     smoothingLabel = 'smoothed (precomputed)';
% else
%     smoothingLabel = 'unsmoothed';
% end
%% Plot
figure('Position', [100 100 1100 650]);

% --- Odd laps heatmap ---
ax_heat_odd = subplot(2, 2, 1);
imagesc(normOdd(sortIdx, :));
clim([0 1]); colormap(flipud(gray));
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
xline(40, 'k--', 'LineWidth', 2.5);
xline(80, 'k--', 'LineWidth', 2.5);
xline(120, 'k--', 'LineWidth', 2.5);
xline(160, 'k--', 'LineWidth', 2.5);
xticks([1 40 80 120 160]);
xticklabels({'1', '40', '80', '120', '160'});
ylabel('ROIs');
title([sessionFileInfo.animal_name ' - Odd laps sorted (' signalToUse ')']);
colorbar; ylabel(colorbar, 'Activity (normalised)');

% --- Even laps heatmap ---
ax_heat_even = subplot(2, 2, 2);
imagesc(normEven(sortIdx, :));
clim([0 1]); colormap(flipud(gray));
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
xline(40, 'k--', 'LineWidth', 2.5);
xline(80, 'k--', 'LineWidth', 2.5);
xline(120, 'k--', 'LineWidth', 2.5);
xline(160, 'k--', 'LineWidth', 2.5);
xticks([1 40 80 120 160]);
xticklabels({'1', '40', '80', '120', '160'});
ylabel('ROI');
title([sessionFileInfo.animal_name ' - Even laps sorted (by odd) (' signalToUse ')']);
colorbar; ylabel(colorbar, 'Activity (normalised)');

% --- Odd laps population mean (bottom left) ---
ax_mean_odd = subplot(2, 2, 3);
hold on;
x_patch = [posVec, fliplr(posVec)];
y_patch = [(popMeanOdd + popSemOdd), fliplr(popMeanOdd - popSemOdd)];
nanMask = isnan(x_patch) | isnan(y_patch);
x_patch(nanMask) = []; y_patch(nanMask) = [];
fill(x_patch, y_patch, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(posVec, popMeanOdd, 'k', 'LineWidth', 2);
xline(40, 'k--', 'LineWidth', 1.5);
xline(80, 'k--', 'LineWidth', 1.5);
xline(120, 'k--', 'LineWidth', 1.5);
xline(160, 'k--', 'LineWidth', 1.5);
xticks([1 40 80 120 160]);
xticklabels({'1', '40', '80', '120', '160'});
xlabel('Position (cm)');
ylabel('Pop. mean (norm.)');
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12);
hold off;

% --- Even laps population mean (bottom right) ---
ax_mean_even = subplot(2, 2, 4);
hold on;
x_patch = [posVec, fliplr(posVec)];
y_patch = [(popMeanEven + popSemEven), fliplr(popMeanEven - popSemEven)];
nanMask = isnan(x_patch) | isnan(y_patch);
x_patch(nanMask) = []; y_patch(nanMask) = [];
fill(x_patch, y_patch, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(posVec, popMeanEven, 'k', 'LineWidth', 2);
xline(40, 'k--', 'LineWidth', 1.5);
xline(80, 'k--', 'LineWidth', 1.5);
xline(120, 'k--', 'LineWidth', 1.5);
xline(160, 'k--', 'LineWidth', 1.5);
xticks([1 40 80 120 160]);
xticklabels({'1', '40', '80', '120', '160'});
xlabel('Position (cm)');
ylabel('Pop. mean (norm.)');
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12);
hold off;

% Link x-axes so zooming/panning stays aligned between heatmap and mean trace
linkaxes([ax_heat_odd, ax_mean_odd], 'x');
linkaxes([ax_heat_even, ax_mean_even], 'x');
xlim(ax_mean_odd, [1, numPosBins]);
xlim(ax_mean_even, [1, numPosBins]);

%% Save
set(gcf, 'PaperUnits', 'inches', ...
'PaperPosition', [0 0 11 8.5], ...
'PaperOrientation', 'landscape');
print(gcf, filename, '-dpng', '-r300');
end
