function plotSortedPopulationResponse_OddEven(sessionFileInfo, response, applySmoothing)

%   Plots normalised population response heatmaps (odd vs even laps),
%   sorted by peak response location in odd laps.
%   Optionally applies temporal smoothing across position bins at plot time.
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

if nargin < 3
    applySmoothing = false;
end

%% Output path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name response.stimName 'ROISandNonROIS_SortedbyOdd_' response.signalUsed '.png']);

%% Extract activity matrix
lapActivity = response.lapPositionActivity;

% Optional plot-time smoothing
if applySmoothing
    w = gausswin(9); w = w / sum(w);

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
meanOdd = squeeze(median(oddLaps, 2, 'omitnan'));
meanEven = squeeze(median(evenLaps, 2, 'omitnan'));

% Normalize across position bins
normOdd = normalize(meanOdd, 2, 'range');
normEven = normalize(meanEven, 2, 'range');

% Sort cells by peak location in odd lap average
[~, peakIdx] = max(normOdd, [], 2);
[~, sortIdx] = sort(peakIdx);

%% Determine smoothing label for figure
if applySmoothing
    smoothingLabel = 'smoothed (plot-time)';
elseif isfield(response, 'smoothingApplied') && response.smoothingApplied
    smoothingLabel = 'smoothed (precomputed)';
else
    smoothingLabel = 'unsmoothed';
end

%% Plot
figure('Position', [100 100 1100 500]);

% --- Odd laps ---
subplot(1, 2, 1);
imagesc(normOdd(sortIdx, :));
caxis([0 1]); colormap(flipud(gray));
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
xline(50, 'k--', 'LineWidth', 1.5);
xline(70, 'k--', 'LineWidth', 1.5);
xline(90, 'k--', 'LineWidth', 1.5);
xline(110, 'k--', 'LineWidth', 1.5);
xticks([0 50 70 90 110 140]);
xticklabels({'0', '50', '70', '90', '110', '140'});
xlabel('Position (cm)');
ylabel('ROIs');
title([sessionFileInfo.animal_name ' - Odd laps sorted (' response.signalUsed ', ' smoothingLabel ')']);
colorbar; ylabel(colorbar, 'Activity (normalized)');

% --- Even laps ---
subplot(1, 2, 2);
imagesc(normEven(sortIdx, :));
caxis([0 1]); colormap(flipud(gray));
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
xline(50, 'k--', 'LineWidth', 2.5);
xline(70, 'k--', 'LineWidth', 2.5);
xline(90, 'k--', 'LineWidth', 2.5);
xline(110, 'k--', 'LineWidth', 2.5);
xticks([0 50 70 90 110 140]);
xticklabels({'0', '50', '70', '90', '110', '140'});
xlabel('Position (cm)');
ylabel('ROI');
title([sessionFileInfo.animal_name ' - Even laps sorted (by odd) (' response.signalUsed ', ' smoothingLabel ')']);
colorbar; ylabel(colorbar, 'Activity (normalized)');

%% Save
set(gcf, 'PaperUnits', 'inches', ...
         'PaperPosition', [0 0 11 8.5], ...
         'PaperOrientation', 'landscape');
print(gcf, filename, '-dpng', '-r300');
end
