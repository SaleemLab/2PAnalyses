function plotSortedPopulationResponseAll(sessionFileInfo, response, applySmoothing)
%   Plots normalised population response heatmaps (all laps combined),
%   sorted by peak response location across all laps.
%   Optionally applies temporal smoothing across position bins at plot time.
%
% Inputs:
%   sessionFileInfo : struct
%   response        : struct with response.lapPositionActivity and response.signalUsed
%   applySmoothing  : logical (optional, default = false)

if nargin < 3
    applySmoothing = false;
end

%% Output path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_SortedAcrossAll_deltaFoverF_smoothed' response.signalUsed '.pdf']);

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

%% Combine all laps
meanAll = squeeze(mean(lapActivity, 2, 'omitnan'));

% Normalize across position bins
normAll = normalize(meanAll, 2, 'range');

% Sort cells by peak location in the mean
[~, peakIdx] = max(normAll, [], 2);
[~, sortIdx] = sort(peakIdx);

%% Determine smoothing label for figure
if applySmoothing
    smoothingLabel = 'smoothed';
elseif isfield(response, 'smoothingApplied') && response.smoothingApplied
    smoothingLabel = 'smoothed(precomputed)';
else
    smoothingLabel = 'unsmoothed';
end

%% Plot
figure('Position', [100 100 700 500]);
imagesc(normAll(sortIdx, :));
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
title([sessionFileInfo.animal_name ' - All laps sorted (' response.signalUsed ', ' smoothingLabel ')']);
colorbar; ylabel(colorbar, 'Activity (normalized)');

%% Save
set(gcf, 'PaperUnits', 'inches', ...
         'PaperPosition', [0 0 11 8.5], ...
         'PaperOrientation', 'landscape');
print(gcf, filename, '-dpdf', '-r300');
end
