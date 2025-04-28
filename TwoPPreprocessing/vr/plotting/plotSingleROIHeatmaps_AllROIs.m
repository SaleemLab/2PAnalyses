function plotSingleROIHeatmaps_AllROIs(sessionFileInfo, response, applySmoothing)
%   Plots ROI activity per lap using heatmaps (lap x position bin).
%   One figure per ROI (cells and non-cells), based on response.lapPositionActivity.
%
% Inputs:
%   sessionFileInfo : struct with session and save info
%   response : struct
%       Must include:
%           - lapPositionActivity (ROIs x laps x position bins)
%           - signalUsed (string for labeling)
%   applySmoothing : logical (optional)

if nargin < 3
    applySmoothing = false;
end

%% Output path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'PerROIHeatmaps_AllROIs');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

%% Extract activity
lapActivityFull = response.lapPositionActivity;
nROIs = size(lapActivityFull, 1);

%% Optional smoothing
if applySmoothing
    w = gausswin(9); w = w / sum(w);
    for iCell = 1:nROIs
        for iLap = 1:size(lapActivityFull, 2)
            trace = squeeze(lapActivityFull(iCell, iLap, :));
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace);
            trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace);
            smoothed(nanMask) = NaN;
            lapActivityFull(iCell, iLap, :) = smoothed;
        end
    end
end

%% Loop through each ROI
for roiIdx = 1:nROIs
    roiActivity = squeeze(lapActivityFull(roiIdx, :, :)); % laps x bins

    % Skip if all NaN
    if all(isnan(roiActivity), 'all')
        continue;
    end

    % Normalize per lap
    normActivity = normalize(roiActivity, 2, 'range');

    % Plot
    figure('Position', [100 100 800 300]);
    imagesc(normActivity);
    caxis([0 1]); colormap(flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
    xline(50, 'k--', 'LineWidth', 1.5);
    xline(70, 'k--', 'LineWidth', 1.5);
    xline(90, 'k--', 'LineWidth', 1.5);
    xline(110, 'k--', 'LineWidth', 1.5);
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    xlabel('Position (cm)');
    ylabel('Lap #');
    title(sprintf('%s - ROI %d (%s)', ...
        sessionFileInfo.animal_name, roiIdx, response.signalUsed));
    colorbar; ylabel(colorbar, 'Activity (normalised)');

    % Save figure
    figName = sprintf('%s_%s_ROI%d_%s.pdf', ...
        sessionFileInfo.animal_name, ...
        sessionFileInfo.session_name, ...
        roiIdx, ...
        response.signalUsed);
    figPath = fullfile(figSaveDir, figName);
    set(gcf, 'PaperUnits', 'inches', ...
             'PaperPosition', [0 0 11 4], ...
             'PaperOrientation', 'landscape');
    print(gcf, figPath, '-dpdf', '-r300');
    close(gcf);
end
end
