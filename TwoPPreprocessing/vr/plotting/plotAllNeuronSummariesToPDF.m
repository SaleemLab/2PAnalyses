function plotAllNeuronSummariesToPDF(sessionFileInfo, response, applySmoothing)
%   Plots summary for all ROIs (cells and non-cells) and exports all
%   plots into one multipage PDF (one page per ROI).
%
%   Each page includes:
%     - Mean ± SEM trace across laps
%     - Lap-by-position heatmap (normalized per lap)
%
% Inputs:
%   sessionFileInfo : struct
%   response : struct
%       Must include:
%           - lapPositionActivity (ROIs x laps x bins)
%           - signalUsed
%   applySmoothing : logical

if nargin < 3
    applySmoothing = false;
end

%% Output path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

pdfName = sprintf('%s_%s_AllROIs_%s_%s_MeansAndLapPositionHeatMap.pdf', ...
    sessionFileInfo.animal_name, ...
    sessionFileInfo.session_name, ...
    response.stimName,...
    response.signalUsed);
pdfPath = fullfile(figSaveDir, pdfName);

%% Extract activity
lapActivityFull = response.lapPositionActivity;
nROIs = size(lapActivityFull, 1);

fprintf('Generating ROI summary plots for %d ROIs...\n', nROIs);

%% Loop through all ROIs
for neuronIdx = 1:nROIs
    roiActivity = squeeze(lapActivityFull(neuronIdx, :, :)); % laps x bins

    if all(isnan(roiActivity), 'all')
        continue;
    end

    % Optional smoothing
    if applySmoothing
        w = gausswin(9); w = w / sum(w);
        for iLap = 1:size(roiActivity, 1)
            trace = roiActivity(iLap, :);
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace);
            trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace);
            smoothed(nanMask) = NaN;
            roiActivity(iLap, :) = smoothed;
        end
    end

    % Compute mean ± SEM
    meanActivity = mean(roiActivity, 1, 'omitnan');
    semActivity = std(roiActivity, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(roiActivity), 1));

    % Normalize lap heatmap
    normLapActivity = normalize(roiActivity, 2, 'range');

    % Create figure
    fig = figure('Visible', 'off', 'Position', [100 100 1200 400]);

    % --- Line plot ---
    subplot(1, 2, 1); hold on;
    x = 1:size(meanActivity, 2);
    fill([x fliplr(x)], ...
         [meanActivity + semActivity, fliplr(meanActivity - semActivity)], ...
         [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(x, meanActivity, 'k', 'LineWidth', 2);
    xline(50, 'k--'); xline(70, 'k--'); xline(90, 'k--'); xline(110, 'k--');
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    xlabel('Position (cm)'); ylabel(sprintf('Mean %s', response.signalUsed));
    title(sprintf('%s - ROI %d (%s)', ...
        sessionFileInfo.animal_name, neuronIdx, response.signalUsed));
    box off; set(gca, 'FontSize', 12);

    % --- Heatmap ---
    subplot(1, 2, 2);
    imagesc(normLapActivity);
    caxis([0 1]); colormap(flipud(gray));
    xline(50, 'k--'); xline(70, 'k--'); xline(90, 'k--'); xline(110, 'k--');
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    xlabel('Position (cm)'); ylabel('Lap #');
    title('Lap-by-position activity (normalised)');
    colorbar; ylabel(colorbar, sprintf('Activity (%s)', response.signalUsed));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');

    % Save current page to multi-page PDF
    exportgraphics(fig, pdfPath, 'Append', true);
    close(fig);
end

fprintf('Done. All ROI plots saved to:\n%s\n', pdfPath);
end
