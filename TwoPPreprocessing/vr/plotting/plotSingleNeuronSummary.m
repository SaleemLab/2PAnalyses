function plotSingleNeuronSummary(sessionFileInfo, response, neuronIdx, applySmoothing)
%   Plots a summary of one neuron's activity across laps:
%   - Mean ± SEM across laps
%   - Heatmap of lap-by-position activity
%
% Inputs:
%   sessionFileInfo : struct
%   response : struct
%       Must include:
%           - lapPositionActivity (ROIs x laps x bins)
%           - signalUsed
%   neuronIdx : index of the neuron (row in lapPositionActivity)
%   applySmoothing : logical (optional)

if nargin < 4
    applySmoothing = false;
end

%% Output path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'SingleNeuronSummaries');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

%% Extract activity
lapActivityFull = response.lapPositionActivity;

if neuronIdx < 1 || neuronIdx > size(lapActivityFull, 1)
    error('Invalid neuron index: %d', neuronIdx);
end

roiActivity = squeeze(lapActivityFull(neuronIdx, :, :)); % laps x bins

if all(isnan(roiActivity), 'all')
    warning('Neuron %d has all-NaN activity. Skipping plot.', neuronIdx);
    return;
end

%% Optional smoothing
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

%% Compute mean and SEM
meanActivity = mean(roiActivity, 1, 'omitnan');
semActivity = std(roiActivity, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(roiActivity), 1));

%% Normalize lap heatmap
normLapActivity = normalize(roiActivity, 2, 'range');

%% Plot
figure('Position', [100 100 1200 400]);

% --- Line plot (mean ± SEM) ---
subplot(1, 2, 1);
hold on;
x = 1:size(meanActivity, 2);
fill([x fliplr(x)], ...
     [meanActivity + semActivity, fliplr(meanActivity - semActivity)], ...
     [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, meanActivity, 'k', 'LineWidth', 2);
xline(50, 'k--'); xline(70, 'k--'); xline(90, 'k--'); xline(110, 'k--');
xticks([0 50 70 90 110 140]);
xticklabels({'0', '50', '70', '90', '110', '140'});
xlabel('Position (cm)');
ylabel('Mean activity');
title(sprintf('%s - ROI %d (%s)', ...
    sessionFileInfo.animal_name, neuronIdx, response.signalUsed));
legend({'SEM', 'Mean'}, 'Location', 'northeast');
box off; set(gca, 'FontSize', 12);

% --- Heatmap ---
subplot(1, 2, 2);
imagesc(normLapActivity);
caxis([0 1]); colormap(flipud(gray));
xline(50, 'k--'); xline(70, 'k--'); xline(90, 'k--'); xline(110, 'k--');
xticks([0 50 70 90 110 140]);
xticklabels({'0', '50', '70', '90', '110', '140'});
xlabel('Position (cm)');
ylabel('Lap #');
title('Lap-by-position activity (normalized)');
colorbar; ylabel(colorbar, 'Activity (normalized)');
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');

%% Save figure
figName = sprintf('%s_%s_ROI%d_%s_Summary.png', ...
    sessionFileInfo.animal_name, ...
    sessionFileInfo.session_name, ...
    neuronIdx, ...
    response.signalUsed);
figPath = fullfile(figSaveDir, figName);
set(gcf, 'PaperUnits', 'inches', ...
         'PaperPosition', [0 0 11 4], ...
         'PaperOrientation', 'landscape');
print(gcf, figPath, '-dpng', '-r300');
close(gcf);
end


