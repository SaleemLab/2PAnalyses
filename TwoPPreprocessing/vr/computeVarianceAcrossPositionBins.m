% function [ratioVarToTuningRange, ratioVarToTuningVar, meanTuning]= computeVarianceAcrossPositionBins(sessionFileInfo, response, signalToUse)
% 
% % calculates metrics of spatial selectivity and lap-to-lap reliability 
% % for each ROI by comparing the variance of the average
% % spatial tuning curve against the average activity 
% % variance within individual position bins .
% 
% %% Handle optional arguments 
% if nargin < 3; signalToUse = 'dFFNeuropilCorrected'; end
% 
% %% Get data 
% if response.signalsZScored
%     disp('Using zScored-lapPositionActivity for computations..')
% else
%     disp('Using lapPositionActivity (without zscoring) for computations..')
% end 
% 
% lapPositionActivity = response.lapPositionActivity.(signalToUse);
% % Mean tuning across laps 
% meanTuning   = squeeze(nanmean(lapPositionActivity,2));
% % Measure variance across lap-averaged position bins
% tuningVar = var(meanTuning,0,2);
% 
% % Tuning range 
% % Max - Min Mean Activity. This is the "signal magnitude." 
% % It measures the difference between the peak and trough of the 
% % average spatial tuning curve. It represents the maximum possible 
% % change in the signal due to position.
% tuningRange = (max(meanTuning')-min(meanTuning'))';
% 
% % Measure mean variance across all laps
% % lap-to-lap variability?
% % It measures the average variance of the activity within a single position
% % bin across all laps. High posAllVar means the cell's activity at any 
% % given location is inconsistent from one lap to the next.
% posAllVar = nanmean(nanvar(lapPositionActivity,0,2),3);
% 
% 
% %% How noisy is the roi's activity compared to the total span/magnitude of its signal
% % LOW = the cell's activity across laps is stable and precise relative 
% % to the magnitude of its spatial change. The position-dependent change
% %  (signal) is much larger than the lap-to-lap fluctuation (noise).
% 
% % HIGH = he cell's activity is highly variable and that the lap-to-lap 
% % inconsistency (noise) is a substantial fraction of the cell's total 
% % firing range (signal).
% ratioVarToTuningRange = posAllVar./tuningRange;
% 
% %% tuning quality or spatial selectivity relative to noise. 
% % It compares the inherent noise to the actual structure of the spatial tuning curve.
% % LOW = strongly-tuned neuron. 
% % It means the variance of the average spatial tuning curve 
% % is much larger than the lap-to-lap variability (low noise). 
% % The cell is selective for a position and reliably fires at that position.
% 
% % HIGH = lap-to-lap variability is high.
% % suggests that the lap-to-lap variability is high relative to the 
% % strength of the spatial tuning. The cell might have a weak or
% % poor place field structure that is easily obscured by the noise.
% ratioVarToTuningVar   = posAllVar./tuningVar;
% 
% %% 
% tuningCurveVariance.ratioVarToTuningVar = ratioVarToTuningVar;
% tuningCurveVariance.ratioVarToTuningRange = ratioVarToTuningRange; 
% 
% 
% %% Check if sessionROIData exists and append 
% if exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
%     save(sessionFileInfo.otherSessFilePaths.sessionROIData, ...
%        "tuningCurveVariance", ...
%          '-append')
% else
%     warning('sessionROIData file not found at: %s. Cannot append variance data.', ...
%         sessionFileInfo.otherSessFilePaths.sessionROIData);
% end
% 
% end 
% 
% 
% % threshold = 10;
% % selective_ROI_indices = find(ratioVarToTuningVar <= threshold);
% % selective_tuning_curves = meanTuning(selective_ROI_indices, :);
% % figure;
% % num_selective_rois = length(selective_ROI_indices);
% % for i = 1:num_selective_rois
% %     current_roi_number = selective_ROI_indices(i);
% %     clf; 
% %     current_tuning_curve = meanTuning(current_roi_number, :);
% %     plot(current_tuning_curve, 'LineWidth', 2, 'Color', 'b');
% %     xlabel('Position Bin');
% %     ylabel('Mean Activity');
% %     title_text = ['Selective ROI #', num2str(current_roi_number), ...
% %                   ' (', num2str(i), ' of ', num2str(num_selective_rois), ')'];
% %     title(title_text);
% %     pause(0.5); 
% % 
% % end
function [ratioVarToTuningRange, ratioVarToTuningVar, meanTuning]= computeVarianceAcrossPositionBins(sessionFileInfo, response, signalToUse)
% calculates metrics of spatial selectivity and lap-to-lap reliability 
% for each ROI by comparing the variance of the average spatial tuning curve
% against the average activity variance within individual position bins.
% Calculates both ratios but plots only the Noise/TuningVar ratio (Selectivity).
%% Handle optional arguments 
if nargin < 3; signalToUse = 'dFF'; end

%% Save figure 
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' signalToUse '_TuningVarFigures.png']);
%% Get data 
if response.signalsZScored
    disp('Using zScored-lapPositionActivity for computations..')
else
    disp('Using lapPositionActivity (without zscoring) for computations..')
end 
lapPositionActivity = response.lapPositionActivity.(signalToUse);

% Mean tuning across laps 
meanTuning   = squeeze(nanmean(lapPositionActivity,2));

% Measure variance across lap-averaged position bins (Signal Structure)
% Measures how much the mean tuning curve ITSELF varies across the track.
tuningVar = var(meanTuning,0,2);

% Tuning range (Signal Magnitude)
% Max - Min Mean Activity.
tuningRange = (max(meanTuning')-min(meanTuning'))';

% Measure mean variance across all laps (Noise)
% Measures the average variance of the activity within a single position
% bin across all laps.
posAllVar = nanmean(nanvar(lapPositionActivity,0,2),3);

%% Calculate Both Ratios
% Noise vs. Signal Magnitude (Reliability)
ratioVarToTuningRange = posAllVar./tuningRange;

%Noise vs. Signal Structure Variance (Selectivity)
ratioVarToTuningVar   = posAllVar./tuningVar;

%% Store Both Results
tuningCurveVariance.ratioVarToTuningVar = ratioVarToTuningVar; 
tuningCurveVariance.ratioVarToTuningRange = ratioVarToTuningRange; 

% Filter out NaNs for the ratio we are plotting
validRatioVar   = ratioVarToTuningVar(~isnan(ratioVarToTuningVar) & isfinite(ratioVarToTuningVar));

if isempty(validRatioVar)
    warning('No valid ROIs found for plotting variance ratios.');
    % Skip plotting, but proceed to saving the calculated metrics
else

    figure('Position', [100 100 800 600]); % Adjusted figure width for 2 tiles
    t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact'); % 1 row, 2 columns

    % --- Tile 1: Histogram of RatioVarToTuningVar ---
    ax1 = nexttile;
    histogram(ax1, validRatioVar, 50, 'EdgeColor', 'k', 'FaceColor', [0.9 0.7 0.5]);
    title(ax1, 'Ratio: Noise / Tuning Var');
    xlabel(ax1, 'Var / Tuning Var (\leftarrow Selective)');
    ylabel(ax1, 'Count of ROIs');

    % Corrected xline syntax
    xline(ax1, nanmedian(validRatioVar), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5, 'Label', 'Median'); 

    set(ax1, 'YScale', 'log'); 

    % --- Tile 2: Example Tuning Curves (Low vs. High RatioVarToTuningVar) ---
    ax2 = nexttile;

    % Find an example ROI with a LOW ratio (selective/good structure)
    [~, lowRatioIdx] = min(ratioVarToTuningVar);
    % Find an example ROI with a HIGH ratio (unselective/poor structure)
    [~, highRatioIdx] = max(ratioVarToTuningVar);

    % Plot the mean tuning curves
    plot(ax2, meanTuning(lowRatioIdx, :), 'Color', [0 0.5 0], 'LineWidth', 2); hold on;
    plot(ax2, meanTuning(highRatioIdx, :), 'Color', [0.5 0 0], 'LineWidth', 2); hold off;

    % Add legend and labels
    legend(ax2, ...
        sprintf('Low Ratio (Idx %d)', lowRatioIdx), ...
        sprintf('High Ratio (Idx %d)', highRatioIdx), ...
        'Location', 'best', 'Interpreter', 'none');
    title(ax2, 'Example Tuning Curves (Based on Var/TuningVar)');
    xlabel(ax2, 'Position Bins');
    ylabel(ax2, sprintf('Mean %s Activity', signalToUse));

    % Add global title
    sessionTitle = sprintf('%s - %s', sessionFileInfo.animal_name, sessionFileInfo.session_name);
    title(t, {'Tuning Selectivity Metric: Noise/TuningVar Ratio', sessionTitle}, 'Interpreter', 'none');

    % Optional: Save the figure


    %% Save
    set(gcf, 'PaperUnits', 'inches', ...
             'PaperPosition', [0 0 11 8.5], ...
             'PaperOrientation', 'landscape');
    print(gcf, filename, '-dpng', '-r300');

end

%% Check if sessionROIData exists and append 
if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, ...
       "tuningCurveVariance", ...
         '-append')
else
    warning('sessionROIData file not found at: %s. Cannot append variance data.', ...
        sessionFileInfo.otherSessFilePaths.sessionROIData);
end
end