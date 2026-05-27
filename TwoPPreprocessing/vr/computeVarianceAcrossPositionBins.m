function [ratioVarToTuningRange, ratioVarToTuningVar, meanTuning]= computeVarianceAcrossPositionBins(sessionFileInfo, response, signalToUse, plotFlag)
% calculates metrics of spatial selectivity and lap-to-lap reliability 
% for each ROI by comparing the variance of the average spatial tuning curve
% against the average activity variance within individual position bins.
% Calculates both ratios but plots only the Noise/TuningVar ratio (Selectivity).
%% Handle optional arguments 
if nargin < 3; signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 4; plotFlag = false; end 

%% Save figure 
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' signalToUse '_TuningVarFigures.png']);
%% Get data 
if response.lapPositionActivityZScored
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
if plotFlag
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