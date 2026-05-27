function [r, p, stableIdx] = checkHalvesCorrelation(sessionFileInfo, response, signalToUse, applySpatialSmoothing, excludeFlaggedLaps, plotFlag)
% Calculates the first half vs. second half laps spatial correlation for
% each ROI.
%
%   r         - [nROIs x 1] vector of Pearson correlation coefficients (r-value)
%   p         - [nROIs x 1] vector of p-values for each correlation
%   stableIdx - [nStableROIs x 1] indices of ROIs with r > 0.4
%% Handle optional inputs
if nargin < 3; signalToUse = 'dFFNeuropilCorrected'; end %changed to dff Jan 2026
if nargin < 4; applySpatialSmoothing = true; end
if nargin < 5; excludeFlaggedLaps = false; end
if nargin < 6; plotFlag = false; end

%% Save figure save path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end
filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' signalToUse '_halvesCorr_SortedbyOdd.png']);

%% Get data
rawActivity = response.lapPositionActivity.(signalToUse);
nTotalLaps = size(rawActivity, 2);

% first optionally remove the flagged trials
if excludeFlaggedLaps && isfield(response, 'flaggedLaps')
    disp('Excluding flagged laps before computing this metric.. ')
    unflaggedMask = true(1, nTotalLaps);
    unflaggedMask(response.flaggedLaps) = false;
else
    unflaggedMask = true(1, nTotalLaps);
end

% then only use the base laps to compute correlations
baseTrialsMask = false(1, nTotalLaps);
baseTrialsMask(response.trialIndicesByCondition.Baseline) = true;

% use the keep lap index and the baseline to find interset and carry on
% with the rest of the aanalysis..
finalValidLapsMask = unflaggedMask & baseTrialsMask;
lapPositionActivity = rawActivity(:, finalValidLapsMask, :);

% lapPositionActivity is (ROI x Laps x Position)
nLaps = size(lapPositionActivity, 2);

% --- Modification for Halves Split ---
% Calculate the split point
midpoint = floor(nLaps / 2);

% Handle the case with very few laps (e.g., less than 2)
if midpoint < 1
    warning('Not enough laps (%d) to split into two halves. Returning NaNs.', nLaps);
    nROIs = size(lapPositionActivity, 1);
    r = NaN(nROIs, 1);
    p = NaN(nROIs, 1);
    stableIdx = [];
    return;
end
% --------------------------------------

%% Optional spatial smoothning before computing correlations
if applySpatialSmoothing
    fprintf('Applying spatial smoothing...\n');
    w = gausswin(10); % 10-bin Gaussian window
    w = w / sum(w);
    for iCell = 1:size(lapPositionActivity, 1)
        for iLap = 1:size(lapPositionActivity, 2)
            trace = squeeze(lapPositionActivity(iCell, iLap, :));
            if all(isnan(trace)), continue; end

            % Handle NaNs for filtering
            nanMask = isnan(trace);
            trace(nanMask) = 0; % Temporarily set NaNs to 0 for filtfilt

            smoothed = filtfilt(w, 1, trace);
            smoothed(nanMask) = NaN; % Restore NaNs
            lapPositionActivity(iCell, iLap, :) = smoothed;
        end
    end
end

%% Split into first and second half laps
firstHalfLaps = lapPositionActivity(:, 1:midpoint, :);
secondHalfLaps = lapPositionActivity(:, (midpoint + 1):end, :);

% Average tuning curves for each half
% Results (meanFirst, meanSecond) are (ROI x Position)
meanFirst = squeeze(mean(firstHalfLaps, 2, 'omitnan'));
meanSecond = squeeze(mean(secondHalfLaps, 2, 'omitnan'));

%% Normalise tuning curves
% This scales each cell's tuning curve to the [0, 1] range.
normFirst = normalize(meanFirst, 2, 'range');
normSecond = normalize(meanSecond, 2, 'range');

%% Calculate Pearson correlation
% Correlate the (Position x ROI) matrices.
% correlation between cell i (first half) and cell j (second half)
[r_matrix, p_matrix] = corr(normFirst', normSecond', 'rows', 'pairwise');
% We only care about the diagonal: corr(cell_i_first, cell_i_second)
r = diag(r_matrix);
p = diag(p_matrix);

%% Identify "stable" cells (0.4 threshold)
stableThresh = 0.6;
stableIdx = find(r > stableThresh);

%% Save results
lapCorr_Halves.rho = r;
lapCorr_Halves.p = p;
lapCorr_Halves.stableThreshold = stableThresh;
lapCorr_Halves.stableIdx = stableIdx;
lapCorr_Halves.firstHalfLapsCount = size(firstHalfLaps, 2);
lapCorr_Halves.secondHalfLapsCount = size(secondHalfLaps, 2);

% Check if the file path exists and save the variables
if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
    disp(['Saving First/Second-Half Lap significance results to: ', sessionFileInfo.otherSessFilePaths.sessionROIData]);
    % Save the new structure name (lapCorr_Halves)
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, "lapCorr_Halves", '-append')
elseif isfield(sessionFileInfo, 'otherSessFilePaths')
    warning('sessionROIData file not found at: %s. Cannot append peak significance data.', ...
        sessionFileInfo.otherSessFilePaths.sessionROIData);
else
    warning('sessionFileInfo.otherSessFilePaths field not found. Cannot save peak significance data.');
end

fprintf('Found %d stable cells (r > %.2f)...\n', length(stableIdx), stableThresh);

%% Plotting
if plotFlag
    fprintf('Plotting %d stable cells (r > %.2f)...\n', length(stableIdx), stableThresh);

    % Get the tuning curves for stable cells
    normFirstStable = normFirst(stableIdx, :);
    normSecondStable = normSecond(stableIdx, :);

    % Sort them by the peak of the first-half tuning curve
    [~, peakIdx] = max(normFirstStable, [], 2);
    [~, sortIdx] = sort(peakIdx);

    fig1 = figure('Name', 'First-Half vs. Second-Half Lap Correlations');

    % First Half Laps
    subplot(121)
    imagesc(normFirstStable(sortIdx,:));
    clim([0 1]); colormap(flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
    xline(40, 'k--', 'LineWidth', 1.5);
    xline(80, 'k--', 'LineWidth', 1.5);
    xline(120, 'k--', 'LineWidth', 1.5);
    xline(160, 'k--', 'LineWidth', 1.5);
    xticks([1 40 80 120 160 200]);
    xticklabels({'1', '40', '80', '120', '160', '200'});
    title(sprintf('First Half Laps (n=%d)', size(firstHalfLaps, 2)));
    xlabel('Position (cm)');
    ylabel('Stable ROI (Sorted)');

    % Second Half Laps
    subplot(122)
    imagesc(normSecondStable(sortIdx,:));
    clim([0 1]); colormap(flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
    xline(40, 'k--', 'LineWidth', 1.5);
    xline(80, 'k--', 'LineWidth', 1.5);
    xline(120, 'k--', 'LineWidth', 1.5);
    xline(160, 'k--', 'LineWidth', 1.5);
    xticks([1 40 80 120 160 200]);
    xticklabels({'1', '40', '80', '120', '160', '200'});
    title(sprintf('Second Half Laps (n=%d)', size(secondHalfLaps, 2)));
    xlabel('Position (cm)');
    ylabel('Stable ROIs (Sorted)');

    %% Save
    set(gcf, 'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 11 8.5], ...
        'PaperOrientation', 'landscape');
    print(gcf, filename, '-dpng', '-r300');
end
end