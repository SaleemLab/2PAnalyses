function [r, p, stableIdx] = checkOddEvenCorrelation(sessionFileInfo, response, signalToUse, applySpatialSmoothing, excludeFlaggedLaps, plotFlag)
% Calculates the odd vs. even laps spatial correlation for roi;
% Use only the baseline condition laps;
%
%   r         - [nROIs x 1] vector of Pearson correlation coefficients (r-value)
%   p         - [nROIs x 1] vector of p-values for each correlation
%   stableIdx - [nStableROIs x 1] indices of ROIs with r > 0.5

%% Handle optional inputs
if nargin < 3; signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 4; applySpatialSmoothing = true; end
if nargin < 5; excludeFlaggedLaps = false; end
if nargin < 6; plotFlag = true; end

stableThresh = 0.6;

figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end
filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' signalToUse '_oddEvenCorr_SortedbyOdd.png']);

rawActivity = response.lapPositionActivity.(signalToUse);
nTotalLaps = size(rawActivity, 2);

if excludeFlaggedLaps && isfield(response, 'flaggedLaps')
    disp('Excluding flagged laps before computing this metric.. ')
    unflaggedMask = true(1, nTotalLaps);
    unflaggedMask(response.flaggedLaps) = false;
else
    unflaggedMask = true(1, nTotalLaps);
end

baseTrialsMask = false(1, nTotalLaps);
baseTrialsMask(response.trialIndicesByCondition.Baseline) = true;

finalValidLapsMask = unflaggedMask & baseTrialsMask;
lapPositionActivity = rawActivity(:, finalValidLapsMask, :);

if applySpatialSmoothing
    fprintf('Applying spatial smoothing...\n');
    w = gausswin(10);
    w = w / sum(w);
    for iCell = 1:size(lapPositionActivity, 1)
        for iLap = 1:size(lapPositionActivity, 2)
            trace = squeeze(lapPositionActivity(iCell, iLap, :));
            if all(isnan(trace)), continue; end

            nanMask = isnan(trace);
            trace(nanMask) = 0;

            smoothed = filtfilt(w, 1, trace);
            smoothed(nanMask) = NaN;
            lapPositionActivity(iCell, iLap, :) = smoothed;
        end
    end
end

oddLaps = lapPositionActivity(:, 1:2:end, :);
evenLaps = lapPositionActivity(:, 2:2:end, :);

meanOdd = squeeze(mean(oddLaps, 2, 'omitnan'));
meanEven = squeeze(mean(evenLaps, 2, 'omitnan'));

normOdd = normalize(meanOdd, 2, 'range');
normEven = normalize(meanEven, 2, 'range');

[r_matrix, p_matrix] = corr(normOdd', normEven', 'rows', 'pairwise');
r = diag(r_matrix);
p = diag(p_matrix);


stableIdx = find(r > stableThresh);

lapCorr_OddEven.rho = r;
lapCorr_OddEven.p = p;
lapCorr_OddEven.stableThreshold = stableThresh;
lapCorr_OddEven.stableIdx = stableIdx;

if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
    disp(['Saving Odd-Even Lap significance results to: ', sessionFileInfo.otherSessFilePaths.sessionROIData]);
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, "lapCorr_OddEven", '-append')
elseif isfield(sessionFileInfo, 'otherSessFilePaths')
    warning('sessionROIData file not found at: %s. Cannot append peak significance data.', ...
        sessionFileInfo.otherSessFilePaths.sessionROIData);
else
    warning('sessionFileInfo.otherSessFilePaths field not found. Cannot save peak significance data.');
end

if plotFlag
    fprintf('Plotting %d stable cells (r > %.2f)...\n', length(stableIdx), stableThresh);

    normOddStable = normOdd(stableIdx, :);
    normEvenStable = normEven(stableIdx, :);

    [~, peakIdx] = max(normOddStable, [], 2);
    [~, sortIdx] = sort(peakIdx);

    fig1 = figure('Name', 'Odd-Even Lap Correlations');

    subplot(121)
    imagesc(normOddStable(sortIdx,:));
    clim([0 1]); colormap(flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
    xline(40, 'k--', 'LineWidth', 1.5);
    xline(80, 'k--', 'LineWidth', 1.5);
    xline(120, 'k--', 'LineWidth', 1.5);
    xline(160, 'k--', 'LineWidth', 1.5);
    xticks([1 40 80 120 160 200]);
    xticklabels({'1', '40', '80', '120', '160', '200'});
    title(sprintf('Odd Laps (n=%d)', size(oddLaps, 2)));
    xlabel('Position (cm)');
    ylabel('Stable ROI (Sorted)');

    subplot(122)
    imagesc(normEvenStable(sortIdx,:));
    clim([0 1]); colormap(flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
    xline(40, 'k--', 'LineWidth', 1.5);
    xline(80, 'k--', 'LineWidth', 1.5);
    xline(120, 'k--', 'LineWidth', 1.5);
    xline(160, 'k--', 'LineWidth', 1.5);
    xticks([1 40 80 120 160 200]);
    xticklabels({'1', '40', '80', '120', '160', '200'});
    title(sprintf('Even Laps (n=%d)', size(evenLaps, 2)));
    xlabel('Position (cm)');
    ylabel('Stable ROIs (Sorted)');

    set(gcf, 'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 11 8.5], ...
        'PaperOrientation', 'landscape');
    print(gcf, filename, '-dpng', '-r300');
end
end