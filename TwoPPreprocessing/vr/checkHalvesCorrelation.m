function [r, p, stableIdx] = checkHalvesCorrelation(sessionFileInfo, response, signalToUse, applySpatialSmoothing, plotFlag)
% Calculates the first half vs. second half laps spatial correlation for
% each ROI.
%
%   r         - [nROIs x 1] vector of Pearson correlation coefficients (r-value)
%   p         - [nROIs x 1] vector of p-values for each correlation
%   stableIdx - [nStableROIs x 1] indices of ROIs with r > 0.4
%% Handle optional inputs
if nargin < 3; signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 4; applySpatialSmoothing = true; end
if nargin < 5; plotFlag = true; end

%% Get data
% lapPositionActivity is (ROI x Laps x Position)
lapPositionActivity = response.lapPositionActivity.(signalToUse);

% Get total number of laps
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

%% Optional spatial smoothing before computing correlations
if applySpatialSmoothing
    fprintf('Applying spatial smoothing...\n');
    w = gausswin(5); % 10-bin Gaussian window
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
stableThresh = 0.4;
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
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, ...
       "lapCorr_Halves", ...
         '-append')
         
elseif isfield(sessionFileInfo, 'otherSessFilePaths')
    warning('sessionROIData file not found at: %s. Cannot append peak significance data.', ...
        sessionFileInfo.otherSessFilePaths.sessionROIData);
else
    warning('sessionFileInfo.otherSessFilePaths field not found. Cannot save peak significance data.');
end

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
    caxis([0 1]); colormap(flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
    xline(50, 'k--', 'LineWidth', 1.5);
    xline(70, 'k--', 'LineWidth', 1.5);
    xline(90, 'k--', 'LineWidth', 1.5);
    xline(110, 'k--', 'LineWidth', 1.5);
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    title(sprintf('First Half Laps (n=%d)', size(firstHalfLaps, 2)));
    xlabel('Position (cm)');
    ylabel('Stable ROI (Sorted)');
    
    % Second Half Laps
    subplot(122)
    imagesc(normSecondStable(sortIdx,:));
    caxis([0 1]); colormap(flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
    xline(50, 'k--', 'LineWidth', 1.5);
    xline(70, 'k--', 'LineWidth', 1.5);
    xline(90, 'k--', 'LineWidth', 1.5);
    xline(110, 'k--', 'LineWidth', 1.5);
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    title(sprintf('Second Half Laps (n=%d)', size(secondHalfLaps, 2)));
    xlabel('Position (cm)');
    ylabel('Stable ROIs (Sorted)');
    
    % Save functionality commented out for generality
    % saveas(fig1,'\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M25041\Analysis\20250416\Figures\StableROIs_Halves.png')
end
end