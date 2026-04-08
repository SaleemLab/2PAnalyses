function [tuningCurveVariance, sigIdx] = computeVarianceAcrossPositionBins_withShuffle(sessionFileInfo, response, signalToUse, zThreshold)
% Calculates Shuffled-Corrected Variance Ratios and plots 
% Odd vs. Even Heatmaps for Significant Boutons.
%
%   zThreshold: Default is 2 (Structure must be >2SD better than shuffle)
%   Aman and Sonali March 2026

%% 1. Setup and Data Extraction
if nargin < 3 || isempty(signalToUse); signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 4 || isempty(zThreshold); zThreshold = 3; end

% Get Real Activity and the Shuffle Matrix
lapActivity = response.lapPositionActivity.(signalToUse);
shuffMatrix = response.lapPositionActivity_ShuffleMatrix.(signalToUse); 
[numROIs, numLaps, numBins] = size(lapActivity);

% Define Save Path for Diagnostic Figure
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'); mkdir(figSaveDir); end
filename = fullfile(figSaveDir, [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' signalToUse '_VarianceShuffle_SanityCheck.png']);

%% 2. Calculate Real Variance Metrics
% Mean tuning across all laps
meanTuning = squeeze(nanmean(lapActivity, 2));

% Signal Structure: Variance of the mean tuning curve across position bins
tuningVar = var(meanTuning, 0, 2);

% Noise: Average lap-to-lap variance within individual position bins
posAllVar = nanmean(nanvar(lapActivity, 0, 2), 3);

% Real Selectivity Ratio (Lower = More Spatial/Structured)
ratioVarToTuningVar = posAllVar ./ (tuningVar + eps);

%% 3. Calculate Shuffled Null Distribution (Bootstrapping)
numShuffles = size(shuffMatrix, 3);
shuffRatioDist = nan(numROIs, numShuffles);

fprintf('Bootstrapping Variance for %d ROIs using %d shuffles...\n', numROIs, numShuffles);
for s = 1:numShuffles
    % Extract the mean tuning curve for this specific shuffle
    shuffMeanTC = shuffMatrix(:, :, s);
    
    % Calculate the spatial structure of the noise (shuffled tuning var)
    shuffTuningVar = var(shuffMeanTC, 0, 2);
    
    % Ratio for this shuffle (Noise remains constant as it is a property of the ROI)
    shuffRatioDist(:, s) = posAllVar ./ (shuffTuningVar + eps);
end

% Calculate Z-Selectivity: (Mean Shuffle Ratio - Real Ratio) / Std Shuffle Ratio
% We subtract Real from Mean because a SMALLER ratio is BETTER.
mu_shuff = mean(shuffRatioDist, 2);
sigma_shuff = std(shuffRatioDist, 0, 2);
zSelectivity = (mu_shuff - ratioVarToTuningVar) ./ (sigma_shuff + eps);

% Identify Significant Boutons based on Z-Score
sigIdx = find(zSelectivity > zThreshold);
nSig = length(sigIdx);

%% 4. Prepare Diagnostic Plots (Odd/Even Heatmaps)
if nSig > 0
    oddIdx = 1:2:numLaps;
    evenIdx = 2:2:numLaps;
    
    % Extract significant ROI data
    meanOdd = squeeze(mean(lapActivity(sigIdx, oddIdx, :), 2, 'omitnan'));
    meanEven = squeeze(mean(lapActivity(sigIdx, evenIdx, :), 2, 'omitnan'));
    
    % Normalize [0 1.2] for visualization (95th percentile scaling)
    normOdd = nan(size(meanOdd)); normEven = nan(size(meanEven));
    for r = 1:nSig
        refMax = prctile(meanOdd(r,:), 95) + eps;
        normOdd(r,:) = meanOdd(r,:) ./ refMax;
        normEven(r,:) = meanEven(r,:) ./ refMax;
    end
    
    % Sort by Odd lap peaks
    [~, peakPos] = max(normOdd, [], 2);
    [~, sortIdx] = sort(peakPos);
    
    %% 5. Plotting
    fig1 = figure('Position', [50 50 1200 850], 'Color', 'w');
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Panel 1: Histogram of Z-Selectivity
    ax1 = nexttile;
    histogram(ax1, zSelectivity, 50, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none');
    xline(ax1, zThreshold, 'r--', 'LineWidth', 2, 'Label', 'Sig Threshold');
    title('Population Spatial Structure (Z-Score)');
    xlabel('Z-Selectivity Index'); ylabel('ROI Count');
    
    % Panel 2: Signal vs Noise Scatter
    ax2 = nexttile;
    scatter(ax2, tuningVar, posAllVar, 20, zSelectivity, 'filled', 'MarkerFaceAlpha', 0.6);
    colormap(ax2, parula); cb = colorbar(ax2); ylabel(cb, 'Z-Selectivity');
    xlabel('Tuning Variance (Signal)'); ylabel('Mean Bin Var (Noise)');
    title('Selectivity Landscape');
    
    % Panel 3: Heatmap Odd Laps
    ax3 = nexttile;
    imagesc(ax3, normOdd(sortIdx, :));
    title(sprintf('Significant Odd Laps (n=%d)', nSig));
    
    % Panel 4: Heatmap Even Laps (Cross-Validation)
    ax4 = nexttile;
    imagesc(ax4, normEven(sortIdx, :));
    title('Significant Even Laps (Sorted by Odd)');
    
    % Apply colormap and landmark lines to heatmaps
    for ax = [ax3, ax4]
        colormap(ax, flipud(gray));
        set(ax, 'CLim', [0 1.2], 'TickDir', 'out', 'box', 'off', 'YDir', 'normal');
        hold(ax, 'on');
        for xpos = [40 80 120 160], xline(ax, xpos, 'r:', 'LineWidth', 1.2); end
    end
    
    title(t, sprintf('%s | %s | %s | Variance Shuffle Analysis', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, signalToUse), 'FontSize', 14);
    
    exportgraphics(fig1, filename, 'Resolution', 300);
    fprintf('Diagnostic plot saved to: %s\n', filename);
else
    warning('No boutons passed the Z > %0.1f threshold. No heatmaps generated.', zThreshold);
end

%% 6. Save Data to Struct
tuningCurveVariance.zSelectivity = zSelectivity;
tuningCurveVariance.ratioVarToTuningVar = ratioVarToTuningVar;
tuningCurveVariance.sigIdx = sigIdx;
tuningCurveVariance.zThresholdUsed = zThreshold;

% % Append to sessionROIData
% if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
%     save(sessionFileInfo.otherSessFilePaths.sessionROIData, "tuningCurveVariance", '-append');
%     fprintf('Saved updated variance metrics to sessionROIData.\n');
% end

end