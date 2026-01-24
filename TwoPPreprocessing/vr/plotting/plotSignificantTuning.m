function plotSignificantTuning(response, nullDist_PeakTuningMetric, signalToUse)
% Generates two figures: 1) All ROIs and 2) ROIs significant > 95th percentile.
% Uses the odd/even lap split to demonstrate spatial stability.

%% 1. Prepare Data
lapActivity = response.lapPositionActivity.(signalToUse);
numROIs = size(lapActivity, 1);

% Split odd and even laps
oddLaps = lapActivity(:, 1:2:end, :); 
evenLaps = lapActivity(:, 2:2:end, :);

% Average and Normalize [0 to 1]
meanOdd = squeeze(nanmean(oddLaps, 2));
meanEven = squeeze(nanmean(evenLaps, 2));
normOdd = normalize(meanOdd, 2, 'range');
normEven = normalize(meanEven, 2, 'range');

% Get Significance Index
sigIdx = nullDist_PeakTuningMetric.isSignificantByPeakShuffling.(signalToUse);
significantROIs = find(sigIdx == 1);

%% 2. Generate Figure 1: All ROIs (Unfiltered)
[~, peakIdxAll] = max(normOdd, [], 2);
[~, sortIdxAll] = sort(peakIdxAll);

fig1 = figure('Name', 'Spatial Tuning: All ROIs', 'Position', [100 100 1000 500]);
t1 = tiledlayout(1, 2, 'Padding', 'compact');
title(t1, 'All Recorded ROIs (Unfiltered)');

% Plot Odd
ax1 = nexttile; imagesc(normOdd(sortIdxAll,:));
title('Odd Laps'); ylabel('ROI #'); colormap(flipud(gray));
applyFormatting(ax1);

% Plot Even
ax2 = nexttile; imagesc(normEven(sortIdxAll,:));
title('Even Laps'); colormap(flipud(gray));
applyFormatting(ax2);

%% 3. Generate Figure 2: Significant ROIs Only (> 95th Percentile)
if ~isempty(significantROIs)
    normOddSig = normOdd(significantROIs, :);
    normEvenSig = normEven(significantROIs, :);
    
    [~, peakIdxSig] = max(normOddSig, [], 2);
    [~, sortIdxSig] = sort(peakIdxSig);
    
    fig2 = figure('Name', 'Significant Spatial Tuning', 'Position', [150 150 1000 500]);
    t2 = tiledlayout(1, 2, 'Padding', 'compact');
    title(t2, sprintf('Significant ROIs (>95th Percentile Shuffle, n=%d/%d)', ...
        length(significantROIs), numROIs));
    
    % Plot Odd Significant
    ax3 = nexttile; imagesc(normOddSig(sortIdxSig,:));
    title('Odd Laps (Significant)'); ylabel('Significant ROI #'); colormap(flipud(gray));
    applyFormatting(ax3);
    
    % Plot Even Significant
    ax4 = nexttile; imagesc(normEvenSig(sortIdxSig,:));
    title('Even Laps (Significant)'); colormap(flipud(gray));
    applyFormatting(ax4);
else
    warning('No significant ROIs found to plot.');
end

%% Internal Helper: Apply your specific formatting parameters
    function applyFormatting(ax)
        set(ax, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
        xline(ax, [50 70 90 110], 'k--', 'LineWidth', 1.2); % Markers for track features
        xticks(ax, [0 50 70 90 110 140]);
        xticklabels(ax, {'0', '50', '70', '90', '110', '140'});
        xlabel(ax, 'Position (cm)');
    end
end