function plotSingleROI_BaselineStack(sessionFileInfo, response, neuronIdx, applySmoothing, signalToUse)
% plotSingleROI_BaselineStack: Normalized stability check with Even scaled to Odd.

if nargin < 4, applySmoothing = false; end
if nargin < 5, signalToUse = 'dFFNeuropilCorrected'; end

%% 1. Data Prep
data = response.lapPositionActivity.(signalToUse);
roiActivity = squeeze(data(neuronIdx, :, :));
conds = fieldnames(response.trialIndicesByCondition);
baseIdx = find(contains(lower(conds), 'baseline'), 1);
allBaseIDs = response.trialIndicesByCondition.(conds{baseIdx});

oddIDs = allBaseIDs(1:2:end);
evenIDs = allBaseIDs(2:2:end);

if applySmoothing
    w = gausswin(5); w = w / sum(w); 
    for iL = 1:size(roiActivity, 1)
        trace = roiActivity(iL, :);
        if all(isnan(trace)), continue; end
        nanMask = isnan(trace); trace(nanMask) = 0;
        smoothed = filtfilt(w, 1, trace); smoothed(nanMask) = NaN;
        roiActivity(iL, :) = smoothed;
    end
end

muOddRaw = mean(roiActivity(oddIDs, :), 1, 'omitnan');
minVal = min(muOddRaw);
maxVal = max(muOddRaw);
rangeVal = maxVal - minVal;
if rangeVal > 0
    roiActivity = (roiActivity - minVal) / rangeVal;
else
    roiActivity = roiActivity - minVal; 
end

muOdd = mean(roiActivity(oddIDs, :), 1, 'omitnan');
[~, peakBin] = max(muOdd);
[~, sortIdxOdd] = sort(roiActivity(oddIDs, peakBin), 'descend');
[~, sortIdxEven] = sort(roiActivity(evenIDs, peakBin), 'descend');

splitData(1).name = 'Baseline Odd';
splitData(1).ids = oddIDs(sortIdxOdd);
splitData(1).color = [0.2 0.2 0.2];
splitData(2).name = 'Baseline Even';
splitData(2).ids = evenIDs(sortIdxEven);
splitData(2).color = [0.85 0.33 0.1];

nBins = size(roiActivity, 2);
xPos = 1:nBins;
tickLocs = [40, 80, 120, 160, 200]; % Added 200 for clarity
tickLabels = {'40', '80', '120', '160', '200'};

%% 2. Limits and Figure
muEven = mean(roiActivity(evenIDs, :), 1, 'omitnan');
yLimit = max([1, max(muEven)]) * 1.2;
fig = figure('Color', 'w', 'Position', [100 50 1100 750]);
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % Changed padding to compact
title(t, sprintf('%s | ROI %d ', ...
    sessionFileInfo.session_name, neuronIdx), 'FontSize', 14, 'FontWeight', 'bold');

%% 3. Plotting Loop
for iR = 1:2
    currIDs = splitData(iR).ids;
    currColor = splitData(iR).color;
    
    % --- Raster ---
    nexttile; hold on;
    imagesc(xPos, 1:length(currIDs), roiActivity(currIDs, :));
    colormap(gca, flipud(gray));
    set(gca, 'YDir', 'reverse', 'Box', 'off', 'XTick', tickLocs);
    
    % FIX: Explicitly set X-limits to remove empty side space
    xlim([1 nBins]); 
    axis tight; % Ensure data fills the box
    
    clim([0 1]); 
    colorbar;
    ylabel('Trial #', 'FontSize', 12);
    
    text(-0.2, 0.5, splitData(iR).name, 'Units', 'normalized', 'Rotation', 90, ...
        'FontWeight', 'bold', 'Color', currColor, 'HorizontalAlignment', 'center');
    
    if iR == 2, set(gca, 'XTickLabel', tickLabels); xlabel('Position (cm)', 'FontSize', 12);
    else, set(gca, 'XTickLabel', []); end

    % --- Tuning Curve ---
    nexttile; hold on;
    mu = mean(roiActivity(currIDs, :), 1, 'omitnan');
    sem = std(roiActivity(currIDs, :), 0, 1, 'omitnan') ./ sqrt(length(currIDs));
    
    fill([xPos, fliplr(xPos)], [mu+sem, fliplr(mu-sem)], currColor, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(xPos, mu, 'Color', currColor, 'LineWidth', 3);
    
    for p = tickLocs, xline(p, 'k:', 'Alpha', 0.2); end
    
    % FIX: Explicitly set X-limits here too
    xlim([1 nBins]); 
    set(gca, 'Box', 'off', 'XTick', tickLocs, 'ylim', [-0.2 yLimit]);
    ylabel('Norm. Activity', 'FontSize', 12);
    
    if iR == 2, set(gca, 'XTickLabel', tickLabels); xlabel('Position (cm)', 'FontSize', 12);
    else, set(gca, 'XTickLabel', []); end
end

% Save
saveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'ROI_Summaries');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
saveas(fig, fullfile(saveDir, sprintf('ROI_%d_Normalised.png', neuronIdx)));

end