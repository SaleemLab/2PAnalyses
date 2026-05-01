function plotSingleROI_ConditionStack(sessionFileInfo, response, neuronIdx, signalToUse)
% plotSingleROI_ConditionStack: Robust vertical stack with dynamic x-labeling.

if nargin < 4, signalToUse = 'dFFNeuropilCorrected'; end

%% 1. Data Prep
data = response.lapPositionActivity.(signalToUse);
roiActivity = squeeze(data(neuronIdx, :, :));
conds = fieldnames(response.trialIndicesByCondition);

% Filter for conditions that actually have trials to avoid empty rows
validConds = {};
for i = 1:length(conds)
    if ~isempty(response.trialIndicesByCondition.(conds{i}))
        validConds{end+1} = conds{i};
    end
end
numValid = length(validConds);

nBins = size(roiActivity, 2);
xPos = 1:nBins;
tickLocs = [40, 80, 120, 160]; 
tickLabels = {'40', '80', '120', '160'};
colors = {[0 0 0], [0 0.45 0.74], [0.85 0.33 0.1], [0.2 0.6 0.2], [0.5 0.2 0.5], [0.3 0.7 0.9]};

%% 2. Optimized Y-Axis Limit (Based on Means)
maxMean = 0;
for iC = 1:numValid
    mu = mean(roiActivity(response.trialIndicesByCondition.(validConds{iC}), :), 1, 'omitnan');
    maxMean = max(maxMean, max(mu));
end
yLimit = maxMean * 1; 
if yLimit <= 0 || isnan(yLimit), yLimit = 1; end

%% 3. Figure Setup
fig = figure('Color', 'w', 'Position', [100 50 1000 950]);
t = tiledlayout(numValid, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
title(t, sprintf('%s | ROI %d | %s', sessionFileInfo.session_name, neuronIdx, signalToUse), ...
    'FontSize', 14, 'FontWeight', 'bold');

% Reference Baseline
baseIdx = find(contains(lower(validConds), 'baseline'), 1);
if isempty(baseIdx), baseIdx = 1; end
baseMu = mean(roiActivity(response.trialIndicesByCondition.(validConds{baseIdx}), :), 1, 'omitnan');

%% 4. Vertical Stack Loop
for iC = 1:numValid
    currCond = validConds{iC};
    trialIDs = response.trialIndicesByCondition.(currCond);
    currColor = colors{mod(iC-1, length(colors))+1};
    
    % --- LEFT COLUMN: RASTER ---
    axL = nexttile; hold on;
    imagesc(xPos, 1:length(trialIDs), roiActivity(trialIDs, :));
    colormap(axL, flipud(gray));
    set(gca, 'YDir', 'reverse', 'Box', 'off', 'XTick', tickLocs);
    
    ylabel(sprintf('%s\n(Laps)', strrep(currCond, '_', ' ')), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', currColor);
    
    % Only add X-labels and Tick labels to the very last row
    if iC == numValid
        set(gca, 'XTickLabel', tickLabels);
        xlabel('Position (cm)', 'FontWeight', 'bold');
    else
        set(gca, 'XTickLabel', []);
    end

    % --- RIGHT COLUMN: TUNING CURVE ---
    axR = nexttile; hold on;
    plot(xPos, baseMu, 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'LineStyle', '--');
    
    mu = mean(roiActivity(trialIDs, :), 1, 'omitnan');
    sem = std(roiActivity(trialIDs, :), 0, 1, 'omitnan') ./ sqrt(length(trialIDs));
    fill([xPos, fliplr(xPos)], [mu+sem, fliplr(mu - sem)], currColor, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(xPos, mu, 'Color', currColor, 'LineWidth', 3);
    
    for p = tickLocs, xline(p, 'k:', 'Alpha', 0.2); end
    set(gca, 'Box', 'off', 'XTick', tickLocs, 'ylim', [-1 yLimit]);
    
    % Only add X-labels and Tick labels to the very last row
    if iC == numValid
        set(gca, 'XTickLabel', tickLabels);
        xlabel('Position (cm)', 'FontWeight', 'bold');
    else
        set(gca, 'XTickLabel', []);
    end
    
    if iC == 1, title('Mean \Delta F/F \pm SEM', 'FontSize', 10); end
end

%% 5. Save
saveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'ROI_Summaries');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
saveas(fig, fullfile(saveDir, sprintf('ROI_%d_Summary_Stacked.png', neuronIdx)));

end