function figHandle = plotSingleROI_ConditionStack(sessionFileInfo, response, neuronIdx, signalToUse)
% plotSingleROI_ConditionStack: Robust vertical stack with sequential lap counts.
% FIXED: Handles explicitly passed to all graphics functions to ensure 
% defaultAxesProperties and offsetAxes are accurately targeted and executed.

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
    nLapsInBlock = length(trialIDs);
    currColor = colors{mod(iC-1, length(colors))+1};
    
    % --- LEFT COLUMN: RASTER ---
    axL = nexttile(t); hold(axL, 'on'); % Explicit tile destination assignment
    
    imagesc(axL, xPos, 1:nLapsInBlock, roiActivity(trialIDs, :));
    colormap(axL, flipud(gray));
    
    % FIXED: Switched from gca to explicit axL handle targeting
    set(axL, 'YDir', 'reverse', 'Box', 'off', 'XTick', tickLocs, ...
             'YTick', [1, nLapsInBlock], 'YTickLabel', {'1', num2str(nLapsInBlock)});
    
    ylabel(axL, sprintf('%s\n(Laps)', strrep(currCond, '_', ' ')), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', currColor);
    
    if iC == numValid
        set(axL, 'XTickLabel', tickLabels);
        xlabel(axL, 'Position (cm)', 'FontWeight', 'bold');
    else
        set(axL, 'XTickLabel', []);
    end
    
    % Call properties on the handle directly
    defaultAxesProperties(axL, true);
    offsetAxes(axL)
    
    % --- RIGHT COLUMN: TUNING CURVE ---
    axR = nexttile(t); hold(axR, 'on'); % Explicit tile destination assignment
    plot(axR, xPos, baseMu, 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'LineStyle', '--');
    
    mu = mean(roiActivity(trialIDs, :), 1, 'omitnan');
    sem = std(roiActivity(trialIDs, :), 0, 1, 'omitnan') ./ sqrt(nLapsInBlock);
    fill(axR, [xPos, fliplr(xPos)], [mu+sem, fliplr(mu - sem)], currColor, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(axR, xPos, mu, 'Color', currColor, 'LineWidth', 3);
    
    for p = tickLocs
        xline(axR, p, 'k--', 'Alpha', 0.3); 
    end
    
    % FIXED: Switched from gca to explicit axR handle targeting
    set(axR, 'Box', 'off', 'XTick', tickLocs, 'ylim', [-1 yLimit]);
    
    if iC == numValid
        set(axR, 'XTickLabel', tickLabels);
        xlabel(axR, 'Position (cm)', 'FontWeight', 'bold');
    else
        set(axR, 'XTickLabel', []);
    end
    
    if iC == 1
        %title(axR, 'Mean \Delta F/F \pm SEM', 'FontSize', 10);
        title(axR, 'Mean Activity \pm SEM', 'FontSize', 10);
    end
    
    % Call properties on the handle directly
    defaultAxesProperties(axR, true);
    offsetAxes(axR)
end

%% 5. Save
saveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'ROI_Summaries');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
saveas(fig, fullfile(saveDir, sprintf('ROI_%d_Summary_Stacked.png', neuronIdx)));

if nargout > 0, figHandle = fig; end
end