function figHandle = plotSingleROI_3ConditionStack(sessionFileInfo, response, neuronIdx, signalToUse)
% plotSingleROI_3ConditionStack: Horizontal layout engine for vector editors.
% Row 1: Raster Trials (Top) | Row 2: Mean Tuning Curves (Bottom)
% FIXED: Restored individual numerical Y-tick labels ('1' and 'nLapsInBlock') 
% across ALL raster panels to accurately reflect varying block lengths.

if nargin < 4, signalToUse = 'dFFNeuropilCorrected'; end

%% 1. Data Prep
data = response.lapPositionActivity.(signalToUse);
roiActivity = squeeze(data(neuronIdx, :, :));
conds = {'Baseline', 'Swap_2_3', 'Omit_2', 'Omit_3'};
validConds = {};
for i = 1:length(conds)
    if isfield(response.trialIndicesByCondition, conds{i}) && ...
       ~isempty(response.trialIndicesByCondition.(conds{i}))
        validConds{end+1} = conds{i};
    end
end
numValid = length(validConds);
nBins = size(roiActivity, 2);
xPos = 1:nBins;
tickLocs = [40, 80, 120, 160]; 
tickLabels = {'40', '80', '120', '160'};
colors = { ...
    [0.000, 0.000, 0.000], ... % Baseline (Black)
    [0.000, 0.400, 1.000], ... % Swap_2_3 (Azul)
    [0.541, 0.012, 0.012], ... % Omit_2 (Blood)
    [0.824, 0.016, 0.176]  ... % Omit_3 (Cherry)
};

%% 2. Optimized Y-Axis Limit & Tick Rounding
maxMean = 0;
for iC = 1:numValid
    mu = mean(roiActivity(response.trialIndicesByCondition.(validConds{iC}), :), 1, 'omitnan');
    maxMean = max(maxMean, max(mu));
end
yLimit = maxMean * 1.05; 
if yLimit <= 0 || isnan(yLimit), yLimit = 1; end

% Clean rounding to the nearest whole integer for the upper tick mark
roundedUpperLimit = round(yLimit); 
yTicksCurve = [0, roundedUpperLimit];

% 
plotMaxLimit = roundedUpperLimit * 1.10;
if plotMaxLimit <= 0, plotMaxLimit = 1; end

%% Figure Canvas Setup (Landscape Format)
fig = figure('Color', 'w', 'Position', [100 50 850 360]);
baseIdx = find(contains(lower(validConds), 'baseline'), 1);
if isempty(baseIdx), baseIdx = 1; end
baseMu = mean(roiActivity(response.trialIndicesByCondition.(validConds{baseIdx}), :), 1, 'omitnan');

%% 
leftMargin   = 0.16; 
plotWidth    = (0.82 - leftMargin) / numValid - 0.03; 
plotHeight   = 0.33; 

row1_Bottom  = 0.52;  % Top Row: Rasters
row2_Bottom  = 0.13;  % Bottom Row: Means

xPositions = linspace(leftMargin, 0.95 - plotWidth, numValid);

%% 
for iC = 1:numValid
    currCond = validConds{iC};
    trialIDs = response.trialIndicesByCondition.(currCond);
    nLapsInBlock = length(trialIDs);
    currColor = colors{mod(iC-1, length(colors))+1};
    
    colX = xPositions(iC);
    
    if strcmpi(currCond, 'Baseline')
        cleanLabel = 'Base';
    else
        cleanLabel = strrep(currCond, '_', ' ');
    end
    

    axL = axes('Position', [colX, row1_Bottom, plotWidth, plotHeight]);
    hold(axL, 'on');
    
    imagesc(axL, xPos, 1:nLapsInBlock, roiActivity(trialIDs, :));
    colormap(axL, flipud(gray));
    set(axL, 'YDir', 'reverse', 'Box', 'off', 'XTick', tickLocs, 'XTickLabel', []);
    
    % Set numerical bounds for every layout iteration to account for variable lap counts
    set(axL, 'YTick', [1, nLapsInBlock]);
    
    % FIXED: Explicitly enabled text label numbers for every unique condition column
    set(axL, 'YTickLabel', {'1', num2str(nLapsInBlock)});
    
    
    title(axL, cleanLabel, 'FontSize', 8, 'Color', currColor, 'FontWeight', 'bold');
    
    if iC == 1
        ylh = ylabel(axL, '# Laps', 'FontSize', 8, 'Color', currColor);
        set(ylh, 'Units', 'normalized', 'Position', [-0.40, 0.5, 0], ...
                 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    end
    
    set(axL, 'FontSize', 8, 'LineWidth', 0.8);
    
    defaultAxesProperties(axL, false); 
    offsetAxes(axL); 
    
    % Strip default background bounding boxes safely
    axL.XAxis.Color = 'none';
    axL.YAxis.Color = 'none';
    
    
    %  TUNING CURVE PANELS (BOTTOM ROW)

    axR = axes('Position', [colX, row2_Bottom, plotWidth, plotHeight]);
    hold(axR, 'on');
    
    plot(axR, xPos, baseMu, 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'LineStyle', '--');
    
    mu = mean(roiActivity(trialIDs, :), 1, 'omitnan');
    sem = std(roiActivity(trialIDs, :), 0, 1, 'omitnan') ./ sqrt(nLapsInBlock);
    fill(axR, [xPos, fliplr(xPos)], [mu+sem, fliplr(mu - sem)], currColor, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(axR, xPos, mu, 'Color', currColor, 'LineWidth', 1.2);
    
    for p = tickLocs
        xline(axR, p, 'k:', 'Alpha', 0.15); 
    end
    
    set(axR, 'Box', 'off', 'XTick', tickLocs, 'XTickLabel', tickLabels, ...
             'YTick', yTicksCurve, 'ylim', [-0.05 plotMaxLimit]);
    xlabel(axR, 'Position (cm)', 'FontSize', 8);
    
    if iC == 1
        set(axR, 'YTickLabel', {'0', num2str(roundedUpperLimit)});
        
        ylh_curve = ylabel(axR, 'Mean Activity', 'FontSize', 8);
        set(ylh_curve, 'Units', 'normalized', 'Position', [-0.40, 0.5, 0], ...
                       'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
                   
%         title(axR, 'Mean Activity \pm SEM', 'FontSize', 8, 'FontAngle', 'italic');
    else
        set(axR, 'YTickLabel', {'', ''});
    end
    
    set(axR, 'FontSize', 8, 'LineWidth', 0.8);
    defaultAxesProperties(axR, true); 
    offsetAxes(axR); 
end

%% 
headerText = sprintf('%s | ROI %d | %s', sessionFileInfo.session_name, neuronIdx, signalToUse);
annotation(fig, 'textbox', [0.1, 0.92, 0.8, 0.06], 'String', headerText, ...
    'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

%%  Save
saveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'ROI_Summaries_FORAS');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
exportgraphics(fig, fullfile(saveDir, sprintf('ROI_%d_Summary_Horizontal.pdf', neuronIdx)), 'ContentType', 'vector');

if nargout > 0, figHandle = fig; end
end