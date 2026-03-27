function plotPopulationActivityAcrossConditions_HalvesStableROIsOnly(sessionFileInfo, response, signalToUse, applySmoothing)
% plotPopulationActivityAcrossConditions - Only plots ROIs identified as stable
if nargin < 3; signalToUse = 'dFF'; end
if nargin < 4; applySmoothing = true; end

%% Load Stability Indices
% Load the stability data saved by checkHalvesCorrelation
try
    vars = load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'lapCorr_Halves');
    stableIdx = vars.lapCorr_Halves.stableIdx;
    fprintf('Found %d stable ROIs. Filtering data...\n', length(stableIdx));
catch
    warning('Could not find lapCorr_Halves in sessionROIData. Plotting all ROIs instead.');
    stableIdx = 1:size(response.lapPositionActivity.(signalToUse), 1);
end

%% output
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
filename = fullfile(figSaveDir, sprintf('%s_%s_%s_%s_HalvesStable_Pop_Conditions.png', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, signalToUse, response.stimName));

%% activity & Filtering
% Subset the activity matrix to only include stable ROIs immediately
lapActivity = response.lapPositionActivity.(signalToUse)(stableIdx, :, :);
[nROIs, ~, nBins] = size(lapActivity);
conds = fieldnames(response.trialIndicesByCondition);

if applySmoothing
    w = gausswin(10); w = w / sum(w);
    for iCell = 1:nROIs
        for iLap = 1:size(lapActivity, 2)
            trace = squeeze(lapActivity(iCell, iLap, :));
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace); trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace); smoothed(nanMask) = NaN;
            lapActivity(iCell, iLap, :) = smoothed;
        end
    end
end

%% Conditions and colours
warmColors = [0.9 0.2 0.2; 1 0.6 0]; 
coolColors = [0 0.45 0.74; 0 0.8 0.8]; 
dataMatrix = {}; labels = {}; nLapsPerCond = []; titleColors = {}; condTypes = {}; rawNames = {};

baseIdx = find(contains(lower(conds), 'baseline') | contains(lower(conds), 'norm'), 1);
if isempty(baseIdx), baseIdx = 1; end

% Add baseline odd/Even
baseLaps = response.trialIndicesByCondition.(conds{baseIdx});
for split = 1:2
    laps = baseLaps(split:2:end);
    dataMatrix{end+1} = squeeze(mean(lapActivity(:, laps, :), 2, 'omitnan'));
    labels{end+1} = iff(split == 1, 'Baseline Odd', 'Baseline Even');
    nLapsPerCond(end+1) = length(laps);
    titleColors{end+1} = [0 0 0]; 
    condTypes{end+1} = 'baseline';
    rawNames{end+1} = conds{baseIdx};
end

% Add others
omitCount = 1; swapCount = 1;
otherNames = conds(setdiff(1:length(conds), baseIdx));
for i = 1:length(otherNames)
    cName = otherNames{i};
    laps = response.trialIndicesByCondition.(cName);
    if isempty(laps), continue; end
    
    dataMatrix{end+1} = squeeze(mean(lapActivity(:, laps, :), 2, 'omitnan'));
    labels{end+1} = strrep(cName, '_', ' ');
    nLapsPerCond(end+1) = length(laps);
    rawNames{end+1} = cName;
    
    nameLow = lower(cName);
    if contains(nameLow, 'omit')
        titleColors{end+1} = warmColors(mod(omitCount-1, size(warmColors,1))+1, :);
        condTypes{end+1} = 'omit'; omitCount = omitCount + 1;
    elseif contains(nameLow, 'swap')
        titleColors{end+1} = coolColors(mod(swapCount-1, size(coolColors,1))+1, :);
        condTypes{end+1} = 'swap'; swapCount = swapCount + 1;
    else
        titleColors{end+1} = [0.5 0.5 0.5]; condTypes{end+1} = 'other';
    end
end

%% sort using baseline odd trials
normRef = normalize(dataMatrix{1}, 2, 'range');
[~, peakIdx] = max(normRef, [], 2);
[~, sortIdx] = sort(peakIdx);

%% plotting
nPlots = length(dataMatrix);
fig = figure('Color', 'w', 'Position', [50 100 350*nPlots 450], 'Visible', 'on');
t = tiledlayout(1, nPlots, 'TileSpacing', 'none', 'Padding', 'compact'); 
thickLine = 2.0; 
thinLine = 0.8;

for i = 1:nPlots
    ax = nexttile;
    normMean = normalize(dataMatrix{i}, 2, 'range');
    imagesc(normMean(sortIdx, :));
    colormap(ax, flipud(gray)); caxis([0 1]); axis tight;
    set(ax, 'Box', 'on', 'LineWidth', 0.6, 'XColor', 'k', 'YColor', 'k');
    
    tH = title(sprintf('%s\n(n=%d)', labels{i}, nLapsPerCond(i)), 'FontSize', 10);
    tH.Color = titleColors{i}; 
    set(gca, 'YDir', 'normal', 'TickDir', 'out');
    
    % Landmark logic
    landmarks = [40 80 120 160];
    currentName = lower(rawNames{i});
    
    for lIdx = 1:4
        p = landmarks(lIdx);
        isTarget = false;
        if ~strcmp(condTypes{i}, 'baseline') && contains(currentName, num2str(lIdx))
            isTarget = true;
        end
        if strcmp(condTypes{i}, 'omit') && ~any(regexp(currentName, '[1-4]')) && p == 120
            isTarget = true;
        end
        if isTarget
            xline(p, '--', 'Color', titleColors{i}, 'LineWidth', thickLine); 
        else
            xline(p, '--', 'Color', [0 0 0], 'LineWidth', thinLine); 
        end
    end
    
    if i == 1
        ylabel('Stable ROIs (sorted by baseline Odd laps)'); 
        yticks([1 nROIs]);
    else
        set(gca, 'YTickLabel', []);
    end
    xticks([1 40 80 120 160 200]); xlabel('Position (cm)');
    if i == nPlots, cb = colorbar; cb.Label.String = 'Norm. Activity'; end
end

header = sprintf('%s | %s | %s (Halves-Stable Only [r>0.4])', sessionFileInfo.animal_name, sessionFileInfo.session_name, signalToUse);
title(t, header, 'FontSize', 12, 'FontWeight', 'bold');

%% Export
exportgraphics(fig, filename, 'Resolution', 300);

    function out = iff(condition, trueVal, falseVal)
        if condition, out = trueVal; else, out = falseVal; end
    end
end