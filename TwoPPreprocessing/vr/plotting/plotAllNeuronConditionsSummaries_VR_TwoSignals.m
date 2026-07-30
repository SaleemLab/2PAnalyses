function plotAllNeuronConditionsSummaries_VR_TwoSignals(sessionFileInfo, response, applySmoothing, signalToUse1, signalToUse2)
% plotAllNeuronConditionsSummaries_VR_TwoSignals: One page per ROI, two rows.
% Row 1 = signalToUse1 (Omits | Swaps | Heatmap)
% Row 2 = signalToUse2 (Omits | Swaps | Heatmap)
% RF panel removed.

if nargin < 3, applySmoothing = true; end
if nargin < 4, signalToUse1 = 'dFFNeuropilCorrected'; end
if nargin < 5, signalToUse2 = 'dFFNeuropilCorrected'; end

%% Load isCell and find indices
vrStimNames = {'LandManipCorridor', 'BaselineCorridor', 'VRCorr'};

if contains(response.stimName, 'CombinedRuns')
    allStimNames = {sessionFileInfo.stimFiles.name};
    isVRMatch = cellfun(@(x) any(contains(x, vrStimNames)), allStimNames);
    isNotCombined = cellfun(@(x) ~contains(x, 'CombinedRuns'), allStimNames);
    validIndices = find(isVRMatch & isNotCombined);
    if ~isempty(validIndices)
        stimIdx = validIndices(2);
    else
        stimIdx = find(isNotCombined, 1);
    end
else
    stimIdx = find(strcmp(response.stimName, {sessionFileInfo.stimFiles.name}), 1);
end

load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'iscell');

%% Create figure directory
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
pdfPath = fullfile(figSaveDir, sprintf('%s_ConditionsSummary_TwoSignals_%s_vs_%s.pdf', ...
    response.stimName, signalToUse1, signalToUse2));

data1 = response.lapPositionActivity.(signalToUse1);
data2 = response.lapPositionActivity.(signalToUse2);
[nROIs, nRows, nBins] = size(data1);
conds = fieldnames(response.trialIndicesByCondition);

%% Color Definitions - Updated for better differentiation
% Omit conditions: Crimson and Vivid Orange
warmColors = [0.85 0.08 0.23;
              1.00 0.50 0.00];

% Swap conditions: Deep Blue and Bright Cyan
coolColors = [0.00 0.45 0.74;
              0.00 0.80 0.80];

colorMap = struct();
omitCount = 1; swapCount = 1;
for iC = 1:length(conds)
    name = lower(conds{iC});
    if contains(name, 'baseline') || contains(name, 'norm'), colorMap.(conds{iC}) = [0 0 0];
    elseif contains(name, 'omit')
        colorMap.(conds{iC}) = warmColors(mod(omitCount-1, size(warmColors,1))+1, :);
        omitCount = omitCount + 1;
    elseif contains(name, 'swap')
        colorMap.(conds{iC}) = coolColors(mod(swapCount-1, size(coolColors,1))+1, :);
        swapCount = swapCount + 1;
    else, colorMap.(conds{iC}) = [0.5 0.5 0.5];
    end
end

%% Plot per ROI
for neuronIdx = 1:nROIs
    roiActivity1 = squeeze(data1(neuronIdx, :, :));
    roiActivity2 = squeeze(data2(neuronIdx, :, :));
    if all(isnan(roiActivity1), 'all') && all(isnan(roiActivity2), 'all'), continue; end

    if applySmoothing
        w = gausswin(10); w = w / sum(w);
        for iL = 1:nRows
            trace1 = roiActivity1(iL, :);
            if ~all(isnan(trace1))
                nanMask1 = isnan(trace1); trace1(nanMask1) = 0;
                smoothed1 = filtfilt(w, 1, trace1); smoothed1(nanMask1) = NaN;
                roiActivity1(iL, :) = smoothed1;
            end

            trace2 = roiActivity2(iL, :);
            if ~all(isnan(trace2))
                nanMask2 = isnan(trace2); trace2(nanMask2) = 0;
                smoothed2 = filtfilt(w, 1, trace2); smoothed2(nanMask2) = NaN;
                roiActivity2(iL, :) = smoothed2;
            end
        end
    end

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [10 50 1400 900]);
    t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Color code title by iscell
    if iscell(neuronIdx, 1) == 1
        txtColor = [0 0.6 0]; % Green
    else
        txtColor = [1 0 0];   % Red
    end
    title(t, sprintf('%s : %s | ROI %d', sessionFileInfo.animal_name, sessionFileInfo.session_name, neuronIdx), ...
        'FontWeight', 'bold', 'Color', txtColor);

    %% Row 1: signalToUse1
    plotSignalRow(roiActivity1, conds, response, colorMap, nBins, signalToUse1);

    %% Row 2: signalToUse2
    plotSignalRow(roiActivity2, conds, response, colorMap, nBins, signalToUse2);

    exportgraphics(fig, pdfPath, 'Append', true);
    close(fig);
end
end

%% Helper function: plot one row of 3 panels (Omits | Swaps | Heatmap) for a given signal
function plotSignalRow(roiActivity, conds, response, colorMap, nBins, signalName)
    nRows = size(roiActivity, 1);

    %% Panel 1: Omit vs Baseline
    nexttile; hold on;
    renderConditionWithSEM(conds, {'baseline', 'omit'}, response, roiActivity, colorMap);
    title('Omits'); ylabel(signalName); xlabel('Position (cm)');
    set(gca, 'XTick', [1 40 80 120 160 200]);

    %% Panel 2: Swap vs Baseline
    nexttile; hold on;
    renderConditionWithSEM(conds, {'baseline', 'swap'}, response, roiActivity, colorMap);
    title('Swaps'); xlabel('Position (cm)');
    set(gca, 'XTick', [1 40 80 120 160 200]);

    %% Panel 3: Heatmap
    axH = nexttile; hold on;
    normAct = normalize(roiActivity, 2, 'range');
    imagesc(1:nBins, 1:nRows, normAct);
    colormap(axH, flipud(gray));

    gutterX = -15;
    for iC = 1:length(conds)
        IDs = response.trialIndicesByCondition.(conds{iC});
        if ~isempty(IDs)
            s = scatter(repmat(gutterX, size(IDs)), IDs, 25, colorMap.(conds{iC}), 'filled');
            s.Clipping = 'off';
        end
    end
    xlim([gutterX-5 nBins]); set(gca, 'YDir', 'normal');
    xlabel('Position (cm)'); ylabel('Lap #'); title('Lap Activity');
    for p = [40 80 120 160], xline(p, 'k--'); end
    cbH = colorbar; cbH.Label.String = 'Norm. Activity';
end

%% Helper function: render with SEM
function renderConditionWithSEM(conds, keywords, response, roiActivity, colorMap)
    lgdLines = []; lgdNames = {};
    x = 1:size(roiActivity, 2);
    for iC = 1:length(conds)
        name = lower(conds{iC});
        if any(cellfun(@(k) contains(name, k), keywords))
            rowIdx = response.trialIndicesByCondition.(conds{iC});
            if ~isempty(rowIdx)
                rowIdx = rowIdx(rowIdx <= size(roiActivity, 1));
                mu = mean(roiActivity(rowIdx, :), 1, 'omitnan');

                if length(rowIdx) > 1
                    sem = std(roiActivity(rowIdx, :), 0, 1, 'omitnan') ./ sqrt(length(rowIdx));
                    fill([x fliplr(x)], [mu + sem, fliplr(mu - sem)], colorMap.(conds{iC}), ...
                        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                end

                l = plot(x, mu, 'Color', colorMap.(conds{iC}), 'LineWidth', 2.5);
                lgdLines(end+1) = l;
                lgdNames{end+1} = sprintf('%s (n=%d)', strrep(conds{iC},'_',' '), length(rowIdx));
            end
        end
    end
    for p = [40 80 120 160], xline(p, 'k--'); end
    if ~isempty(lgdLines), legend(lgdLines, lgdNames, 'Location', 'northeast', 'FontSize', 7); end
end