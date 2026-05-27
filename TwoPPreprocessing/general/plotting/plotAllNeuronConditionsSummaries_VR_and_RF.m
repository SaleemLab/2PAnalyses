function plotAllNeuronConditionsSummaries_VR_and_RF(sessionFileInfo, response, applySmoothing, signalToUse)
% plotAllNeuronSummaries_VR_and_RF: Unified summary for VR and Receptive Fields.
% Corrected for BonVision: Row 1 = Bottom of screen. No data flip, just YDir normal.

if nargin < 3, applySmoothing = true; end
if nargin < 4, signalToUse = 'dFFNeuropilCorrected'; end 

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
pdfPath = fullfile(figSaveDir, sprintf('%s_ConditionsSummary_PlusRFs_%s.pdf', ...
    response.stimName , signalToUse));

% Load RF data if available
rfDataAvailable = false;
rfType = ''; 
if isfield(sessionFileInfo.otherSessFilePaths, 'sessionROIData') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file')
    vars = who('-file', sessionFileInfo.otherSessFilePaths.sessionROIData);
    if ismember('RFMapping', vars) 
        load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'RFMapping', 'RFMappingMetadata', 'roiInfo');
        rfDataAvailable = true;
        rfType = 'grid';
    elseif ismember('sparseNoiseRF', vars)
        % load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'sparseNoiseRF', 'roiInfo');
        rfDataAvailable = true;
        rfType = 'sparse';
    end
end

data = response.lapPositionActivity.(signalToUse);
[nROIs, nRows, nBins] = size(data);
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
    roiActivity = squeeze(data(neuronIdx, :, :));
    if all(isnan(roiActivity), 'all'), continue; end
    
    if applySmoothing
        w = gausswin(10); w = w / sum(w);
        for iL = 1:nRows
            trace = roiActivity(iL, :);
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace); trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace); smoothed(nanMask) = NaN;
            roiActivity(iL, :) = smoothed;
        end
    end
    
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [10 50 1800 450]);
    t = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Color code title by iscell
    if iscell(neuronIdx, 1) == 1
        txtColor = [0 0.6 0]; % Green
    else
        txtColor = [1 0 0];   % Red
    end
    
    % title(t, sprintf('%s : %s | ROI %d ', sessionFileInfo.animal_name, sessionFileInfo.session_name, neuronIdx,'FontWeight', 'bold', 'Color', txtColor);

    %% Column 1: Omit vs Baseline
    nexttile; hold on;
    renderConditionWithSEM(conds, {'baseline', 'omit'}, response, roiActivity, colorMap);
    title('Omits'); ylabel(signalToUse); xlabel('Position (cm)');
    set(gca, 'XTick', [1 40 80 120 160 200]); 
    
    %% Column 2: Swap vs Baseline
    nexttile; hold on; 
    renderConditionWithSEM(conds, {'baseline', 'swap'}, response, roiActivity, colorMap);
    title('Swaps'); xlabel('Position (cm)');
    set(gca, 'XTick', [1 40 80 120 160 200]); 
    
    %% Column 3: Heatmap
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
    
    %% Column 4: RF
    axRF = nexttile; hold on;
    if rfDataAvailable
        if strcmp(rfType, 'grid')
            renderGridRF(axRF, RFMapping(neuronIdx), RFMappingMetadata, applySmoothing);
        else
            rfRaw = flipud(sparseNoiseRF.initMap{neuronIdx}(:, :, end, 4));
            imagesc(linspace(-70, 20, size(rfRaw, 2)), linspace(-20, 40, size(rfRaw, 1)), imgaussfilt(rfRaw, 1));
            set(gca, 'YDir', 'normal'); 
            colormap(axRF, 'parula'); colorbar;
            xline(0, 'k:', 'Alpha', 0.5); yline(0, 'k:', 'Alpha', 0.5);
            xlabel('Azimuth (\circ)'); ylabel('Elevation (\circ)');
            title('Receptive field (SVD)');
        end
    else
        text(0.5, 0.5, 'No RF Data', 'HorizontalAlignment', 'center');
        axis off;
    end
    
    exportgraphics(fig, pdfPath, 'Append', true);
    close(fig);
end
end

%% Helper function: RF grid superimposed with line plots
function renderGridRF(ax, roiRF, metadata, applySmoothing)
    uAz = metadata.uAz; uEl = metadata.uEl;
    imagesc(ax, uAz, uEl, imgaussfilt(roiRF.meanGridResponse, 0.8));
    colormap(ax, flipud(gray));
    set(ax, 'YDir', 'normal');
    cb = colorbar(ax); cb.Label.String = '\Delta F/F';
    
    % Use a Deep Orange for RF traces to distinguish from Omit Crimson
    rfLineColor = [1.0 0.4 0.0]; 
    maxAmp = roiRF.peakAmplitude;
    dEl = abs(uEl(1)-uEl(2)); dAz = abs(uAz(2)-uAz(1));
    vS = dEl * 0.4; hS = dAz * 0.85; 
    tVec = metadata.timeVector;
    
    for r = 1:size(roiRF.meanGridResponse, 1)
        for c = 1:size(roiRF.meanGridResponse, 2)
            tr = roiRF.meanTemporalResponse(:,r,c);
            trB = roiRF.meanBlankResponse;
            
            if applySmoothing
                tr = smoothdata(tr, 'gaussian', 3);
                trB = smoothdata(trB, 'gaussian', 3);
            end
            
            tNorm = (tVec - min(tVec)) / (max(tVec) - min(tVec));
            sX = (tNorm - 0.5) * hS + uAz(c);
            sY_B = (trB / maxAmp * vS) + uEl(r);
            sY_S = (tr / maxAmp * vS) + uEl(r);
            
            plot(ax, sX, sY_B, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'LineStyle', ':');
            plot(ax, sX, sY_S, 'Color', rfLineColor, 'LineWidth', 1.2);
        end
    end
    title(ax, 'Receptive field'); 
    xlabel(ax, 'Azimuth (\circ)'); ylabel('Elevation (\circ)');
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