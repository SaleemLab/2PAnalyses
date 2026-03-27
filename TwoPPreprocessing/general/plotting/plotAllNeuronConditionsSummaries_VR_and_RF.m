function plotAllNeuronConditionsSummaries_VR_and_RF(sessionFileInfo, response, applySmoothing, signalToUse)
% plotAllNeuronSummaries_VR_and_RF: Unified summary for VR and Receptive Fields.
% Corrected for BonVision: Row 1 = Bottom of screen. No data flip, just YDir normal.
% This version checks if sparsenoise (V1 somas) or RF mapping (RSP boutons) exists and runs one or the other. 


if nargin < 3, applySmoothing = true; end
if nargin < 4, signalToUse = 'dFFNeuropilCorrected'; end 

%% Load isCell from any of the stim files: 
vrStimNames = {'LandManipCorridor', 'BaselineCorridor', 'VRCorr'};
%allStimNames = {sessionFileInfo.stimFiles.name};

%Use combined runs as the combined pipeline does not merge the

if contains(response.stimName, 'CombinedRuns')

    % Find indices that match any of the vrStimNames
    allStimNames = {sessionFileInfo.stimFiles.name};
    isVRMatch = cellfun(@(x) any(contains(x, vrStimNames)), allStimNames);
    isNotCombined = cellfun(@(x) ~contains(x, 'CombinedRuns'), allStimNames);

    validIndices = find(isVRMatch & isNotCombined);

    if ~isempty(validIndices)
        stimIdx = validIndices(2); % Take the first one found
    else
        % Fallback: if no specific VR match, just take the first
        % non-combined file; This is temporary and needs to change. Only
        % M25132 has this issue on day 1 of experience 
        stimIdx = find(isNotCombined, 1);
    end
else
    % Standard match for individual runs
    stimIdx = find(strcmp(response.stimName, {sessionFileInfo.stimFiles.name}), 1);
end
% Load classification data (Column 1: 1=cell, 0=not a cell)
load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'iscell');

%% Create figure 
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
pdfPath = fullfile(figSaveDir, sprintf('%s_%s_ConditionsSummary_PlusRFs_%s.pdf', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, signalToUse));

% Load RF data if available (Check for SparseNoise vs Grid Mapping)
rfDataAvailable = false;
rfType = ''; 

if isfield(sessionFileInfo.otherSessFilePaths, 'sessionROIData') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file')
    vars = who('-file', sessionFileInfo.otherSessFilePaths.sessionROIData); % who lists variables in workspace 
    if ismember('RFMapping', vars) 
        % plane identity will say 0 for all boutons (z-motion corrected; this is not the selected plane) 
        load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'RFMapping', 'RFMappingMetadata', 'roiInfo');
        rfDataAvailable = true;
        rfType = 'grid';
    elseif ismember('sparseNoiseRF', vars)
        load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'sparseNoiseRF', 'roiInfo');
        rfDataAvailable = true;
        rfType = 'sparse';
    end
end

data = response.lapPositionActivity.(signalToUse);
[nROIs, nRows, nBins] = size(data);
conds = fieldnames(response.trialIndicesByCondition);


% Color Definitions
warmColors = [0.9 0.2 0.2; 1 0.6 0]; 
coolColors = [0 0.45 0.74; 0 0.8 0.8]; 
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
    
    % Smoothing
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
    
    % Color code title by iscell: Green (1) or Red (0)
    if iscell(neuronIdx, 1) == 1
        txtColor = [0 0.6 0]; % Greenish? 
    else
        txtColor = [1 0 0]; % Red
    end
    
    title(t, sprintf('%s : %s | ROI %d | Plane %d', sessionFileInfo.animal_name, sessionFileInfo.session_name, neuronIdx, roiInfo.roiPlaneIdentity(neuronIdx)), ...
        'FontWeight', 'bold', 'Color', txtColor);

    %% Omit vs Baseline
    nexttile; hold on;
    renderConditionWithSEM(conds, {'baseline', 'omit'}, response, roiActivity, colorMap);
    title('Omits'); ylabel(signalToUse); xlabel('Position (cm)');
    set(gca, 'XTick', [1 40 80 120 160 200]); 
    
    %% Swap vs Baseline
    nexttile; hold on; 
    renderConditionWithSEM(conds, {'baseline', 'swap'}, response, roiActivity, colorMap);
    title('Swaps'); xlabel('Position (cm)');
    set(gca, 'XTick', [1 40 80 120 160 200]); 
    
    %% Heatmap
    axH = nexttile; hold on;
    normAct = normalize(roiActivity, 2, 'range');
    imagesc(1:nBins, 1:nRows, normAct);
    colormap(axH, flipud(gray));
    
    gutterX = -15;
    for iC = 1:length(conds)
        IDs = response.trialIndicesByCondition.(conds{iC});
        % With the unified indexing, absIDs directly match the row indices
        if ~isempty(IDs)
            s = scatter(repmat(gutterX, size(IDs)), IDs, 25, colorMap.(conds{iC}), 'filled');
            s.Clipping = 'off'; 
        end
    end
    xlim([gutterX-5 nBins]); set(gca, 'YDir', 'normal');
    xlabel('Position (cm)'); ylabel('Lap #'); title('Lap Activity');
    for p = [40 80 120 160], xline(p, 'k--'); end
    cbH = colorbar; cbH.Label.String = 'Norm. Activity';
    
    %% RF
    axRF = nexttile; hold on;
    if rfDataAvailable
        if strcmp(rfType, 'grid')
            % Render Grid RF using the helper function
            renderGridRF(axRF, RFMapping(neuronIdx), RFMappingMetadata, applySmoothing);
        else
            % SNT (SVD)
            rfRaw = flipud(sparseNoiseRF.initMap{neuronIdx}(:, :, end, 4)); 
            imagesc(linspace(-70, 20, size(rfRaw, 2)), linspace(-20, 40, size(rfRaw, 1)), imgaussfilt(rfRaw, 1));
            set(gca, 'YDir', 'normal'); 
            colormap(axRF, 'parula'); colorbar;
            xline(0, 'k:', 'Alpha', 0.5); yline(0, 'k:', 'Alpha', 0.5);
            xlabel('Azimuth (\circ)'); 
            ylabel('Elevation (\circ)');
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

%% Helper function: RF grid superimposed with line plots (including blank trials) 
function renderGridRF(ax, roiRF, metadata, applySmoothing)
    uAz = metadata.uAz; uEl = metadata.uEl;
    % Plot the smoothed grid heatmap background
    imagesc(ax, uAz, uEl, imgaussfilt(roiRF.meanGridResponse, 0.8));
    colormap(ax, flipud(gray));
    set(ax, 'YDir', 'normal');
    cb = colorbar(ax); cb.Label.String = 'dFF';
    
    % this 
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
            
            % Overlay blank (gray) and solid (red) 
            plot(ax, sX, sY_B, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'LineStyle', ':');
            plot(ax, sX, sY_S, 'r', 'LineWidth', 1);
        end
    end
    title(ax, 'Receptive field'); 
    xlabel(ax, 'Azimuth (\circ)'); 
    ylabel(ax, 'Elevation (\circ)');
end

%%  get mean plus sem
function renderConditionWithSEM(conds, keywords, response, roiActivity, colorMap)
    % Removed 'compLaps' from arguments as it is no longer needed
    lgdLines = []; lgdNames = {};
    x = 1:size(roiActivity, 2);
    for iC = 1:length(conds)
        name = lower(conds{iC});
        if any(cellfun(@(k) contains(name, k), keywords))
            % absIDs are now the DIRECT row indices and match lap-position-activity idx (e.g., [1, 5, 10...])
            rowIdx = response.trialIndicesByCondition.(conds{iC});
            
            if ~isempty(rowIdx)
                % rowIdx doesn't exceed the number of rows in activity
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
    % 
    for p = [40 80 120 160], xline(p, 'k--'); end
    if ~isempty(lgdLines), legend(lgdLines, lgdNames, 'Location', 'northeast', 'FontSize', 7); end
end