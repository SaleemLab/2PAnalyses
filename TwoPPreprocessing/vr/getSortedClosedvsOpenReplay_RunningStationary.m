% function [responseCL, responseOL] = getSortedClosedvsOpenReplay_RunningStationary(sessionFileInfo, CLStimName, OLStimName, signalToUse, doPlot)
% if nargin < 4; signalToUse = 'dFFNeuropilCorrected'; end
% if nargin < 5; doPlot = true; end
% 
% runThresh = 1; 
% stableThresh = 0.6;
% 
% %% 
% allPaths = {sessionFileInfo.stimFiles.Response}; 
% idxCL = find(contains(allPaths, CLStimName));
% idxOL = find(contains(allPaths, OLStimName));
% 
% if isempty(idxCL); error('CL path not found.'); end
% if isempty(idxOL); error('OL path not found.'); end
% 
% %% 
% L = load(allPaths{idxCL(1)}, 'lapPositionActivity', 'lapPositionRunningSpeed');
% responseCL.lapPositionActivity.(signalToUse) = L.lapPositionActivity.(signalToUse);
% responseCL.lapPositionRunningSpeed = L.lapPositionRunningSpeed;
% 
% L = load(allPaths{idxOL(1)}, 'lapPositionActivity', 'lapPositionRunningSpeed');
% responseOL.lapPositionActivity.(signalToUse) = L.lapPositionActivity.(signalToUse);
% responseOL.lapPositionRunningSpeed = L.lapPositionRunningSpeed;
% clear L; 
% 
% %% filter for stable rOIs
% stabilityMetricUsed = 'lapCorr_Halves';
% try
% 
%     vars = load(sessionFileInfo.otherSessFilePaths.sessionROIData,stabilityMetricUsed);
%     stableIdx = find(vars.lapCorr_Halves.rho >= stableThresh);
% catch
%     stableIdx = 1:size(responseCL.lapPositionActivity.(signalToUse), 1);
% end
% 
% actCL_stable = responseCL.lapPositionActivity.(signalToUse)(stableIdx, :, :);
% actOL_stable = responseOL.lapPositionActivity.(signalToUse)(stableIdx, :, :);
% 
% %% masking running and stationary
% maskRunCL = responseCL.lapPositionRunningSpeed > runThresh;
% maskRunOL = responseOL.lapPositionRunningSpeed > runThresh;
% 
% responseCL.stateProfiles.Running = getMaskedMean(actCL_stable, maskRunCL);
% responseCL.stateProfiles.Stationary = getMaskedMean(actCL_stable, responseCL.lapPositionRunningSpeed <= runThresh);
% responseOL.stateProfiles.Running = getMaskedMean(actOL_stable, maskRunOL);
% responseOL.stateProfiles.Stationary = getMaskedMean(actOL_stable, responseOL.lapPositionRunningSpeed <= runThresh);
% 
% if doPlot
%     nLapsCL = size(actCL_stable, 2);
%     oddIdx = 1:2:nLapsCL;
%     evenIdx = 2:2:nLapsCL;
%     w = gausswin(15); w = w / sum(w);
% 
%     % speed for plotting 
%     maxSpd = max([max(responseCL.lapPositionRunningSpeed(:)), max(responseOL.lapPositionRunningSpeed(:))]);
%     speedLimit = ceil(maxSpd / 5) * 5; 
% 
%     % Smooth Raw Means for all conditions
%     sOdd  = smoothData(squeeze(mean(actCL_stable(:, oddIdx, :), 2, 'omitnan')), w);
%     sEven = smoothData(squeeze(mean(actCL_stable(:, evenIdx, :), 2, 'omitnan')), w);
%     sRunOL  = smoothData(responseOL.stateProfiles.Running, w);
%     sStatOL = smoothData(responseOL.stateProfiles.Stationary, w);
% 
%     %  REFERENCE NORMALIZATION: Use sOdd to define the scaling factors
%     minOdd = min(sOdd, [], 2);          % Minimum of each ROI in VR Odd
%     maxOdd = max(sOdd, [], 2);          % Maximum of each ROI in VR Odd
%     rangeOdd = maxOdd - minOdd;         % Dynamic range per ROI
% 
%     minEven = min(sEven, [], 2);        % Same as above for even [for figure]
%     maxEven = max(sEven, [], 2);
%     rangeEven = maxEven - minEven;
% 
%     % Avoid division by zero for inactive ROIs
%     rangeOdd(rangeOdd == 0) = 1; 
% 
%     % Apply VR even scaling to Run and Stationary
%     normOdd  = (sOdd - minOdd) ./ rangeOdd;
%     normEven = (sEven - minEven) ./ rangeEven;
%     normRunOL  = (sRunOL - minEven) ./ rangeEven;
%     normStatOL = (sStatOL - minEven) ./ rangeEven;
% 
%     % use peak position in odd to sort 
%     [~, peakPos] = max(normOdd, [], 2);
%     [~, sortIdx] = sort(peakPos);
% 
%     % 
%     hFig = figure('Position', [50 50 1800 950], 'Color', 'w');
%     t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'loose');
% 
%     % 
%     neuralCells = {normOdd, normEven, normRunOL, normStatOL};
%     nTitles = {'Closed-Loop: Odd (Ref)', 'Closed-Loop: Even', 'Open-Loop Replay: Running', 'Open-Loop Replay: Stationary'};
%     for i = 1:4
%         nexttile(i);
%         imagesc(neuralCells{i}(sortIdx, :));
%         formatHeatmap(nTitles{i});
%     end
% 
%     % Row 2: Speed 
%     speedDataRaw = {responseCL.lapPositionRunningSpeed(oddIdx,:), ...
%                     responseCL.lapPositionRunningSpeed(evenIdx,:), ...
%                     responseOL.lapPositionRunningSpeed, ...
%                     responseOL.lapPositionRunningSpeed};
%     sTitles = {'Speed: VR (Odd)', 'Speed: VR (Even)', 'Speed: Replay Run', 'Speed: Replay Stat'};
% 
%     for i = 1:4
%         nexttile(i+4);
%         imagesc(speedDataRaw{i});
%         formatSpeedMap(sTitles{i}, speedLimit, runThresh);
%     end
% 
%     title(t, sprintf('Animal: %s | Signal: %s\nTop: Neural | Bottom: Speed', ...
%         sessionFileInfo.animal_name, signalToUse), 'Interpreter', 'none');
% 
%     % Save
%     figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
%     if ~exist(figSaveDir, 'dir'); mkdir(figSaveDir); end
%     savePath = fullfile(figSaveDir, [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_NeuralBehav_ClosedvsOpenLoop']);
%     exportgraphics(hFig, [savePath '.png'], 'Resolution', 300);
% end
% 
% %% save
% responseCL.stateProfiles.signalUsed = signalToUse;
% responseCL.stateProfiles.stableIdx = stableIdx;
% responseCL.stateProfiles.stabilityMetricUsed = stabilityMetricUsed;
% responseCL.stateProfiles.stableThresh = stableThresh;
% 
% responseOL.stateProfiles.signalUsed = signalToUse;
% responseOL.stateProfiles.stableIdx = stableIdx;
% responseOL.stateProfiles.stabilityMetricUsed = stabilityMetricUsed;
% responseOL.stateProfiles.stableThresh = stableThresh;
% 
% resCL_toSave.stateProfiles = responseCL.stateProfiles;
% resOL_toSave.stateProfiles = responseOL.stateProfiles;
% 
% save(allPaths{idxCL(1)}, '-struct', 'resCL_toSave', '-append');
% save(allPaths{idxOL(1)}, '-struct', 'resOL_toSave', '-append');
% end
% 
% %% helper functions [move to utility]
% 
% function meanAct = getMaskedMean(activity, mask)
%     [nCells, ~, nBins] = size(activity);
%     meanAct = nan(nCells, nBins);
%     for c = 1:nCells
%         for b = 1:nBins
%             validLaps = mask(:, b);
%             if any(validLaps)
%                 meanAct(c, b) = mean(squeeze(activity(c, validLaps, b)), 'omitnan');
%             end
%         end
%     end
% end
% 
% function smoothed = smoothData(data, w)
%     smoothed = data;
%     for i = 1:size(data, 1)
%         trace = data(i, :);
%         nanMask = isnan(trace); trace(nanMask) = 0;
%         s = filtfilt(w, 1, trace); s(nanMask) = NaN;
%         smoothed(i, :) = s;
%     end
% end
% 
% function formatHeatmap(titleStr)
%     clim([0 1]); 
%     colormap(gca, flipud(gray)); cb = colorbar; cb.Label.String = 'Norm. \Delta F/F';
%     set(gca, 'TickDir', 'out', 'YDir', 'normal', 'FontSize', 11, 'Box', 'off');
%     xlabel('Position (cm)'); ylabel('ROIs'); title(titleStr);
%     for x = 40:40:160, xline(x, 'k--', 'LineWidth', 2); end
% end
% 
% % gemini [set <1cm to black]
% function formatSpeedMap(titleStr, speedLimit, thresh)
%     nColors = 256;
%     idxThresh = max(1, round((thresh / speedLimit) * nColors));
% 
%     % Use Parula as the base (the blue/yellow you like)
%     baseMap = parula(nColors);
% 
%     % Force 0 to thresh to be PURE BLACK
%     baseMap(1:idxThresh, :) = 0; 
% 
%     % Force the start of the color range to be a visible Blue/Cyan
%     % This skips the 'dark' part of parula so it doesn't look muddy
%     colorStartIdx = round(0.15 * nColors); % Start 15% into the parula map
%     if idxThresh < nColors
%         % Create a gradient from vibrant blue to yellow for the remaining space
%         newColors = parula(nColors - idxThresh);
%         baseMap(idxThresh+1:end, :) = newColors;
%     end
% 
%     colormap(gca, baseMap);
%     clim([0 speedLimit]);
%     cb = colorbar; cb.Label.String = 'Speed (cm/s)';
%     set(gca, 'TickDir', 'out', 'YDir', 'normal', 'FontSize', 11, 'Box', 'off');
%     xlabel('Position (cm)'); ylabel('Laps'); title(titleStr);
%     for x = 40:40:160, xline(x, 'w--', 'LineWidth', 1.2, 'Alpha', 0.4); end
% end

function [responseCL, responseOL] = getSortedClosedvsOpenReplay_RunningStationary(sessionFileInfo, CLStimName, OLStimName, signalToUse, doPlot)
% function [responseCL, responseOL] = getSortedClosedvsOpenReplay_RunningStationary(sessionFileInfo, CLStimName, OLStimName, signalToUse, doPlot)
% Extracts Running vs Stationary profiles for CL and OL, saving raw means for future alignment.

if nargin < 4; signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 5; doPlot = true; end

runThresh = 1; 
stableThresh = 0.4;

%% 1. Find Paths
allPaths = {sessionFileInfo.stimFiles.Response}; 
% idxCL = find(contains(allPaths, CLStimName));
% idxOL = find(contains(allPaths, OLStimName));

idxCL = find(cellfun(@(x) ischar(x) && contains(x, CLStimName), allPaths));
idxOL = find(cellfun(@(x) ischar(x) && contains(x, OLStimName), allPaths));


if isempty(idxCL); error('CL path not found.'); end
if isempty(idxOL); error('OL path not found.'); end

%% 2. Load Data
L = load(allPaths{idxCL(1)}, 'lapPositionActivity', 'lapPositionRunningSpeed', 'trialIndicesByCondition');
responseCL.lapPositionActivity.(signalToUse) = L.lapPositionActivity.(signalToUse);
responseCL.lapPositionRunningSpeed = L.lapPositionRunningSpeed;

L = load(allPaths{idxOL(1)}, 'lapPositionActivity', 'lapPositionRunningSpeed');
responseOL.lapPositionActivity.(signalToUse) = L.lapPositionActivity.(signalToUse);
responseOL.lapPositionRunningSpeed = L.lapPositionRunningSpeed;
clear L; 

%% 3. Filter for Stable ROIs
stabilityMetricUsed = 'lapCorr_Halves';
try
    vars = load(sessionFileInfo.otherSessFilePaths.sessionROIData, stabilityMetricUsed);
    stableIdx = find(vars.lapCorr_Halves.rho >= stableThresh);
catch
    stableIdx = 1:size(responseCL.lapPositionActivity.(signalToUse), 1);
end

% Subset to stable ROIs for profile calculation
actCL_stable = responseCL.lapPositionActivity.(signalToUse)(stableIdx, :, :);
actOL_stable = responseOL.lapPositionActivity.(signalToUse)(stableIdx, :, :);


%% 4. Masking Running and Stationary (Raw Means)
maskRunCL = responseCL.lapPositionRunningSpeed > runThresh;
maskRunOL = responseOL.lapPositionRunningSpeed > runThresh;

% We save raw means to stateProfiles (Normalization happens in the next function)
responseCL.stateProfiles.Running = getMaskedMean(actCL_stable, maskRunCL);
responseCL.stateProfiles.Stationary = getMaskedMean(actCL_stable, responseCL.lapPositionRunningSpeed <= runThresh);

responseOL.stateProfiles.Running = getMaskedMean(actOL_stable, maskRunOL);
responseOL.stateProfiles.Stationary = getMaskedMean(actOL_stable, responseOL.lapPositionRunningSpeed <= runThresh);

%% 5. Plotting (Internal Audit)
if doPlot
    nLapsCL = size(actCL_stable, 2);
    oddIdx = 1:2:nLapsCL;
    evenIdx = 2:2:nLapsCL;
    w = gausswin(15); w = w / sum(w);
    
    maxSpd = max([max(responseCL.lapPositionRunningSpeed(:)), max(responseOL.lapPositionRunningSpeed(:))]);
    speedLimit = ceil(maxSpd / 5) * 5; 
    
    % Temporary Reference Normalization for THIS figure only
    sOdd  = smoothData(squeeze(mean(actCL_stable(:, oddIdx, :), 2, 'omitnan')), w);
    sEven = smoothData(squeeze(mean(actCL_stable(:, evenIdx, :), 2, 'omitnan')), w);
    sRunOL  = smoothData(responseOL.stateProfiles.Running, w);
    sStatOL = smoothData(responseOL.stateProfiles.Stationary, w);
    
    % Calculate scaling factors for Odd
    minOdd = min(sOdd, [], 2);          
    maxOdd = max(sOdd, [], 2);
    rangeOdd = maxOdd - minOdd;
    rangeOdd(rangeOdd == 0) = 1; % Prevent NaN if Odd is flat

    % Calculate scaling factors for Even
    minEven = min(sEven, [], 2);      
    maxEven = max(sEven, [], 2);
    rangeEven = maxEven - minEven;
    rangeEven(rangeEven == 0) = 1; % Prevent NaN if Even is flat

    %
    % Odd is normalized to its own Max/Min (for clean sorting)
    normOdd  = (sOdd - minOdd) ./ rangeOdd;
    
    % Even is normalized to its own Max/Min
    normEven = (sEven - minEven) ./ rangeEven;
    
    %  OL States are scaled to the EVEN range 
    normRunOL  = (sRunOL - minEven) ./ rangeEven;
    normStatOL = (sStatOL - minEven) ./ rangeEven;

    % Use peak position in ODD to define the ROI sorting order
    [~, peakPos] = max(normOdd, [], 2);
    [~, sortIdx] = sort(peakPos);

    hFig = figure('Position', [50 50 1800 950], 'Color', 'w');
    t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'loose');
    
    % Row 1: Neural
    neuralCells = {normOdd, normEven, normRunOL, normStatOL};
    nTitles = {'CL: Odd (Ref)', 'CL: Even', 'OL Replay: Running', 'OL Replay: Stationary'};
    for i = 1:4
        nexttile(i);
        imagesc(neuralCells{i}(sortIdx, :));
        formatHeatmap(nTitles{i});
    end
    
    % Row 2: Speed
    speedDataRaw = {responseCL.lapPositionRunningSpeed(oddIdx,:), ...
                    responseCL.lapPositionRunningSpeed(evenIdx,:), ...
                    responseOL.lapPositionRunningSpeed, ...
                    responseOL.lapPositionRunningSpeed};
    sTitles = {'Speed: CL (Odd)', 'Speed: CL (Even)', 'Speed: Replay Run', 'Speed: Replay Stat'};
    
    for i = 1:4
        nexttile(i+4);
        imagesc(speedDataRaw{i});
        formatSpeedMap(sTitles{i}, speedLimit, runThresh);
    end
    
    title(t, sprintf('Animal: %s | Signal: %s\n Neural (Ref-Norm to CL Odd) | Bottom: Speed', ...
        sessionFileInfo.animal_name, signalToUse), 'Interpreter', 'none');
    
    figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
    if ~exist(figSaveDir, 'dir'); mkdir(figSaveDir); end
    savePath = fullfile(figSaveDir, [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_NeuralBehav_OpenvsCLosedLoop']);
    exportgraphics(hFig, [savePath '.png'], 'Resolution', 300);
end

%% 6. Save Append
% We save the metadata and the raw profiles
responseCL.stateProfiles.signalUsed = signalToUse;
responseCL.stateProfiles.stableIdx = stableIdx;
responseCL.stateProfiles.stabilityMetricUsed = stabilityMetricUsed;
responseCL.stateProfiles.stableThresh = stableThresh;

responseOL.stateProfiles.signalUsed = signalToUse;
responseOL.stateProfiles.stableIdx = stableIdx;
responseOL.stateProfiles.stabilityMetricUsed = stabilityMetricUsed;
responseOL.stateProfiles.stableThresh = stableThresh;

resCL_toSave.stateProfiles = responseCL.stateProfiles;
resOL_toSave.stateProfiles = responseOL.stateProfiles;

save(allPaths{idxCL(1)}, '-struct', 'resCL_toSave', '-append');
save(allPaths{idxOL(1)}, '-struct', 'resOL_toSave', '-append');
end

%% --- Helper Functions ---
function meanAct = getMaskedMean(activity, mask)
    [nCells, ~, nBins] = size(activity);
    meanAct = nan(nCells, nBins);
    for c = 1:nCells
        for b = 1:nBins
            validLaps = mask(:, b);
            if any(validLaps)
                meanAct(c, b) = mean(squeeze(activity(c, validLaps, b)), 'omitnan');
            end
        end
    end
end

function smoothed = smoothData(data, w)
    smoothed = data;
    for i = 1:size(data, 1)
        trace = data(i, :);
        nanMask = isnan(trace); trace(nanMask) = 0;
        s = filtfilt(w, 1, trace); s(nanMask) = NaN;
        smoothed(i, :) = s;
    end
end

function formatHeatmap(titleStr)
    clim([0 1]); colormap(gca, flipud(gray)); 
    cb = colorbar; cb.Label.String = 'Norm. \Delta F/F';
    set(gca, 'TickDir', 'out', 'YDir', 'normal', 'FontSize', 11, 'Box', 'off');
    xlabel('Position (cm)'); ylabel('ROIs'); title(titleStr);
    for x = 40:40:160, xline(x, 'k--', 'LineWidth', 2); end
end

function formatSpeedMap(titleStr, speedLimit, thresh)
    nColors = 256;
    idxThresh = max(1, round((thresh / speedLimit) * nColors));
    baseMap = parula(nColors);
    baseMap(1:idxThresh, :) = 0; 
    if idxThresh < nColors
        newColors = parula(nColors - idxThresh);
        baseMap(idxThresh+1:end, :) = newColors;
    end
    colormap(gca, baseMap); clim([0 speedLimit]);
    cb = colorbar; cb.Label.String = 'Speed (cm/s)';
    set(gca, 'TickDir', 'out', 'YDir', 'normal', 'FontSize', 11, 'Box', 'off');
    xlabel('Position (cm)'); ylabel('Laps'); title(titleStr);
    for x = 40:40:160, xline(x, 'w--', 'LineWidth', 1.2, 'Alpha', 0.4); end
end