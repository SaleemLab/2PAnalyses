function [responseCL, responseOP] = getSortedClosedvsOpenReplay_RunningStationaryV2(sessionFileInfo, CLStimName, OLStimName, signalToUse, doPlot)
if nargin < 4; signalToUse = 'spks'; end
if nargin < 5; doPlot = true; end
runThresh = 1; 
stableThresh = 0.6;

allPaths = {sessionFileInfo.stimFiles.Response}; 
idxCL = find(contains(allPaths, CLStimName));
idxOP = find(contains(allPaths, OLStimName));
if isempty(idxCL); error('CL path not found.'); end
if isempty(idxOP); error('OL path not found.'); end

L = load(allPaths{idxCL(1)}, 'lapPositionActivity', 'lapPositionRunningSpeed');
responseCL.lapPositionActivity.(signalToUse) = L.lapPositionActivity.(signalToUse);
responseCL.lapPositionRunningSpeed = L.lapPositionRunningSpeed;
L = load(allPaths{idxOP(1)}, 'lapPositionActivity', 'lapPositionRunningSpeed');
responseOP.lapPositionActivity.(signalToUse) = L.lapPositionActivity.(signalToUse);
responseOP.lapPositionRunningSpeed = L.lapPositionRunningSpeed;
clear L; 

try
    vars = load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'lapCorr_Halves');
    stableIdx = find(vars.lapCorr_Halves.rho >= stableThresh);
catch
    stableIdx = 1:size(responseCL.lapPositionActivity.(signalToUse), 1);
end

actCL_stable = responseCL.lapPositionActivity.(signalToUse)(stableIdx, :, :);
actOP_stable = responseOP.lapPositionActivity.(signalToUse)(stableIdx, :, :);

responseCL.stateProfiles.Running = getMaskedMean(actCL_stable, responseCL.lapPositionRunningSpeed > runThresh);
responseCL.stateProfiles.Stationary = getMaskedMean(actCL_stable, responseCL.lapPositionRunningSpeed <= runThresh);
responseCL.stateProfiles.RawMean = squeeze(mean(actCL_stable, 2, 'omitnan'));

responseOP.stateProfiles.Running = getMaskedMean(actOP_stable, responseOP.lapPositionRunningSpeed > runThresh);
responseOP.stateProfiles.Stationary = getMaskedMean(actOP_stable, responseOP.lapPositionRunningSpeed <= runThresh);

if doPlot
    nLapsCL = size(actCL_stable, 2);
    oddIdx = 1:2:nLapsCL;
    w = gausswin(15); w = w / sum(w);
    
    maskRunCL_Odd = responseCL.lapPositionRunningSpeed(oddIdx, :) > runThresh;
    meanRunCL_Odd = getMaskedMean(actCL_stable(:, oddIdx, :), maskRunCL_Odd);
    
    normOddRef = normalize(smoothData(meanRunCL_Odd, w), 2, 'range');
    [~, peakPos] = max(normOddRef, [], 2);
    [~, sortIdx] = sort(peakPos);
    
    normRunCL = normalize(smoothData(responseCL.stateProfiles.Running, w), 2, 'range');
    normRunOP = normalize(smoothData(responseOP.stateProfiles.Running, w), 2, 'range');
    normStatCL = normalize(smoothData(responseCL.stateProfiles.Stationary, w), 2, 'range');
    normStatOP = normalize(smoothData(responseOP.stateProfiles.Stationary, w), 2, 'range');
    
    hFig = figure('Position', [50 100 1200 900], 'Color', 'w');
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
    
    % Data ordered for: [Run VR, Run Replay; Stat VR, Stat Replay]
    dataCell = {normRunCL, normRunOP, normStatCL, normStatOP};
    titles = {'Running: VR', 'Running: Replay', 'Stationary: VR', 'Stationary: Replay'};
    
    for i = 1:4
        nexttile;
        imagesc(dataCell{i}(sortIdx, :));
        formatHeatmap(titles{i});
    end
    
    title(t, sprintf('Animal: %s | Signal: %s\nRows: State | Cols: Loop Condition', ...
        sessionFileInfo.animal_name, signalToUse), 'Interpreter', 'none');
    
    figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
    if ~exist(figSaveDir, 'dir'); mkdir(figSaveDir); end
    savePath = fullfile(figSaveDir, [sessionFileInfo.animal_name '_' OLStimName '_StateComparison_Heatmap']);
    exportgraphics(hFig, [savePath '.png'], 'Resolution', 300);
end

resCL_toSave.stateProfiles = responseCL.stateProfiles;
resCL_toSave.stableIdx = stableIdx;
resOP_toSave.stateProfiles = responseOP.stateProfiles;
resOP_toSave.stableIdx = stableIdx;

disp(['Appending stateProfiles to ', allPaths{idxCL(1)}]);
save(allPaths{idxCL(1)}, '-struct', 'resCL_toSave', '-append');

disp(['Appending stateProfiles to ', allPaths{idxOP(1)}]);
save(allPaths{idxOP(1)}, '-struct', 'resOP_toSave', '-append');

end

function meanAct = getMaskedMean(activity, mask)
    [nCells, ~, nBins] = size(activity);
    meanAct = nan(nCells, nBins);
    for c = 1:nCells
        for b = 1:nBins
            data = squeeze(activity(c, :, b));
            validLaps = mask(:, b);
            if any(validLaps)
                meanAct(c, b) = mean(data(validLaps), 'omitnan');
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
    clim([0 1]); colormap(flipud(gray));
    set(gca, 'TickDir', 'out', 'YDir', 'normal', 'FontSize', 9, 'Box', 'off');
    xlabel('Position (cm)'); ylabel('ROIs'); title(titleStr);
    for x = 40:40:160, xline(x, 'k--', 'LineWidth', 2); end
    xticks([1 40 80 120 160 200]);
end