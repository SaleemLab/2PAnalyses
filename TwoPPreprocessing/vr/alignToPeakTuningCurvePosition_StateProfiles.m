function [alignedFull, response] = alignToPeakTuningCurvePosition_StateProfiles(sessionFileInfo, CLStimName, OLStimName, plotFlag)
if nargin < 4, plotFlag = 1; end
w = gausswin(15); w = w / sum(w);

%% 1. Find and Load Data
allPaths = {sessionFileInfo.stimFiles.Response};
idxCL = find(contains(allPaths, CLStimName));
idxOL = find(contains(allPaths, OLStimName));
if isempty(idxCL) || isempty(idxOL); error('CL or OL paths not found.'); end

L_CL = load(allPaths{idxCL(1)}, 'lapPositionActivity', 'trialIndicesByCondition', 'stateProfiles');
L_OL = load(allPaths{idxOL(1)}, 'stateProfiles');
numBins = size(L_CL.stateProfiles.Running, 2);
stableIdx = L_CL.stateProfiles.stableIdx; 
numStable = length(stableIdx);

%% 2. Reconstruct Reference (CL Odd) and Scaling (CL Even)
actCL = L_CL.lapPositionActivity.(L_CL.stateProfiles.signalUsed)(stableIdx, :, :);
conds = fieldnames(L_CL.trialIndicesByCondition);
baseIdx = find(contains(lower(conds), 'baseline') | contains(lower(conds), 'norm'), 1);
if isempty(baseIdx), baseIdx = 1; end
baseLaps = L_CL.trialIndicesByCondition.(conds{baseIdx});
oddLaps = baseLaps(1:2:end);
evenLaps = baseLaps(2:2:end);

sOdd  = smoothData(squeeze(mean(actCL(:, oddLaps, :), 2, 'omitnan')), w);
sEven = smoothData(squeeze(mean(actCL(:, evenLaps, :), 2, 'omitnan')), w);
sRunOL  = smoothData(L_OL.stateProfiles.Running, w);
sStatOL = smoothData(L_OL.stateProfiles.Stationary, w);

%% 3. Normalization and Peak Detection (With Constraints)
startBin = 10;
endBin = 190;

% Reference Odd for sorting
minOdd = min(sOdd, [], 2); maxOdd = max(sOdd, [], 2); rangeOdd = maxOdd - minOdd;
rangeOdd(rangeOdd == 0) = 1;
normOdd  = (sOdd - minOdd) ./ rangeOdd;
normEven = (sEven - minOdd) ./ rangeOdd;
normRunOL  = (sRunOL - minOdd) ./ rangeOdd;
normStatOL = (sStatOL - minOdd) ./ rangeOdd;

% Anchor Peaks constrained to startBin/endBin
anchorPeaks = zeros(numStable, 1);
for i = 1:numStable
    [~, p] = max(normOdd(i, startBin:endBin));
    anchorPeaks(i) = p + (startBin - 1);
end

%% 4. Alignment
dataCells = {normEven, normRunOL, normStatOL};
titles = {'CL_Even_SelfRef', 'OL_Running', 'OL_Stationary'};
displayTitles = {'CL: Even', 'OL: Running', 'OL: Stationary'};
centerIdx = numBins + 1;
alignedFull = struct();

for s = 1:3
    thisData = dataCells{s};
    alignedData = nan(numStable, numBins * 2 + 1);
    
    for i = 1:numStable
        p = anchorPeaks(i);
        targetStart = centerIdx - (p - 1);
        targetEnd = targetStart + numBins - 1;
        % Ensure indices are valid to avoid out-of-bounds errors
        if targetStart >= 1 && targetEnd <= size(alignedData, 2)
            alignedData(i, targetStart:targetEnd) = thisData(i, :);
        end
    end
    
    alignedFull.(titles{s}).data = alignedData;
    alignedFull.(titles{s}).anchorPeaks = anchorPeaks;
    
    if plotFlag
        fig = figure('Color', 'w', 'Position', [100 100 550 800]);
        fullShiftBins = (-numBins : numBins);
        [~, snakeSortIdx] = sort(anchorPeaks, 'descend');
        plotMat = alignedData(snakeSortIdx, :);
        
        subplot('Position', [0.15 0.45 0.7 0.45]);
        h = imagesc(fullShiftBins, 1:numStable, plotMat);
        set(h, 'AlphaData', ~isnan(plotMat)); 
        colormap(flipud(gray));
        clim([0 1]); 
        hold on;
        xline(-80, 'r--', 'LineWidth', 1.5, 'Alpha', 0.5);
        xline(80, 'b--', 'LineWidth', 1.5, 'Alpha', 0.5);
        xline(0, 'y-', 'LineWidth', 2);
        xlim([-100, 100]); title(displayTitles{s});
        axis off;
        
        subplot('Position', [0.15 0.12 0.7 0.25]);
        popMean = mean(alignedData, 1, 'omitnan');
        popSEM = std(alignedData, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(alignedData), 1));
        fill([fullShiftBins, fliplr(fullShiftBins)], [popMean+popSEM, fliplr(popMean-popSEM)], ...
             [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        hold on; plot(fullShiftBins, popMean, 'k', 'LineWidth', 2);
        xline(-80, 'r--', 'LineWidth', 1.5); xline(80, 'b--', 'LineWidth', 1.5); xline(0, 'y-', 'LineWidth', 2);
        xlabel('Dist. from CL Master Peak (cm)'); ylabel('Mean Resp.');
        xlim([-100, 100]); xticks([-80, -40, 0, 40, 80]);
        set(gca, 'Box', 'off', 'TickDir', 'out');
        
        figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
        if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
        exportgraphics(fig, fullfile(figSaveDir, sprintf('%s_%s_%s_MasterAligned.png', ...
            sessionFileInfo.animal_name, sessionFileInfo.session_name, titles{s})), 'Resolution', 300);
        close(fig);
    end
end
%% Save back
response = L_OL.stateProfiles; 
response.peakAlignedToMaster = alignedFull;
stateProfiles = response; 
save(allPaths{idxOL(1)}, 'stateProfiles', '-append');
disp(['Master aligned profiles saved to: ', allPaths{idxOL(1)}]);
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


% function [alignedFull, response] = alignToPeakTuningCurvePosition_StateProfiles(sessionFileInfo, VRStimName, plotFlag)
% % alignToPeakTuningCurvePosition_StateProfiles - Aligns to independent peaks using original plotting style.
% if nargin < 3, plotFlag = 1; end
% 
% %% Load data
% stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
% if isempty(stimIdx), error('specified vrstimname not found.'); end
% 
% responsefilePath = sessionFileInfo.stimFiles(stimIdx).Response;
% data = load(responsefilePath, 'stateProfiles');
% response = data.stateProfiles; 
% 
% % Matrix dimensions
% numBins = size(response.Running, 2); 
% centerIdx = numBins + 1;
% states = {'Running', 'Stationary'};
% alignedFull = struct();
% 
% %% Align independently 
% for s = 1:length(states)
%     thisState = states{s};
%     rawProfile = response.(thisState); 
%     nROIs = size(rawProfile, 1);
% 
%     % Find independent peaks for this specific state
%     [~, statePeaks] = max(rawProfile, [], 2);
% 
%     % Initialize padded matrix
%     alignedData = nan(nROIs, numBins * 2 + 1);
% 
%     for i = 1:nROIs
%         peakPos = statePeaks(i);
%         profile = rawProfile(i, :);
% 
%         targetStart = centerIdx - (peakPos - 1);
%         targetEnd = targetStart + numBins - 1;
% 
%         alignedData(i, targetStart:targetEnd) = profile;
%     end
% 
%     alignedFull.(thisState).data = alignedData;
%     alignedFull.(thisState).peaks = statePeaks;
% 
%     %% Plotting (Original Style)
%     if plotFlag
%         fig = figure('Color', 'w', 'Position', [100 100 550 800]);
%         fullShiftBins = (-numBins : numBins);
%         [~, snakeSortIdx] = sort(statePeaks, 'descend');
% 
%         % 1. Heatmap
%         subplot('Position', [0.15 0.45 0.7 0.45]);
%         imagesc(fullShiftBins, 1:nROIs, alignedData(snakeSortIdx, :));
%         colormap(flipud(gray));
%         clim([0, 1]); 
% 
%         hold on;
%         xline(-80, 'r', 'LineWidth', 1.5);
%         xline(80, 'b', 'LineWidth', 1.5);
%         xline(0, 'y', 'LineWidth', 1.5);
% 
%         xlim([-100, 100]); 
%         title(sprintf('%s Aligned Population (n=%d)', thisState, nROIs));
%         axis off;
% 
%         % 2. Mean Profile
%         subplot('Position', [0.15 0.12 0.7 0.25]);
%         popMean = mean(alignedData, 1, 'omitnan');
%         numSamples = sum(~isnan(alignedData), 1);
%         popSEM = std(alignedData, 0, 1, 'omitnan') ./ sqrt(numSamples);
% 
%         fill([fullShiftBins, fliplr(fullShiftBins)], [popMean+popSEM, fliplr(popMean-popSEM)], ...
%              [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%         hold on;
%         plot(fullShiftBins, popMean, 'k', 'LineWidth', 2);
% 
%         xline(-80, 'r', 'LineWidth', 1.5);
%         xline(80, 'b', 'LineWidth', 1.5);
%         xline(0, 'y', 'LineWidth', 1.5);
% 
%         xlabel('Distance from peak (cm)');
%         ylabel('Mean \deltaF/F');
%         xlim([-100, 100]); 
%         xticks([-80, -40, 0, 40, 80]);
%         ylim([min(popMean)*0.9, max(popMean)*1.2]);
% 
%         set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 10);
% 
%         % Save
%         figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
%         if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
%         filename = fullfile(figSaveDir, sprintf('%s_%s_%s_%s_AlignedPeak_%s.png', ...
%             sessionFileInfo.animal_name, sessionFileInfo.session_name, response.signalUsed, VRStimName, thisState));
%         exportgraphics(fig, filename, 'Resolution', 300);
%     end
% end
% 
% %% Save to file
% response.peakAlignedIndependent = alignedFull;
% stateProfiles = response; % Restore variable name for storage
% save(responsefilePath, 'stateProfiles', '-append');
% disp(['Saving independent aligned tuning curves to ', responsefilePath]);
% 
% end