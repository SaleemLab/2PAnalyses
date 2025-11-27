function [R_sorted, normOdd, normEven, sortIdx, numCells, dataScope, signalToUse] = SpatialTuningAndBoutonCorrelationMatrix(mouseID, session, targetArea, checkVRRecOnly, checkNonVRRecOnly, onlyIncludeROIs, signalToUse, applySmoothing, plotFlag)
%   Plots a 3-panel figure:
%   1. Pairwise ROI correlation matrix
%   2. Sorted odd lap spatial tuning heatmap
%   3. Sorted even lap spatial tuning heatmap
%
%   All three panels are sorted based on the peak response location from the
%   odd laps (Panel 2).
%% Default Arguments
if nargin < 6; onlyIncludeROIs = false; end
if nargin < 7; signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 8; applySmoothing = true; end
if nargin < 9; plotFlag = false; end
% if nargin < 9; plotFlag = true; end
%% Load SessionFileInfo
rootDir = 'Z:\ibn-vision\DATA\SUBJECTS';
fileName = sprintf('%s_%s_sessionFileInfo.mat', mouseID, session);
sessionFilePath = fullfile(rootDir, mouseID, 'Analysis', session, fileName);
if ~exist(sessionFilePath, 'file')
    fprintf('Session file not found: %s\n', sessionFilePath);
    return;
end
fprintf('Loading session file: %s\n', sessionFilePath);
load(sessionFilePath);
% sessionFileInfo is now loaded
%% Get F and ROIs
if checkVRRecOnly
    fprintf('checkVRRecOnly=true. Searching for VRCorr stim files...\n');
    isVRStim = find( contains({sessionFileInfo.stimFiles.name}, 'VRCorr') & ~contains({sessionFileInfo.stimFiles.name}, 'CombinedRuns') );
    if isempty(isVRStim)
        fprintf('Error: No VRCorr stim file found.\n');
        return;
    elseif length(isVRStim) > 1 % This is a short-term fix
        fprintf('Multiple VRCorr files found. Selecting last one (Index: %d).\n', isVRStim(end));
        isVRStim = isVRStim(end); % Select the last one
    else
        fprintf('Found VRCorr file (Index: %d).\n', isVRStim);
    end
    % Load for VR stimulus
    twoPDataPath = sessionFileInfo.stimFiles(isVRStim).mergedBonsai2PSuite2pData;
    fprintf('Loading 2P data file: %s\n', twoPDataPath);
    load(twoPDataPath); % Loads twoPData
    isROI = twoPData.iscell;
    F = twoPData.F;
    dataScope = 'VRStimulusOnly';

elseif checkNonVRRecOnly
    fprintf('checkNonVRRecOnly=true. Searching for stim files excluding VRCorr\n');
    isNotVRStim = find(~contains({sessionFileInfo.stimFiles.name}, 'VRCorr'));
    if length(isNotVRStim) > 1
        % Load the first TwoPData
        firstTwoPDataPath = sessionFileInfo.stimFiles(isNotVRStim(1)).mergedBonsai2PSuite2pData;
        firstTwoPData = load(firstTwoPDataPath);
        % Initialise the F matrix and isROI array
        F = firstTwoPData.twoPData.F;
        isROI = firstTwoPData.twoPData.iscell;
        % Now loop through the *remaining* stim files (starting from the 2nd)
        for thisOtherStim = 2:length(isNotVRStim)
            thisStimIdx = isNotVRStim(thisOtherStim);
            twoPDataPath = sessionFileInfo.stimFiles(thisStimIdx).mergedBonsai2PSuite2pData;
            twopdata = load(twoPDataPath);
            % Concatenate the F (neuron x time) across the time vector
            F = [F, twopdata.twoPData.F];
            dataScope = 'NonVRStimOnly';
        end
    end
else
    fprintf('checkVRRecOnly=false. Finding and loading Fall.mat...\n');
    FallPath = findFile(sessionFileInfo.suite2pFiles.planes, 'Fall.mat');
    fprintf('Loading Fall.mat file: %s\n', FallPath);
    fAll = load(FallPath); % Loads a struct, e.g., fAll.F, fAll.iscell
    isROI = fAll.iscell;
    F = fAll.F;
    dataScope = 'ConcatenatedStimuli';
end
% Select ROIs based on isROI and the onlyIncludeROIs flag
if onlyIncludeROIs
    ROIs = find(isROI(:, 1) == 1);
    fprintf('Including %d ROIs (iscell == 1)\n', length(ROIs));
else
    % Include all potential ROIs
    ROIs = 1:size(isROI, 1);
    fprintf('Including all %d potential ROIs\n', length(ROIs));
end
numCells = length(ROIs);
if numCells < 2
    fprintf('Fewer than 2 ROIs found. Cannot proceed.\n');
    return;
end
%% Correlation matrix
% Extract F for selected ROIs
fSelected = F(ROIs, :);
% Calculate Pairwise Correlations
fprintf('Calculating %d x %d correlation matrix...\n', numCells, numCells);
R = corrcoef(fSelected'); % transpose!!
% R is now numCells x numCells, corresponding to the rois list
%% Calculate spatial tuning curves
% Find 'CombinedRuns' first
isCombined = find(contains({sessionFileInfo.stimFiles.name}, 'CombinedRuns') );
if ~isempty(isCombined)
    if length(isCombined) > 1
        fprintf('Multiple CombinedRuns files found. Selecting last one (Index: %d).\n', isCombined(end));
        stimFileIndex = isCombined(end);
    else
        fprintf('Found CombinedRuns file (Index: %d).\n', isCombined);
        stimFileIndex = isCombined;
    end
else
    % If no CombinedRuns, look for 'VRCorr'
    fprintf('No CombinedRuns file found. Searching for VRCorr stim files...\n');
    isVRCorr = find( contains({sessionFileInfo.stimFiles.name}, 'VRCorr') & ~contains({sessionFileInfo.stimFiles.name}, 'CombinedRuns') );
    if isempty(isVRCorr)
        fprintf('Error: No CombinedRuns or VRCorr stim file found to load response.mat from.\n');
        return;
    elseif length(isVRCorr) > 1
        fprintf('Multiple VRCorr files found. Selecting last one (Index: %d).\n', isVRCorr(end));
        stimFileIndex = isVRCorr(end);
    else
        fprintf('Found VRCorr file (Index: %d).\n', isVRCorr);
        stimFileIndex = isVRCorr;
    end
end
responseFilePath = fullfile(sessionFileInfo.stimFiles(stimFileIndex).Response);
if ~exist(responseFilePath, 'file')
    fprintf('Response file not found: %s\n', responseFilePath);
    fprintf('Please check the loading path in SpatialTuningAndBoutonCorrelationMatrix.\n');
    return;
end
fprintf('Loading response file: %s\n', responseFilePath);
load(responseFilePath, 'response');
try
    lapActivityAll = response.lapPositionActivity.(signalToUse);
catch
    fprintf('Error: signalToUse "%s" not found in response.lapPositionActivity.\n', signalToUse);
    return;
end
% Filter the lapActivity matrix to match the exact ROIs used for the correlation matrix
if max(ROIs) > size(lapActivityAll, 1)
    fprintf('CRITICAL ERROR: Mismatch in ROI count.\n');
    fprintf('iscell has %d rows, but response.mat has %d rows.\n', size(isROI, 1), size(lapActivityAll, 1));
    return;
end
lapActivity = lapActivityAll(ROIs, :, :);
fprintf('Filtered lap activity to %d matched ROIs.\n', numCells);
if applySmoothing
    w = gausswin(9); w = w / sum(w);
    for iCell = 1:size(lapActivity, 1)
        for iLap = 1:size(lapActivity, 2)
            trace = squeeze(lapActivity(iCell, iLap, :));
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace);
            trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace);
            smoothed(nanMask) = NaN;
            lapActivity(iCell, iLap, :) = smoothed;
        end
    end
end
% Split odd and even laps
oddLaps = lapActivity(:, 1:2:end, :);
evenLaps = lapActivity(:, 2:2:end, :);
% Average across laps
meanOdd = squeeze(mean(oddLaps, 2, 'omitnan'));
meanEven = squeeze(mean(evenLaps, 2, 'omitnan'));
% Normalize across position bins
normOdd = normalize(meanOdd, 2, 'range');
normEven = normalize(meanEven, 2, 'range');
% Sort cells by peak location in odd lap average
[~, peakIdx] = max(normOdd, [], 2);
[~, sortIdx] = sort(peakIdx);
%% Apply sorting to correlation matrix
% Re-order the correlation matrix 'R' using the 'sortIdx'
R_sorted = R(sortIdx, sortIdx);

%% 3-pannel figure
if plotFlag
    hFig = figure('Position', [100 100 1600 600]);

    % Sorted Correlation Matrix
    ax1 = subplot(1, 3, 1);
    imagesc(R_sorted);
    colormap(ax1, 'jet');
    caxis([-0.2 1]);
    axis square;
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 10);
    xlabel('Sorted ROIs');
    ylabel('Sorted ROIs');
    title(sprintf('Correlation Matrix (Sorted)\n%s', dataScope));
    colorbar; ylabel(colorbar, 'Pearson (r)');

    % Odd Laps (The sorting reference)
    ax2 = subplot(1, 3, 2);
    imagesc(normOdd(sortIdx, :));
    caxis([0 1]); colormap(ax2, flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 10, 'YDir', 'normal');
    xline(50, 'k--', 'LineWidth', 1.5);
    xline(70, 'k--', 'LineWidth', 1.5);
    xline(90, 'k--', 'LineWidth', 1.5);
    xline(110, 'k--', 'LineWidth', 1.5);
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    xlabel('Position (cm)');
    ylabel('ROIs (sorted by odd peak)');
    title(sprintf('Odd Laps (Sorted)\n%s', signalToUse));
    colorbar; ylabel(colorbar, 'Activity (normalised)');

    % Even Laps (Sorted by Odd)
    ax3 = subplot(1, 3, 3);
    imagesc(normEven(sortIdx, :));
    caxis([0 1]); colormap(ax3, flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 10, 'YDir', 'normal');
    xline(50, 'k--', 'LineWidth', 1.5);
    xline(70, 'k--', 'LineWidth', 1.5);
    xline(90, 'k--', 'LineWidth', 1.5);
    xline(110, 'k--', 'LineWidth', 1.5);
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    xlabel('Position (cm)');
    ylabel('ROIs (sorted by odd peak)');
    title(sprintf('Even Laps (Sorted by Odd)\n%s', signalToUse));
    colorbar; ylabel(colorbar, 'Activity (normalised)');

    sgtitle(sprintf('%s :: Mouse: %s, Session: %s (n=%d ROIs)', ...
        targetArea, mouseID, session, numCells), 'Interpreter', 'none');

    %% Save Figure
    fprintf('Saving figure...\n');
    savePath = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
    if ~exist(savePath, 'dir')
        mkdir(savePath);
        fprintf('Created save directory: %s\n', savePath);
    end

    fileName = [mouseID '_' session '_TuningAndCorrelationMatrix_' dataScope '.png'];
    fullSavePath = fullfile(savePath, fileName);

    % Use print for saving
    set(hFig, 'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 16 6], ...
        'PaperOrientation', 'landscape');
    print(hFig, fullSavePath, '-dpng', '-r300');

    fprintf('Figure saved to: %s\n', fullSavePath);
end
end