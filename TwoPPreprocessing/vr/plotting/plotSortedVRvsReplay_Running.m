function plotSortedVRvsReplay_Running(sessionFileInfo, responseCL, responseOP, signalToUse)
% Plots: VR (Raw) vs Replay (Run/Stat) for Halves-Stable ROIs only.
% Sorted by VR Raw Odd laps.

if nargin < 4; signalToUse = 'dFFNeuropilCorrected'; end
runThresh = 1; 
stableThresh = 0.6;

%% 1. Load and Filter for Stable ROIs
% We use the CL session to determine which ROIs are stable
try
    vars = load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'lapCorr_Halves');
    stableIdx = find(vars.lapCorr_Halves.rho >= stableThresh);
    fprintf('Found %d stable ROIs (rho >= %.2f)\n', length(stableIdx), stableThresh);
catch
    stableIdx = 1:size(responseCL.lapPositionActivity.(signalToUse), 1);
    warning('Stability data not found. Using all %d ROIs.', length(stableIdx));
end

% Filter both CL and OP data matrices by these stable indices
actCL = responseCL.lapPositionActivity.(signalToUse)(stableIdx, :, :);
actOP = responseOP.lapPositionActivity.(signalToUse)(stableIdx, :, :);

%% 2. Split Laps and Create Speed Masks
nLapsCL = size(actCL, 2);
oddIdx = 1:2:nLapsCL;
evenIdx = 2:2:nLapsCL;

% Replay Speed Masks
maskRunOP = responseOP.lapPositionRunningSpeed > runThresh;
maskStatOP = responseOP.lapPositionRunningSpeed <= runThresh;

%% 3. Calculate Means
% VR: Raw (No speed mask, just mean of stable ROIs)
meanRawCL_Odd  = squeeze(mean(actCL(:, oddIdx, :), 2, 'omitnan'));
meanRawCL_Even = squeeze(mean(actCL(:, evenIdx, :), 2, 'omitnan'));

% Replay: Split by Speed (Masked mean of stable ROIs)
meanRunOP  = getMaskedMean(actOP, maskRunOP);
meanStatOP = getMaskedMean(actOP, maskStatOP);

%% 4. Smoothing, Sorting & Normalization
w = gausswin(15); w = w / sum(w);

% Reference: Raw VR Odd (Stable Only)
normOdd = normalize(smoothData(meanRawCL_Odd, w), 2, 'range');
[~, peakPos] = max(normOdd, [], 2);
[~, sortIdx] = sort(peakPos);

% Apply to others
normEvenCL = normalize(smoothData(meanRawCL_Even, w), 2, 'range');
normRunOP  = normalize(smoothData(meanRunOP, w), 2, 'range');
normStatOP = normalize(smoothData(meanStatOP, w), 2, 'range');

%% 5. Plotting (4 Panels)
hFig = figure('Position', [50 100 1600 500], 'Color', 'w');
t = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'loose');

% Panel 1: VR Odd (Ref)
nexttile; imagesc(normOdd(sortIdx, :)); 
formatHeatmap('Closed-Loop Odd');

% Panel 2: VR Even
nexttile; imagesc(normEvenCL(sortIdx, :)); 
formatHeatmap('Closed-Loop Even');

% Panel 3: Replay Running
nexttile; imagesc(normRunOP(sortIdx, :)); 
formatHeatmap('Open-Loop Replay: Running');

% Panel 4: Replay Stationary
nexttile; imagesc(normStatOP(sortIdx, :)); 
formatHeatmap('Open-Loop Replay: Stationary');

title(t, sprintf('%s: Stable ROIs (rho >= %.1f) | Sorted by VR Raw Odd', ...
    sessionFileInfo.animal_name, stableThresh), 'Interpreter', 'none');

%% 6. Save Logic
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
pngName = fullfile(figSaveDir, [sessionFileInfo.animal_name '_ClosedVsOpenLoop_RunningStationaryFramesAcrossPositionBins.png']);

exportgraphics(hFig, pngName, 'Resolution', 300);
savefig(hFig, strrep(pngName, '.png', '.fig'));
end

%% --- Helper Functions ---

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
    nanMask = isnan(trace);
    trace(nanMask) = 0;
    s = filtfilt(w, 1, trace);
    s(nanMask) = NaN;
    smoothed(i, :) = s;
end
end

function formatHeatmap(titleStr)
clim([0 1]); colormap(flipud(gray));
set(gca, 'TickDir', 'out', 'YDir', 'normal', 'FontSize', 9, 'Box', 'off');
xlabel('Position (cm)'); ylabel('Stable ROIs');
title(titleStr);
xline(40, 'k--', 'LineWidth', 2.5);
xline(80, 'k--', 'LineWidth', 2.5);
xline(120, 'k--', 'LineWidth', 2.5);
xline(160, 'k--', 'LineWidth', 2.5);
xticks([1 40 80 120 160]);
xticklabels({'1', '40', '80', '120', '160', '200'});
end