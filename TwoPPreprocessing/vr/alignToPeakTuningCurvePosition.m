function [alignedFull, response] = alignToPeakTuningCurvePosition(sessionFileInfo, VRStimName, signalToUse, plotFlag)
if nargin < 3, signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 4, plotFlag = 1; end
stableThreshold = 0.7; 
%% Load Data
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('specified vrstimname not found.'); end
responsefilePath = sessionFileInfo.stimFiles(stimIdx).Response;
sessionROIDataFilePath = sessionFileInfo.otherSessFilePaths.sessionROIData;
response = load(responsefilePath, 'lapPositionActivity', 'trialIndicesByCondition');
lapActivity = response.lapPositionActivity.(signalToUse);
conds = fieldnames(response.trialIndicesByCondition);
[numROIs, numLaps, numPosBins] = size(lapActivity);
% Pick stable ROIs using rho
vars = load(sessionROIDataFilePath, 'lapCorr_Halves');
stableMask = vars.lapCorr_Halves.rho >= stableThreshold;
%% Smoothing
w = gausswin(15); w = w / sum(w);
smoothedActivity = lapActivity;
for iCell = 1:numROIs
    for iLap = 1:numLaps
        trace = squeeze(lapActivity(iCell, iLap, :));
        if all(isnan(trace)), continue; end
        nanMask = isnan(trace); trace(nanMask) = 0;
        smoothed = filtfilt(w, 1, trace); smoothed(nanMask) = NaN;
        smoothedActivity(iCell, iLap, :) = smoothed;
    end
end
%% Split Baseline Laps
baseIdx = find(contains(lower(conds), 'baseline') | contains(lower(conds), 'norm'), 1);
if isempty(baseIdx), baseIdx = 1; end
baseLaps = response.trialIndicesByCondition.(conds{baseIdx});
oddBaseline = baseLaps(1:2:end);
evenBaseline = baseLaps(2:2:end);
% Filter for stable activity
stableSmoothed = smoothedActivity(stableMask, :, :);
numStable = size(stableSmoothed, 1);
meanOdd  = squeeze(mean(stableSmoothed(:, oddBaseline, :), 2, 'omitnan'));
meanEven = squeeze(mean(stableSmoothed(:, evenBaseline, :), 2, 'omitnan'));
%% Reference Normalization
minOdd = min(meanOdd, [], 2);
maxOdd = max(meanOdd, [], 2);
rangeOdd = maxOdd - minOdd;
rangeOdd(rangeOdd == 0) = 1; 
normOdd  = (meanOdd - minOdd) ./ rangeOdd;
normEven = (meanEven - minOdd) ./ rangeOdd;
%% Align to Peak
numBins = numPosBins; 
alignedFull = nan(numStable, numBins * 2 + 1);
centerIdx = numBins + 1;
allPeakPos = zeros(numStable, 1); 

% Exclusion boundaries: 30cm to 170cm
startBin = 10;
endBin = 190;

for i = 1:numStable
    % Find peak ONLY within 30-170cm to avoid initial onsets and offsets
    [~, relativePeak] = max(normOdd(i, startBin:endBin));
    peakPos = relativePeak + (startBin - 1); 
    
    allPeakPos(i) = peakPos;
    
    evenProfile = normEven(i, :);
    
    targetStart = centerIdx - (peakPos - 1);
    targetEnd = targetStart + numBins - 1;
    
    alignedFull(i, targetStart:targetEnd) = evenProfile;
end
%% Plotting
if plotFlag
    fig = figure('Color', 'w', 'Position', [100 100 550 800]);
    fullShiftBins = (-numBins : numBins); 
    [~, snakeSortIdx] = sort(allPeakPos, 'descend'); 
    
    % Heatmap
    subplot('Position', [0.15 0.45 0.7 0.45]);
    imagesc(fullShiftBins, 1:numStable, alignedFull(snakeSortIdx, :));
    colormap(flipud(gray));
    clim([0 1]); 
    
    hold on;
    xline(-80, 'r', 'LineWidth', 1.5);
    xline(80, 'b', 'LineWidth', 1.5);
    xline(0, 'y', 'LineWidth', 2);
    
    xlim([-100, 100]); ylabel('Sorted Stable ROIs');
    title(sprintf('Aligned Population (n=%d)', numStable));
    axis off; 
    
    % Mean Plot
    subplot('Position', [0.15 0.12 0.7 0.25]);
    popMean = mean(alignedFull, 1, 'omitnan');
    numSamples = sum(~isnan(alignedFull), 1);
    popSEM = std(alignedFull, 0, 1, 'omitnan') ./ sqrt(numSamples);
    
    fill([fullShiftBins, fliplr(fullShiftBins)], [popMean+popSEM, fliplr(popMean-popSEM)], ...
         [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;
    plot(fullShiftBins, popMean, 'k', 'LineWidth', 2);
    
    xline(-80, 'r', 'LineWidth', 1.5);
    xline(80, 'b', 'LineWidth', 1.5);
    xline(0, 'y', 'LineWidth', 1.5);
    
    xlabel('Distance from peak (cm)'); 
    ylabel('Mean \DeltaF/F');
    xlim([-100, 100]); 
    xticks([-80, -40, 0, 40, 80]);
    ylim()
    %ylim([min(popMean)*0.9, max(popMean)*1.2]);
    
    set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 10);
    
    figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
    if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
    filename = fullfile(figSaveDir, sprintf('%s_%s_%s_%s_PeakAligned.png', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, signalToUse, VRStimName));
    exportgraphics(fig, filename, 'Resolution', 300);
end
response.peakAlignedBaselineTuningCurve = alignedFull;
save(responsefilePath, '-struct', 'response', '-append');
end