function [alignedFull, response] = alignToPeakTuningCurvePosition(sessionFileInfo, VRStimName, signalToUse, plotFlag)
% findlandmarktunedrois - identifies rois that peak at specific landmark positions.
% uses odd laps to find peaks and even laps to validate consistency.
if nargin < 3, signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 4, plotFlag = 1; end
%% load data
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('specified vrstimname not found.'); end

responsefilePath = sessionFileInfo.stimFiles(stimIdx).Response;
sessionROIDataFilePath = sessionFileInfo.otherSessFilePaths.sessionROIData;

% loading lap activity and the trial conditions
response = load(responsefilePath, 'lapPositionActivity', 'trialIndicesByCondition');
lapActivity = response.lapPositionActivity.(signalToUse);
conds = fieldnames(response.trialIndicesByCondition);
[numROIs, numLaps, numPosBins] = size(lapActivity);

% pick stable rois 
vars = load(sessionROIDataFilePath, 'lapCorr_Halves');
stableIdx = vars.lapCorr_Halves.stableIdx;
%% 
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
filename = fullfile(figSaveDir, sprintf('%s_%s_%s_%s_AlignedPeakResponses.png', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, signalToUse, VRStimName));


%% smoothing (spatial filtering)
w = gausswin(15); w = w / sum(w);
smoothedActivity = lapActivity;
for iCell = 1:numROIs
    for iLap = 1:numLaps
        trace = squeeze(lapActivity(iCell, iLap, :));
        if all(isnan(trace)), continue; end
        nanMask = isnan(trace);
        trace(nanMask) = 0;
        smoothed = filtfilt(w, 1, trace);
        smoothed(nanMask) = NaN;
        smoothedActivity(iCell, iLap, :) = smoothed;
    end
end

%% pick out baseline laps and split into odd and even 
baseIdx = find(contains(lower(conds), 'baseline') | contains(lower(conds), 'norm'), 1);
if isempty(baseIdx), baseIdx = 1; end
baseLaps = response.trialIndicesByCondition.(conds{baseIdx});

% cross-validate and split 
oddBaseline = baseLaps(1:2:end);
evenBaseline = baseLaps(2:2:end);
meanOdd  = squeeze(mean(smoothedActivity(:, oddBaseline, :), 2, 'omitnan'));
meanEven = squeeze(mean(smoothedActivity(:, evenBaseline, :), 2, 'omitnan'));

normOdd = normalize(meanOdd, 2, 'range');
normEven = normalize(meanEven, 2, 'range');

%% Look 90cm left/right to see the 80cm twin 

% shiftRange = 90;              
% shiftBins = -shiftRange:shiftRange;

%%align using full track with nan-padding
roisToInclude = stableIdx; 
numBins = numPosBins; % 200 bins for 200cm track

% add nan twice the track size so any peak can be centered at 0 
alignedFull = nan(length(roisToInclude), numBins * 2 + 1);
centerIdx = numBins + 1;
allPeakPos = zeros(length(roisToInclude), 1); % store positions for plotting 

for i = 1:length(roisToInclude)
    iROI = roisToInclude(i);
    
    % Find the anchor peak position on odd laps 
    [~, peakPos] = max(normOdd(iROI, :));
    allPeakPos(i) = peakPos; % Save for the snake-sort later
    
    % Get the Even Profile (the data to be plotted)
    evenProfile = normEven(iROI, :);
    
    % Calculate where the start of the track lands so peakPos hits centerIdx
    targetStart = centerIdx - (peakPos - 1);
    targetEnd = targetStart + numBins - 1;
    
    % Drop the profile into the NaN canvas
    alignedFull(i, targetStart:targetEnd) = evenProfile;
end

%% heatmap + mean 
if plotFlag
    fig = figure('Color', 'w', 'Position', [100 100 550 800]);
    nROIs = size(alignedFull, 1);
    fullShiftBins = (-numBins : numBins); 
    [~, snakeSortIdx] = sort(allPeakPos, 'descend'); 
    
    % heatmap 
    subplot('Position', [0.15 0.45 0.7 0.45]);
    imagesc(fullShiftBins, 1:nROIs, alignedFull(snakeSortIdx, :));
    colormap(flipud(gray));
    clim([0.25, 0.75]); 
    
    hold on;
    xline(-80, 'r', 'LineWidth', 1.5);
    xline(80, 'b', 'LineWidth', 1.5);
    xline(0, 'y', 'LineWidth', 1.5);
    
    xlim([-100, 100]); 
    ylabel('Stable ROIs');
    title(sprintf('Aligned Population (n=%d)', nROIs));
    
    % remove ticks and clean axis 
    axis off;
    % text(-110, nROIs/2, 'Stable ROIs', 'Rotation', 90, 'HorizontalAlignment', 'center', 'FontSize', 12);
    % ylabel('Stable ROIs');

    % Mean across spatially aliged bins 
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
    ylabel('Mean Response');
    xlim([-100, 100]); 
    xticks([-80, -40, 0, 40, 80]);
    ylim([min(popMean)*0.9, max(popMean)*1.2]);
    
    % CLEAN UP BOTTOM AXIS
    set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 10);
    % defaultAxesProperties(gca, true)
    
    exportgraphics(fig, filename, 'Resolution', 300);
end

response.peakAlignedBaselineTuningCurve = alignedFull;

disp(['Saving peak aligned-baseline tuning curves to ', sessionFileInfo.stimFiles(stimIdx).Response]);
save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append');
end

