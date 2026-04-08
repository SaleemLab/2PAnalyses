function conservativeTuningMetric = getConservativeTuning_ZScoreAndSignificance(sessionFileInfo, response, signalToUse, plotFlag)
% Calculates Tuning Significance using the Wei et al. (2026) Lower Bound method.
% Optimized for Z-scored dFF data with 1 cm binning.


if nargin < 3 || isempty(signalToUse); signalToUse = 'dFFNeuropilCorrected'; end 
if nargin < 4 || isempty(plotFlag); plotFlag = true; end

% Paper uses 97.5th percentile for the shuffle test 
targetPercentile = 97.5; 

figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'); mkdir(figSaveDir); end
filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' signalToUse '_WeiMethod_1cm.png']);

%% 
realActivity = response.lapPositionActivity.(signalToUse); 
shuffMatrix = response.lapPositionActivity_ShuffleMatrix.(signalToUse);
[numROIs, numLaps, numBins] = size(realActivity);

%% 
% Calculate mean across laps 
meanTuning = squeeze(mean(realActivity, 2, 'omitnan')); 

% Calculate SEM to establish the Lower Bound 
semTuning = squeeze(std(realActivity, 0, 2, 'omitnan') ./ sqrt(numLaps));
lowerBoundTuning = meanTuning - semTuning; 

isConservativeStable = zeros(numROIs, 1);

for i = 1:numROIs
    % Generate shuffle distribution from the mean of shuffled laps 
    shuffMeanDist = squeeze(mean(shuffMatrix(i, :, :), 2, 'omitnan'));
    threshold_sig = prctile(shuffMeanDist, targetPercentile);
    
    % Lower bound (Mean-SEM) must exceed shuffle threshold 
    sigBins = lowerBoundTuning(i, :) > threshold_sig;
    
    % Paper requires at least 15 cm 
    % Since 1 bin = 1 cm, we require 15 consecutive bins.
    consecutiveSig = diff([0, sigBins, 0]);
    starts = find(consecutiveSig == 1);
    stops = find(consecutiveSig == -1);
    
    % Field width check: (stops - starts) gives number of bins
    if any((stops - starts) >= 15) 
        isConservativeStable(i) = 1;
    end
end

%%
conservativeTuningMetric.lowerBoundTuning = lowerBoundTuning;
conservativeTuningMetric.isConservativeStable = isConservativeStable;
conservativeTuningMetric.targetPercentile = targetPercentile;

if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, "conservativeTuningMetric", '-append');
end

%% 
if plotFlag && ~isempty(realActivity)

    oddIdx = 1:2:numLaps; evenIdx = 2:2:numLaps;
    meanOdd = squeeze(mean(realActivity(:, oddIdx, :), 2, 'omitnan'));
    meanEven = squeeze(mean(realActivity(:, evenIdx, :), 2, 'omitnan'));
    
    % Normalize 0-1 based on 98th percentile to suppress noise spikes
    normOdd = nan(numROIs, numBins); normEven = nan(numROIs, numBins);
    for r = 1:numROIs
        refVal = prctile(meanOdd(r,:), 98) + eps;
        normOdd(r,:) = meanOdd(r,:) ./ refVal;
        normEven(r,:) = meanEven(r,:) ./ refVal;
    end
    
    stableIdx = find(isConservativeStable);
    [~, peakPos] = max(meanOdd, [], 2); [~, sortOrder] = sort(peakPos);
    
    fig1 = figure('Name', 'Wei et al. Method (1cm Bins)', 'Position', [100 100 1200 800], 'Color', 'w');
    t = tiledlayout(2, 2, 'TileSpacing', 'compact');
    
    % Heatmaps for all ROIs
    ax(1) = nexttile; imagesc(normOdd(sortOrder,:)); title('ALL: Odd Laps');
    ax(2) = nexttile; imagesc(normEven(sortOrder,:)); title('ALL: Even (Sorted by Odd)');
    
    % Heatmaps for Significant ROIs
    if ~isempty(stableIdx)
        [~, pSig] = max(meanOdd(stableIdx,:), [], 2); [~, sSig] = sort(pSig);
        ax(3) = nexttile; imagesc(normOdd(stableIdx(sSig),:)); 
        title(sprintf('SIG: Odd Laps (n=%d)', length(stableIdx)));
        ax(4) = nexttile; imagesc(normEven(stableIdx(sSig),:)); 
        title('SIG: Even (Sorted by Odd)');
    end
    
    % Formatting to match paper style
    for i = 1:length(ax)
        if ~isgraphics(ax(i)), continue; end
        colormap(ax(i), flipud(gray)); set(ax(i), 'CLim', [0 1.1], 'YDir', 'normal');
        hold(ax(i), 'on');
        % Mark 40, 80, 120, 160 cm positions
        for x = [40 80 120 160], xline(ax(i), x, 'r:', 'Alpha', 0.4); end
    end
    
    title(t, sprintf('Method: Lower Bound (Mean-SEM) > %0.1fth %%ile | 15cm Min Width', targetPercentile));
    exportgraphics(fig1, filename, 'Resolution', 300);
end
end