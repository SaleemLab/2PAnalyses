function saleemMetric = getSaleemPositionTuning(sessionFileInfo, response, signalToUse, plotFlag)
% Identifies position-tuned ROIs based on Saleem et al. (2018).
% Uses Trial Correlation and Shuffle-based Significance for selection.

%% 1. Setup
if nargin < 3 || isempty(signalToUse); signalToUse = 'dFFNeuropilCorrected'; end 
if nargin < 4 || isempty(plotFlag); plotFlag = true; end

% Target percentile for significance (Saleem typically uses p < 0.05 or 0.01)
targetPercentile = 95; 
% Reliability threshold (Saleem often uses > 0.1 for calcium imaging)
relThreshold = 0.1;

figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'); mkdir(figSaveDir); end
filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' signalToUse '_SaleemTuningSanityCheck.png']);

%% 2. Data Extraction
realActivity = response.lapPositionActivity.(signalToUse); % Z-scored
shuffMatrix = response.lapPositionActivity_ShuffleMatrix.(signalToUse);
[numROIs, numLaps, numBins] = size(realActivity);

%% 3. Calculate Saleem Metrics
meanTuning = squeeze(mean(realActivity, 2, 'omitnan')); 
trialCorr = nan(numROIs, 1);
isSignificant = zeros(numROIs, 1);

for i = 1:numROIs
    % A. Calculate Trial-to-Trial Reliability (Median Pairwise Pearson Correlation)
    roiLaps = squeeze(realActivity(i, :, :)); 
    corrMat = corr(roiLaps', 'Rows', 'pairwise');
    trialCorr(i) = median(corrMat(triu(true(size(corrMat)), 1)), 'omitnan');
    
    % B. Calculate Significance via Shuffle Distribution
    shuffMeanDist = squeeze(mean(shuffMatrix(i, :, :), 2, 'omitnan'));
    threshold_sig = prctile(shuffMeanDist, targetPercentile);
    
    % ROI is significant if its peak mean activity exceeds the shuffle threshold
    if max(meanTuning(i, :)) > threshold_sig
        isSignificant(i) = 1;
    end
end

% Combined Inclusion: Stable (Reliability) AND Non-Random (Significant)
isPositionTuned = (trialCorr > relThreshold) & (isSignificant == 1);

%% 4. Save to sessionROIData
saleemMetric.trialCorr = trialCorr;
saleemMetric.isPositionTuned = isPositionTuned;
if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, "saleemMetric", '-append');
end

%% 5. Plotting (Sanity Check)
if plotFlag && ~isempty(realActivity)
    % Prepare Odd/Even laps for cross-validation
    oddIdx = 1:2:numLaps; evenIdx = 2:2:numLaps;
    meanOdd = squeeze(mean(realActivity(:, oddIdx, :), 2, 'omitnan'));
    meanEven = squeeze(mean(realActivity(:, evenIdx, :), 2, 'omitnan'));
    
    % Normalize for heatmap display
    normOdd = nan(numROIs, numBins); normEven = nan(numROIs, numBins);
    for r = 1:numROIs
        ref = prctile(meanOdd(r,:), 98) + eps;
        normOdd(r,:) = meanOdd(r,:) ./ ref;
        normEven(r,:) = meanEven(r,:) ./ ref;
    end
    
    sigIdx = find(isPositionTuned);
    [~, peakAll] = max(meanOdd, [], 2); [~, sortAll] = sort(peakAll);
    
    fig1 = figure('Name', 'Saleem Position Tuning', 'Position', [100 100 1200 800], 'Color', 'w');
    t = tiledlayout(2, 2, 'TileSpacing', 'compact');
    
    % All ROIs: Shows global structure vs noise
    ax(1) = nexttile; imagesc(normOdd(sortAll,:)); title(sprintf('ALL: Odd Laps (n=%d)', numROIs));
    ax(2) = nexttile; imagesc(normEven(sortAll,:)); title('ALL: Even Laps (Sorted by Odd)');
    
    % Significant ROIs: Shows how many pass and if they align
    if ~isempty(sigIdx)
        [~, peakSig] = max(meanOdd(sigIdx,:), [], 2); [~, sortSig] = sort(peakSig);
        ax(3) = nexttile; imagesc(normOdd(sigIdx(sortSig),:)); 
        title(sprintf('PASSED: r > %0.1f & p < 0.05 (n=%d)', relThreshold, length(sigIdx)));
        ax(4) = nexttile; imagesc(normEven(sigIdx(sortSig),:)); 
        title('PASSED: Even Laps (Sorted by Odd)');
    end
    
    % Global Formatting
    for i = 1:length(ax)
        if ~isgraphics(ax(i)), continue; end
        colormap(ax(i), flipud(gray)); set(ax(i), 'CLim', [0 1.1], 'YDir', 'normal');
        hold(ax(i), 'on');
        for x = [40 80 120 160], xline(ax(i), x, 'r:', 'Alpha', 0.4); end
    end
    
    title(t, sprintf('%s | %s | %0.1f%% Units Spatially Stable', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, (length(sigIdx)/numROIs)*100));
    exportgraphics(fig1, filename, 'Resolution', 300);
end
end