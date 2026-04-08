function siMetric = getSpatialInfo_WithShuffle(sessionFileInfo, response, signalToUse, plotFlag)
% Calculates Spatial Mutual Information (Skaggs et al. 1993) 
% and corrects it using a shuffled null distribution.

%% 1. Setup
if nargin < 3 || isempty(signalToUse); signalToUse = 'dFFNeuropilCorrected'; end 
if nargin < 4 || isempty(plotFlag); plotFlag = true; end

figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'); mkdir(figSaveDir); end

filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' signalToUse '_SpatialInfo_SanityCheck.png']);

%% 2. Data Extraction
realActivity = response.lapPositionActivity.(signalToUse); 
shuffMatrix = response.lapPositionActivity_ShuffleMatrix.(signalToUse);
[numROIs, numLaps, numBins] = size(realActivity);

% Mean tuning curve for SI calculation
meanTuning = squeeze(mean(realActivity, 2, 'omitnan')); % [ROIs x Bins]
% Replace any negative values (from z-scoring) with 0 for SI log math
meanTuning(meanTuning < 0) = 0; 
shuffMatrix(shuffMatrix < 0) = 0;

%% 3. Calculate Spatial Information (SI)
% P is occupancy probability (assuming uniform for fixed-speed/fixed-track)
P = 1 / numBins; 

% Helper function for SI (Skaggs Formula)
calcSI = @(tc) sum(P .* (tc ./ (mean(tc) + eps)) .* log2((tc ./ (mean(tc) + eps)) + eps));

rawSI = nan(numROIs, 1);
shuffSI = nan(numROIs, size(shuffMatrix, 3));

disp('Calculating Mutual Information...');
for i = 1:numROIs
    % Real SI
    rawSI(i) = calcSI(meanTuning(i, :));
    
    % Shuffled SI (The Bootstrapped Null Distribution)
    for j = 1:size(shuffMatrix, 3)
        shuffSI(i, j) = calcSI(shuffMatrix(i, :, j));
    end
end

% Calculate Z-scored SI (The "Lap-Corrected" Information)
mu_s = mean(shuffSI, 2);
sigma_s = std(shuffSI, 0, 2);
zSI = (rawSI - mu_s) ./ (sigma_s + eps);

% Significance: Real SI > 95th percentile of Shuffled SI
siThreshold = prctile(shuffSI, 95, 2);
isSISignificant = rawSI > siThreshold;

%% 4. Plotting Sanity Check
if plotFlag
    % Prep Odd/Even for display
    oddIdx = 1:2:numLaps; evenIdx = 2:2:numLaps;
    normOdd = normalize(squeeze(mean(realActivity(:, oddIdx, :), 2, 'omitnan')), 2, 'range');
    normEven = normalize(squeeze(mean(realActivity(:, evenIdx, :), 2, 'omitnan')), 2, 'range');
    
    % Sort by SI Z-Score
    [sortedZ, sortIdx] = sort(zSI, 'descend');
    sigIdx = find(isSISignificant);
    
    fig1 = figure('Name', 'Spatial Information Check', 'Position', [100 100 1000 800], 'Color', 'w');
    t = tiledlayout(2, 2, 'TileSpacing', 'compact');
    
    % Top Row: Top 50 "Most Informative" ROIs
    topN = min(50, numROIs);
    ax1 = nexttile; imagesc(normOdd(sortIdx(1:topN),:)); title('Top 50 SI: Odd Laps');
    ax2 = nexttile; imagesc(normEven(sortIdx(1:topN),:)); title('Top 50 SI: Even Laps');
    
    % Bottom Row: Significant SI ROIs (Sorted by Peak)
    if ~isempty(sigIdx)
        normOddSig = normOdd(sigIdx, :);
        [~, peaks] = max(normOddSig, [], 2); [~, s] = sort(peaks);
        ax3 = nexttile; imagesc(normOddSig(s,:)); title(sprintf('Sig SI: Odd (n=%d)', length(sigIdx)));
        ax4 = nexttile; imagesc(normEven(sigIdx(s),:)); title('Sig SI: Even (Sorted by Odd)');
    else
        nexttile; axis off; nexttile; axis off;
    end
    
    colormap(flipud(gray));
    allAx = [ax1, ax2, ax3, ax4];
    for i = 1:4
        if isgraphics(allAx(i))
            set(allAx(i), 'CLim', [0 1], 'TickDir', 'out', 'box', 'off', 'YDir', 'normal');
            hold(allAx(i), 'on');
            for xpos = [40 80 120 160], xline(allAx(i), xpos, 'r:', 'LineWidth', 1); end
        end
    end
    
    title(t, sprintf('Spatial Information | %s | %d Significant ROIs', ...
        sessionFileInfo.session_name, length(sigIdx)));
    exportgraphics(fig1, filename, 'Resolution', 300);
end

%% 5. Save
siMetric.rawSI = rawSI;
siMetric.zSI = zSI;
siMetric.isSignificant = isSISignificant;

save(sessionFileInfo.otherSessFilePaths.sessionROIData, 'siMetric', '-append');
end