targetSessionIdx = 4; 
thisSession = RSPData(targetSessionIdx); 

condNames = fieldnames(thisSession.ConditionData);
targetCondition = condNames{1}; 

meanOdd = squeeze(mean(thisSession.ConditionData.(targetCondition).LapActivity(:, 1:2:end, :), 2, 'omitnan'));
minOdd = min(meanOdd, [], 2);
maxOdd = max(meanOdd, [], 2);
rangeOdd = maxOdd - minOdd;
rangeOdd(rangeOdd == 0) = 1; 
normOdd_real = (meanOdd - minOdd) ./ rangeOdd; 

[numROIs, numBins] = size(normOdd_real);
x = 1:numBins;

excludeEdgePeakFlags = thisSession.SMI.ExcludeEdgePeakCells; 
allFilteredKeptROIs = thisSession.FilteredROIs;

edgeExcludedROIsIdx = find(excludeEdgePeakFlags);
keptROIsIdx = allFilteredKeptROIs;

if ~isempty(edgeExcludedROIsIdx)
    numToPlot = min(6, length(edgeExcludedROIsIdx)); 
    figure('Name', sprintf('Audit: Excluded Edge-Peak ROIs (%s)', thisSession.Session), ...
           'Color', [1 1 1], 'Position', [50, 100, 1500, 700]);
       
    for i = 1:numToPlot
        roiID = edgeExcludedROIsIdx(i);
        subplot(2, 3, i); hold on;
        
        plot(x, normOdd_real(roiID, :), 'Color', [0.824, 0.016, 0.176], 'LineWidth', 2); 
        
        patch([0 30 30 0], [0 0 1.2 1.2], [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        patch([numBins-30 numBins numBins numBins-30], [0 0 1.2 1.2], [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        
        gPeak = thisSession.SMI.GlobalPeakBin(roiID);
        if gPeak <= numBins && gPeak > 0
            plot(gPeak, normOdd_real(roiID, gPeak), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
        end
        
        title(sprintf('ROI %d\nGlobal Peak Bin: %d', roiID, gPeak));
        xlabel('Position Bins'); ylabel('Norm Activity'); ylim([0 1.2]); xlim([1 numBins]); grid on;
    end
    sgtitle('Flagged Background/Edge Responsive ROIs (should have peaks in shaded zones)');
else
    disp('No ROIs were excluded by the Edge-Peak filter in this session.');
end

if ~isempty(keptROIsIdx)
    numToPlot = min(6, length(keptROIsIdx)); 
    figure('Name', sprintf('Audit: Kept ROIs (%s)', thisSession.Session), ...
           'Color', [1 1 1], 'Position', [100, 150, 1500, 700]);
       
    for i = 1:numToPlot
        roiID = keptROIsIdx(i);
        subplot(2, 3, i); hold on;
        
        plot(x, normOdd_real(roiID, :), 'Color', [0.000, 0.400, 1.000], 'LineWidth', 2); 
        
        patch([0 30 30 0], [0 0 1.2 1.2], [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        patch([numBins-30 numBins numBins numBins-30], [0 0 1.2 1.2], [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        
        gPeak = thisSession.SMI.GlobalPeakBin(roiID);
        if gPeak <= numBins && gPeak > 0
            plot(gPeak, normOdd_real(roiID, gPeak), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
        end
        
        title(sprintf('ROI %d\nGlobal Peak Bin: %d', roiID, gPeak));
        xlabel('Position Bins'); ylabel('Norm Activity'); ylim([0 1.2]); xlim([1 numBins]); grid on;
    end
    sgtitle('Clean Kept ROIs (Should have peaks out of shaded zones)');
else
    disp('No ROIs survived the complete filter set for this session.');
end