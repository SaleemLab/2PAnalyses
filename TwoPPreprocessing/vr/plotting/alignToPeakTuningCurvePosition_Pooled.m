function [alignedFull, fig] = alignToPeakTuningCurvePosition_Pooled(RegionData, regionName)
    if nargin < 2, regionName = 'Region'; end
    fprintf('Processing %s population alignment...\n', regionName);
    
    if strcmpi(regionName, 'RSP')
        lineColor = 'k';
        cellLabelText = 'boutons';
    else
        lineColor = [0.6 0.6 0.6]; 
        cellLabelText = 'somas';
    end
    allAlignedProfiles = [];
    allPeakPositions   = [];
    numBins = [];
    
    for s = 1:length(RegionData)
        sess = RegionData(s);
        
        if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
           ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
            continue;
        end
        
        lapActivity = sess.ConditionData.Baseline.LapActivity;
        [numROIsTotal, numLaps, numPosBins] = size(lapActivity);
        numBins = numPosBins; 
        
        % --- SPATIAL SMOOTHING ONLY ---
        w_space = gausswin(15); w_space = w_space / sum(w_space);
        smoothedActivity = lapActivity;
        
        for iCell = 1:numROIsTotal
            for iLap = 1:numLaps
                trace = squeeze(lapActivity(iCell, iLap, :));
                if all(isnan(trace)), continue; end
                nanMask = isnan(trace); trace(nanMask) = 0;
                smoothed = filtfilt(w_space, 1, trace); smoothed(nanMask) = NaN;
                smoothedActivity(iCell, iLap, :) = smoothed;
            end
        end
        roisToAnalyze = sess.FilteredROIs;
        roiActivity   = smoothedActivity(roisToAnalyze, :, :);
        numROIs       = length(roisToAnalyze);
        
        if numROIs == 0, continue; end
        
        oddLaps  = 1:2:numLaps;
        evenLaps = 2:2:numLaps;
        
        meanOdd  = squeeze(mean(roiActivity(:, oddLaps, :), 2, 'omitnan'));
        meanEven = squeeze(mean(roiActivity(:, evenLaps, :), 2, 'omitnan'));
        
        if numROIs == 1
            meanOdd  = meanOdd';
            meanEven = meanEven';
        end
        
        minOdd = min(meanOdd, [], 2);
        maxOdd = max(meanOdd, [], 2);
        rangeOdd = maxOdd - minOdd;
        rangeOdd(rangeOdd == 0) = 1; 
        
        normOdd  = (meanOdd - minOdd) ./ rangeOdd;
        normEven = (meanEven - minOdd) ./ rangeOdd;
        
        centerIdx = numBins + 1;
        sessionAligned = nan(numROIs, numBins * 2 + 1);
        sessionPeaks   = zeros(numROIs, 1);
        
        startBin = 30;
        endBin = 170;
        
        for i = 1:numROIs
            [~, relativePeak] = max(normOdd(i, startBin:endBin));
            peakPos = relativePeak + (startBin - 1); 
            sessionPeaks(i) = peakPos;
            
            evenProfile = normEven(i, :);
            targetStart = centerIdx - (peakPos - 1);
            targetEnd   = targetStart + numBins - 1;
            
            sessionAligned(i, targetStart:targetEnd) = evenProfile;
        end
        
        allAlignedProfiles = [allAlignedProfiles; sessionAligned];
        allPeakPositions   = [allPeakPositions; sessionPeaks];
    end
    
    if isempty(allAlignedProfiles)
        error('No valid cells found tracking across the filtered dataset entries.');
    end
    
    alignedFull = allAlignedProfiles;
    numStable = size(alignedFull, 1);
    
    fprintf('Successfully grouped and evaluated %d stable %s for %s.\n', ...
        numStable, lower(cellLabelText), regionName);
    
    %% --- Layout Geometry Setup ---
    fig = figure('Color', 'w', 'Position', [100 100 550 700]);
    fullShiftBins = (-numBins : numBins); 
    [~, snakeSortIdx] = sort(allPeakPositions, 'descend'); 
    
    heatmapPos = [0.15 0.42 0.68 0.50];
    tracePos   = [0.15 0.12 0.68 0.22];
    
    % --- Step 1: Render Bottom Mean Trend Plot ---
    ax2 = axes('Position', tracePos);
    popMean = mean(alignedFull, 1, 'omitnan');
    numSamples = sum(~isnan(alignedFull), 1);
    popSEM = std(alignedFull, 0, 1, 'omitnan') ./ sqrt(numSamples);
    
    hold(ax2, 'on');
    fill(ax2, [fullShiftBins, fliplr(fullShiftBins)], [popMean+popSEM, fliplr(popMean-popSEM)], ...
         [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    plot(ax2, fullShiftBins, popMean, 'Color', lineColor, 'LineWidth', 2);
    
    xline(ax2, -80, '--r', 'LineWidth', 1.2);
    xline(ax2, 80, '--b', 'LineWidth', 1.2);
    xline(ax2, 0, '-y', 'LineWidth', 1.5);
    
    xlabel(ax2, 'Distance from peak (cm)', 'FontName', 'Arial', 'FontSize', 11); 
    ylabel(ax2, 'Number of ROIs', 'FontName', 'Arial', 'FontSize', 11);  
    xlim(ax2, [-100, 100]); 
    xticks(ax2, [-80, -40, 0, 40, 80]);
    ylim(ax2, [min(popMean)*0.95, max(popMean)*1.15]);
    
    hold(ax2, 'off');
    drawnow; 
    
    % MATCHING REFERENCE: Pass false so the layout doesn't override width dimensions
    defaultAxesProperties(ax2, false);
    
    % Force positioning coordinates to match layout boundaries
    set(ax2, 'Position', tracePos, 'Box', 'off', 'TickDir', 'out');
    
    % --- Step 2: Render Top Heatmap Axis ---
    ax1 = axes('Position', heatmapPos);
    displayMatrix = alignedFull(snakeSortIdx, :);
    imagesc(ax1, fullShiftBins, 1:numStable, displayMatrix);
    colormap(ax1, flipud(gray));
    
    set(ax1, 'CLim', [0.25 0.75], 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out');
    set(ax1, 'YTick', [], 'YColor', 'none');
    set(ax1, 'XTick', [], 'XColor', 'none'); 
    set(ax1, 'Position', heatmapPos);
    
    hold(ax1, 'on');
    xline(ax1, -80, '--r', 'LineWidth', 1.2);
    xline(ax1, 80, '--b', 'LineWidth', 1.2);
    xline(ax1, 0, '-y', 'LineWidth', 1.8);
    xlim(ax1, [-100, 100]); 
    title(ax1, sprintf('%s Aligned Population', regionName), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'normal');
    
    % Scale indicators and custom dynamic annotations
    plot(ax1, [-100, -60], [-numStable*0.04, -numStable*0.04], '-k', 'LineWidth', 2.5, 'Clipping', 'off');
    text(ax1, -80, -numStable*0.09, '40 cm', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Clipping', 'off');
    text(ax1, -112, numStable/2, sprintf('%d %s', numStable, lower(cellLabelText)), ...
        'Rotation', 90, 'FontName', 'Arial', 'FontSize', 12, 'FontAngle', 'italic', ...
        'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Clipping', 'off');
    hold(ax1, 'off');
    
    % --- Step 3: Clean and Rebuild Isolated Side Colorbar Layout ---
    delete(findobj(fig, 'Type', 'colorbar'));
    
    cb = colorbar(ax1); 
    cb.Position = [0.86, 0.42, 0.03, 0.50];
    cb.Ticks = [0.25 0.50 0.75]; 
    cb.TickLabels = {'0.25', '0.50', '0.75'}; 
    cb.TickDirection = 'out'; 
    cb.Box = 'off';
    cb.FontName = 'Arial'; 
    cb.FontSize = 10;
    cb.Label.String = 'Activity (norm.)'; 
    cb.Label.FontName = 'Arial'; 
    cb.Label.FontSize = 12; 
    cb.Label.Rotation = 90; 
    cb.Label.Units = 'normalized'; 
    cb.Label.Position = [4, 0.5, 0]; 
    cb.Label.HorizontalAlignment = 'center'; 
    cb.Label.VerticalAlignment = 'bottom';
    
    drawnow;
    
    h_listeners = findobj(fig, '-property', 'MarkedClean');
    for idx = 1:length(h_listeners)
        delete(findobj(h_listeners(idx), '-class', 'event.listener'));
    end
    
    baseFileName = sprintf('%s_PeakAlignedTuning.png', regionName);
    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4\tuningCurved_PeakAlignedTuning_RSP_VISp';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    
    fullSavePath = fullfile(outputDir, baseFileName);
    saveFigureFormats(fig, fullSavePath);
end