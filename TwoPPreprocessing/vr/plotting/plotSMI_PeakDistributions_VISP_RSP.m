function [alignedFull, fig] = plotSMI_PeakDistributions_VISP_RSP(RegionData, regionName)
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
    fill(ax2, [fullShiftBins, fliplr(fullShiftBins)], ...
         [popMean+popSEM, fliplr(popMean-popSEM)], ...
         [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    plot(ax2, fullShiftBins, popMean, 'Color', lineColor, 'LineWidth', 2);
    
    % landmark lines — all 4 distances from peak
    xline(ax2, 0,    '-y',  'LineWidth', 1.5);   % peak
    xline(ax2, -40,  '--r', 'LineWidth', 1.2);   % adjacent landmark
    xline(ax2,  40,  '--r', 'LineWidth', 1.2);
    xline(ax2, -80,  '--b', 'LineWidth', 1.2);   % same type landmark
    xline(ax2,  80,  '--b', 'LineWidth', 1.2);
    xline(ax2, -120, '--r', 'LineWidth', 1.2);   % 3rd landmark
    xline(ax2,  120, '--r', 'LineWidth', 1.2);
    xline(ax2, -160, '--b', 'LineWidth', 1.2);   % 4th landmark
    xline(ax2,  160, '--b', 'LineWidth', 1.2);
    
    xlabel(ax2, 'Distance from peak (cm)', 'FontName', 'Arial', 'FontSize', 11); 
    ylabel(ax2, 'Activity (norm.)', 'FontName', 'Arial', 'FontSize', 11);  
    xlim(ax2, [-160, 160]); 
    xticks(ax2, [-160, -120, -80, -40, 0, 40, 80, 120, 160]);
    ylim(ax2, [min(popMean - popSEM)*0.95, max(popMean + popSEM)*1.15]);
    
    hold(ax2, 'off');
    drawnow; 
    
    defaultAxesProperties(ax2, false);
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
    xline(ax1, 0,    '-y',  'LineWidth', 1.8);   % peak
    xline(ax1, -40,  '--r', 'LineWidth', 1.2);
    xline(ax1,  40,  '--r', 'LineWidth', 1.2);
    xline(ax1, -80,  '--b', 'LineWidth', 1.2);
    xline(ax1,  80,  '--b', 'LineWidth', 1.2);
    xline(ax1, -120, '--r', 'LineWidth', 1.2);
    xline(ax1,  120, '--r', 'LineWidth', 1.2);
    xline(ax1, -160, '--b', 'LineWidth', 1.2);
    xline(ax1,  160, '--b', 'LineWidth', 1.2);
    
    xlim(ax1, [-160, 160]); 
    title(ax1, sprintf('%s Aligned Population', regionName), ...
        'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'normal');
    
    % scale bar and cell count annotation
    plot(ax1, [-160, -120], [-numStable*0.04, -numStable*0.04], ...
        '-k', 'LineWidth', 2.5, 'Clipping', 'off');
    text(ax1, -140, -numStable*0.09, '40 cm', ...
        'FontName', 'Arial', 'FontSize', 10, ...
        'HorizontalAlignment', 'center', 'Clipping', 'off');
    text(ax1, -185, numStable/2, sprintf('%d %s', numStable, lower(cellLabelText)), ...
        'Rotation', 90, 'FontName', 'Arial', 'FontSize', 12, 'FontAngle', 'italic', ...
        'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Clipping', 'off');
    hold(ax1, 'off');
    
    % --- Step 3: Colorbar ---
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
    
    % listener cleanup
    h_listeners = findobj(fig, '-property', 'MarkedClean');
    for idx = 1:length(h_listeners)
        delete(findobj(h_listeners(idx), '-class', 'event.listener'));
    end
    
    % save
    baseFileName = sprintf('%s_PeakAlignedTuning.png', regionName);
    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4\tuningCurved_PeakAlignedTuning_RSP_VISp_160';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    
    fullSavePath = fullfile(outputDir, baseFileName);
    saveFigureFormats(fig, fullSavePath);
end

% function plotSMI_PeakDistributions_VISP_RSP(RSPData, VISpData)
%    
%     fprintf('Processing RSP Peak Positions from saved data...\n');
%     
%     rspPeaks_LowSMI  = [];
%     rspPeaks_HighSMI = [];
%     allRSP_SMI = [];
%     
%     rspColor  = 'k'; 
%     vispColor = [0.6 0.6 0.6]; 
%     
%     % First pass: gather SMI to find the median
%     for s = 1:length(RSPData)
%         sess = RSPData(s);
%         if isfield(sess, 'SMI') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
%             vals = sess.SMI.SMI(sess.FilteredROIs);
%             allRSP_SMI = [allRSP_SMI; vals(~isnan(vals))];
%         end
%     end
%     medianRSP_SMI = median(allRSP_SMI);
%     
%     % Second pass: Extract peaks directly from saved structure fields
%     for s = 1:length(RSPData)
%         sess = RSPData(s);
%         if ~isfield(sess, 'SMI') || ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
%             continue;
%         end
%         
%         smiValues = sess.SMI.SMI;
%         
%         % Direct loading from saved structure
%         if isfield(sess.SMI, 'RpBin')
%             peakBins = sess.SMI.RpBin;
%         elseif isfield(sess.SMI, 'GlobalPeakBin')
%             peakBins = sess.SMI.GlobalPeakBin;
%         else
%             error('Could not find peak location variables (RpBin/GlobalPeakBin) in RSPData(%d).SMI', s);
%         end
%         
%         for i = 1:length(sess.FilteredROIs)
%             roiIdx = sess.FilteredROIs(i);
%             smiVal = smiValues(roiIdx);
%             prefBin = peakBins(roiIdx);
%             
%             if isnan(smiVal) || isnan(prefBin); continue; end
%             
%             % Inclusive boundary tracking (<= 0 and >= median)
%             if smiVal <= 0
%                 rspPeaks_LowSMI = [rspPeaks_LowSMI; prefBin];
%             elseif smiVal >= medianRSP_SMI
%                 rspPeaks_HighSMI = [rspPeaks_HighSMI; prefBin];
%             end
%         end
%     end
%     
%     fprintf('Processing VISp Peak Positions from saved data...\n');
%     
%     vispPeaks_LowSMI  = [];
%     vispPeaks_HighSMI = [];
%     allVISp_SMI = [];
%     
%     % First pass: gather SMI to find the median
%     for s = 1:length(VISpData)
%         sess = VISpData(s);
%         if isfield(sess, 'SMI') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
%             vals = sess.SMI.SMI(sess.FilteredROIs);
%             allVISp_SMI = [allVISp_SMI; vals(~isnan(vals))];
%         end
%     end
%     medianVISp_SMI = median(allVISp_SMI);
%     
%     % Second pass: Extract peaks directly from saved structure fields
%     for s = 1:length(VISpData)
%         sess = VISpData(s);
%         if ~isfield(sess, 'SMI') || ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
%             continue;
%         end
%         
%         smiValues = sess.SMI.SMI;
%         
%         % Direct loading from saved structure
%         if isfield(sess.SMI, 'RpBin')
%             peakBins = sess.SMI.RpBin;
%         elseif isfield(sess.SMI, 'GlobalPeakBin')
%             peakBins = sess.SMI.GlobalPeakBin;
%         else
%             error('Could not find peak location variables (RpBin/GlobalPeakBin) in VISpData(%d).SMI', s);
%         end
%         
%         for i = 1:length(sess.FilteredROIs)
%             roiIdx = sess.FilteredROIs(i);
%             smiVal = smiValues(roiIdx);
%             prefBin = peakBins(roiIdx);
%             
%             if isnan(smiVal) || isnan(prefBin); continue; end
%             
%             % Inclusive boundary tracking (<= 0 and >= median)
%             if smiVal <= 0
%                 vispPeaks_LowSMI = [vispPeaks_LowSMI; prefBin];
%             elseif smiVal >= medianVISp_SMI
%                 vispPeaks_HighSMI = [vispPeaks_HighSMI; prefBin];
%             end
%         end
%     end
%     
%     %% --- Plotting Section ---
%     figHandle = figure('Name', 'Peak Position Count Distribution', ...
%                        'Color', [1 1 1], 'Position', [150, 150, 800, 500]);
%                    
%     binEdges = 1:1:200; 
%     landmarkCentres = [40, 80, 120, 160];
%     
%     % Subplot 1: RSP low SMI (SMI <= 0) - Black Color
%     subplot(2, 2, 1); hold on;
%     histogram(rspPeaks_LowSMI, 'BinEdges', binEdges, 'Normalization', 'count', ...
%               'FaceColor', rspColor, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
%     xlabel('Position (cm)');
%     ylabel('Count');
%     title(sprintf('RSP boutons: SMI \\leq 0\n(n = %d)', length(rspPeaks_LowSMI)));
%     xlim([0 200]);
%     xticks(landmarkCentres);
%     box off; axis square;
%     defaultAxesProperties(gca)
%     offsetAxes(gca)
%     
%     % Subplot 2: RSP high SMI (SMI >= median) - Black Color
%     subplot(2, 2, 2); hold on;
%     histogram(rspPeaks_HighSMI, 'BinEdges', binEdges, 'Normalization', 'count', ...
%               'FaceColor', rspColor, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
%     xlabel('Position (cm)');
%     ylabel('Count');
%     title(sprintf('RSP boutons: SMI \\geq median\n(n = %d)', length(rspPeaks_HighSMI)));
%     xlim([0 200]);
%     xticks(landmarkCentres);
%     box off; axis square;
%     defaultAxesProperties(gca)
%     offsetAxes(gca)
%     
%     % Subplot 3: VISp low SMI (SMI <= 0) - Gray Color
%     subplot(2, 2, 3); hold on;
%     histogram(vispPeaks_LowSMI, 'BinEdges', binEdges, 'Normalization', 'count', ...
%               'FaceColor', vispColor, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
%     xlabel('Position (cm)');
%     ylabel('Count');
%     title(sprintf('VISp somas: SMI \\leq 0\n(n = %d)', length(vispPeaks_LowSMI)));
%     xlim([0 200]);
%     xticks(landmarkCentres);
%     box off; axis square;
%     defaultAxesProperties(gca)
%     offsetAxes(gca)
%     
%     % Subplot 4: VISp high SMI (SMI >= median) - Gray Color
%     subplot(2, 2, 4); hold on;
%     histogram(vispPeaks_HighSMI, 'BinEdges', binEdges, 'Normalization', 'count', ...
%               'FaceColor', vispColor, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
%     xlabel('Position (cm)');
%     ylabel('Count');
%     title(sprintf('VISp somas: SMI \\geq median\n(n = %d)', length(vispPeaks_HighSMI)));
%     xlim([0 200]);
%     xticks(landmarkCentres);
%     box off; axis square;
%     defaultAxesProperties(gca)
%     offsetAxes(gca)
%     
%     % --- SAFE, VERSION-INDEPENDENT LISTENER CLEANUP ---
%     % Find all background listener handles assigned to this figure window and delete them.
%     % This stops offsetAxes background execution dead in its tracks right before saving.
%     figListeners = event.listener.empty;
%     h = findobj(figHandle, '-property', 'MarkedClean');
%     for idx = 1:length(h)
%         delete(findobj(h(idx), '-class', 'event.listener'));
%     end
%     % --------------------------------------------------
%     
% %     % Save execution block
% %     outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4\smiPeakDistribution_RSP_VISP';
% %     if ~exist(outputDir, 'dir')
% %         mkdir(outputDir);
% %     end
%     
%     baseFileName = 'smi_peak_position_distributions';
%     fullSavePath = fullfile(outputDir, baseFileName);
%     
%     saveFigureFormats(figHandle, fullSavePath);
% end