% function saveSpeedPositionActivity_AllROIs(sessionFileInfo, response, options)
% % saveallroistopdf - loops through all rois and saves them into a single pdf.
% % every roi gets its own page.
% 
%     % 1. setup defaults and directory
%     if nargin < 3, options = struct(); end
%     if ~isfield(options, 'applySmoothing'), options.applySmoothing = true; end
%     if ~isfield(options, 'smoothSigma'),    options.smoothSigma = [1.1, 1.5]; end
% 
%     figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'SpeedPositionSummaries');
%     if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
% 
%     % 2. define the master pdf filename
%     pdfName = sprintf('%s_%s_AllROIs_SpeedPosition.pdf', ...
%         sessionFileInfo.animal_name, sessionFileInfo.session_name);
%     fullPDFPath = fullfile(figSaveDir, pdfName);
% 
%     % 3. delete existing pdf if it exists (so you don't keep appending to old files)
%     if exist(fullPDFPath, 'file'), delete(fullPDFPath); end
% 
%     % 4. get total number of rois from the response matrix
%     numROIs = size(response.speedPositionActivity.matrix, 3);
%     fprintf('starting pdf generation for %d rois...\n', numROIs);
% 
%     % 5. the loop
%     for targetROI = 1:numROIs
%         % call your plotting function
%   
%         plotSpeedPositionActivity_ForROI(sessionFileInfo, response, targetROI, ...
%             options.applySmoothing, options.smoothSigma);
%         % 
%         % find the figure that was just created
%         figHandle = gcf; 
%         
%         % append to pdf
%         % the first time this runs, it creates the file; after that, it appends pages.
%         exportgraphics(figHandle, fullPDFPath, 'ContentType', 'vector', 'Append', true);
%         
%         % close figure immediately to save ram
%         close(figHandle);
%         
%         % console progress update
%         if mod(targetROI, 10) == 0 || targetROI == numROIs
%             fprintf('processed %d/%d rois...\n', targetROI, numROIs);
%         end
%     end
% 
%     fprintf('success! all rois saved to: %s\n', fullPDFPath);
% end

function saveSpeedPositionActivity_AllROIs(sessionFileInfo, response, options)
    if nargin < 3, options = struct(); end
    if ~isfield(options, 'applySmoothing'), options.applySmoothing = true; end
    if ~isfield(options, 'smoothSigma'),    options.smoothSigma = [1.1, 1.5]; end

    figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'SpeedPositionSummaries');
    if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end

    pdfName = sprintf('%s_%s_AllROIs_SpeedPosition.pdf', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name);
    fullPDFPath = fullfile(figSaveDir, pdfName);

    if exist(fullPDFPath, 'file'), delete(fullPDFPath); end

    numROIs = size(response.speedPositionActivity.matrix, 3);
    fprintf('Starting master PDF generation for %d ROIs...\n', numROIs);

    for targetROI = 1:numROIs
        figHandle = plotSpeedPositionActivity_ForROI_Internal(sessionFileInfo, response, targetROI, ...
            options.applySmoothing, options.smoothSigma);
        
        exportgraphics(figHandle, fullPDFPath, 'ContentType', 'vector', 'Append', true);
        close(figHandle);
        
        if mod(targetROI, 10) == 0 || targetROI == numROIs
            fprintf('Processed %d/%d ROIs...\n', targetROI, numROIs);
        end
    end
    fprintf('Success! All ROIs saved to: %s\n', fullPDFPath);
end

function figHandle = plotSpeedPositionActivity_ForROI_Internal(sessionFileInfo, response, targetROI, applySmoothing, smoothSigma)


    tuningSurface = response.speedPositionActivity.matrix(:, :, targetROI);
    speedCenters  = response.speedPositionActivity.speedBinCenters;
    speedEdges    = response.speedPositionActivity.speedEdges;
    [numSpeedBins, numPosBins] = size(tuningSurface);
    
    if applySmoothing
        mask = ~isnan(tuningSurface);
        dataZeroed = tuningSurface;
        dataZeroed(~mask) = 0; 
        blurredData = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
        blurredMask = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
        tuningSurface = blurredData ./ blurredMask;
        tuningSurface(isnan(tuningSurface)) = 0;
    else
        tuningSurface(isnan(tuningSurface)) = 0;
    end
    
    figHandle = figure('Name', sprintf('Bouton %d', targetROI), ...
        'Position', [100 100 720 600], 'Color', 'w', 'Visible', 'off');
    
    leftMargin = 0.12;
    plotWidth  = 0.62; 
    
    ax1 = subplot('Position', [leftMargin, 0.68, plotWidth, 0.24]); 
    hold on;
    
    landmarkCentres = [40, 80, 120, 160];
    for c = landmarkCentres
        if c <= numPosBins
            xline(c, ':', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2);
        end
    end
    
    x_pos = 1:numPosBins;
    lowLine  = mean(tuningSurface(1:3, :), 1, 'omitnan');
    medLine  = mean(tuningSurface(4:7, :), 1, 'omitnan');
    highLine = mean(tuningSurface(8:10, :), 1, 'omitnan');
    
    lowThreshold  = speedEdges(4);
    highThreshold = speedEdges(8);
    
    plot(x_pos, lowLine, 'Color', [0, 0.9, 1], 'LineWidth', 2);
    plot(x_pos, medLine, 'Color', 'k', 'LineWidth', 2);
    plot(x_pos, highLine, 'Color', [1, 0, 0.9], 'LineWidth', 2);
    
    ylabel('\DeltaF/F (Neu)');
    set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
    title(sprintf('Bouton %d - Stratified Speed Profile', targetROI), 'FontWeight', 'normal');
    
    labels = { ...
        sprintf('Low (<%.1f cm/s)', lowThreshold), ...
        sprintf('Med (%.1f-%.1f cm/s)', lowThreshold, highThreshold), ...
        sprintf('High (>%.1f cm/s)', highThreshold) ...
    };
    legend(ax1, labels, 'Location', 'northeast', 'Box', 'off', 'FontSize', 9);
    
    ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.50]); 
    imagesc(1:numPosBins, speedCenters, tuningSurface);
    set(ax2, 'YDir', 'normal'); 
    colormap(parula);
    
    activeData = tuningSurface(tuningSurface > 0);
    if ~isempty(activeData)
        maxVal = prctile(activeData, 98.5);
    else
        maxVal = 1;
    end
    set(ax2, 'CLim', [0, maxVal]);
    
    c = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.025, 0.50]);
    c.Label.String = '\DeltaF/F [NeuC]';
    
    hold on;
    for c_vert = landmarkCentres
        if c_vert <= numPosBins
            xline(c_vert, '--w', 'Alpha', 0.3, 'LineWidth', 1.2);
        end
    end
    
    xticks([1 40 80 120 160 200]);
    xticklabels({'1', '40', '80', '120', '160', '200'});
    xlim([1, numPosBins]);
    
    set(ax2, 'YScale', 'log'); 
    set(ax2, 'YMinorTick', 'on', 'TickDir', 'out');
    yticks([2, 5, 10, 20, 30]); 
    yticklabels({'2', '5', '10', '20', '30'});
    minY = min(speedCenters);
    maxY = max(speedCenters);
    ylim([minY, maxY]);
    
    xlabel('Position (cm)');
    ylabel('Running speed (cm/s)');
    set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
    
    axStrat = axes('Position', [leftMargin + plotWidth + 0.11, 0.12, 0.02, 0.50]);
    hold on;
    
    patch([0 1 1 0], [minY minY lowThreshold lowThreshold], [0.68, 0.92, 0.98], 'EdgeColor', 'none');
    patch([0 1 1 0], [lowThreshold lowThreshold highThreshold highThreshold], [0.88, 0.88, 0.88], 'EdgeColor', 'none');
    patch([0 1 1 0], [highThreshold highThreshold maxY maxY], [0.96, 0.64, 0.76], 'EdgeColor', 'none');
    
    set(axStrat, 'YScale', 'log', 'YLim', [minY, maxY], 'XLim', [0, 1]);
    set(axStrat, 'YTick', [mean([minY, lowThreshold]), mean([lowThreshold, highThreshold]), mean([highThreshold, maxY])], ...
                 'YTickLabel', {'Low', 'Med', 'High'}, 'YAxisLocation', 'right', ...
                 'XTick', [], 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
    ylabel(axStrat, 'Stratification Zones', 'FontSize', 10);
    
    linkaxes([ax1, ax2], 'x');
end