function saveSpeedPositionAndTuningSummaries(sessionFileInfo, response, respGray, respDark, options)
    if nargin < 5, options = struct(); end
    if ~isfield(options, 'applySmoothing'), options.applySmoothing = true; end
    if ~isfield(options, 'smoothSigma'),    options.smoothSigma = [1.1, 1.5]; end
    
    figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'ComprehensiveSpeedSummaries');
    if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
    
    pdfName = sprintf('%s_%s_AllROIs_ComprehensiveSpeed.pdf', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name);
    fullPDFPath = fullfile(figSaveDir, pdfName);
    if exist(fullPDFPath, 'file'), delete(fullPDFPath); end
    
    numROIs = size(response.speedPositionActivity.matrix, 3);
    fprintf('Starting master PDF generation for %d ROIs...\n', numROIs);
    
    for targetROI = 1:numROIs
        figHandle = plotComprehensiveSpeed_Internal(sessionFileInfo, response, respGray, respDark, targetROI, ...
            options.applySmoothing, options.smoothSigma);
        
        exportgraphics(figHandle, fullPDFPath, 'ContentType', 'vector', 'Append', true);
        close(figHandle);
        
        if mod(targetROI, 10) == 0 || targetROI == numROIs
            fprintf('Processed %d/%d ROIs...\n', targetROI, numROIs);
        end
    end
    fprintf('Success! All ROIs saved to: %s\n', fullPDFPath);
end

function figHandle = plotComprehensiveSpeed_Internal(sessionFileInfo, response, respGray, respDark, targetROI, applySmoothing, smoothSigma)
    % 1. Extract VR Data Matrices
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
    
    % Check for matching Darkness / Gray inputs
    hasDark = ~isempty(respDark) && isfield(respDark, 'tuningCurve');
    hasGray = ~isempty(respGray) && isfield(respGray, 'tuningCurve');
    
    % Find correct neural signal field name (e.g., dFFNeuropilCorrected, dFF, or spks)
    targetFields = {'dFFNeuropilCorrected', 'dFF', 'spks', 'spikes'};
    signalField = 'dFF'; % Default fallback
    if hasGray
        fNames = fieldnames(respGray.tuningCurve);
        matched = intersect(fNames, targetFields, 'stable');
        if ~isempty(matched), signalField = matched{1}; end
    end
    
    % 2. Setup Figure Layout Window
    figHandle = figure('Name', sprintf('Bouton %d Comprehensive', targetROI), ...
        'Position', [50 50 1100 620], 'Color', 'w', 'Visible', 'off');
    
    % Spatial Geometry Dimensions
    leftMargin = 0.08;
    plotWidth  = 0.50; 
    tuningLeft = 0.72;
    tuningWidth = 0.22;
    
    %% --- Subplot 1: Stratified Position Curves (Top-Left) ---
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
    title(sprintf('Bouton %d - Virtual Reality Spatial Dynamics (Mismatched Gain)', targetROI), 'FontWeight', 'bold');
    
    labels = { ...
        sprintf('Low (<%.1f cm/s)', lowThreshold), ...
        sprintf('Med (%.1f-%.1f cm/s)', lowThreshold, highThreshold), ...
        sprintf('High (>%.1f cm/s)', highThreshold) ...
    };
    legend(ax1, labels, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
    
    %% --- Subplot 2: 2D Position-Speed Heatmap (Bottom-Left) ---
    ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.50]); 
    imagesc(1:numPosBins, speedCenters, tuningSurface);
    set(ax2, 'YDir', 'normal'); 
    colormap(ax2, parula);
    
    activeData = tuningSurface(tuningSurface > 0);
    maxVal = 1; if ~isempty(activeData), maxVal = prctile(activeData, 98.5); end
    set(ax2, 'CLim', [0, maxVal]);
    
    c = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.015, 0.12, 0.018, 0.50]);
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
    
    set(ax2, 'YScale', 'log', 'YMinorTick', 'on', 'TickDir', 'out');
    yticks([2, 5, 10, 20, 30]); yticklabels({'2', '5', '10', '20', '30'});
    minY = min(speedCenters); maxY = max(speedCenters);
    ylim([minY, maxY]);
    xlabel('Position (cm)'); ylabel('VR Virtual Speed (cm/s)');
    set(ax2, 'Box', 'off', 'FontSize', 10);
    
    %% --- Subplot 3: Small Vertical Stratification Zone Strip ---
    axStrat = axes('Position', [leftMargin + plotWidth + 0.075, 0.12, 0.015, 0.50]);
    hold on;
    patch([0 1 1 0], [minY minY lowThreshold lowThreshold], [0.68, 0.92, 0.98], 'EdgeColor', 'none');
    patch([0 1 1 0], [lowThreshold lowThreshold highThreshold highThreshold], [0.88, 0.88, 0.88], 'EdgeColor', 'none');
    patch([0 1 1 0], [highThreshold highThreshold maxY maxY], [0.96, 0.64, 0.76], 'EdgeColor', 'none');
    
    set(axStrat, 'YScale', 'log', 'YLim', [minY, maxY], 'XLim', [0, 1]);
    set(axStrat, 'YTick', [], 'XTick', [], 'Box', 'off', 'TickDir', 'out');
    
    %% --- Subplot 4: Non-Spatial Speed Tuning Curves (Far Right Panel) ---
    axTuning = subplot('Position', [tuningLeft, 0.12, tuningWidth, 0.80]);
    hold on;
    
    % Dynamic Title text generation based on significance shuffles
    sigText = 'Speed Profile';
    if hasGray
        sigG = respGray.tuningCurve.(signalField).isSignificant_999(targetROI);
        colG = '\color[rgb]{0, 0.5, 0}Significant'; if ~sigG, colG = '\color[rgb]{0.7, 0, 0}Not Sig'; end
        sigText = sprintf('Gray: %s', colG);
    end
    if hasDark
        sigD = respDark.tuningCurve.(signalField).isSignificant_999(targetROI);
        colD = '\color[rgb]{0, 0.5, 0}Significant'; if ~sigD, colD = '\color[rgb]{0.7, 0, 0}Not Sig'; end
        sigText = sprintf('%s\\color{black} | Dark: %s', sigText, colD);
    end
    title(axTuning, sigText, 'FontSize', 10, 'Interpreter', 'tex');
    
    % Render Gray Screen Line Layer
    if hasGray
        tcG = respGray.tuningCurve;
        x_valsG = [0.5, tcG.speedBins(1:end-1) + diff(tcG.speedBins)/2];
        y_meanG = [tcG.(signalField).statMean(targetROI); tcG.(signalField).moveMean(targetROI, :)'];
        y_errG  = [tcG.(signalField).statSEM(targetROI);  tcG.(signalField).moveSEM(targetROI, :)'];
        errorbar(axTuning, x_valsG, y_meanG, y_errG, 'Color', 'k', 'LineStyle', '-', ...
            'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    end
    
    % Render Darkness Line Layer
    if hasDark
        tcD = respDark.tuningCurve;
        x_valsD = [0.5, tcD.speedBins(1:end-1) + diff(tcD.speedBins)/2];
        y_meanD = [tcD.(signalField).statMean(targetROI); tcD.(signalField).moveMean(targetROI, :)'];
        y_errD  = [tcD.(signalField).statSEM(targetROI);  tcD.(signalField).moveSEM(targetROI, :)'];
        errorbar(axTuning, x_valsD, y_meanD, y_errD, 'Color', [0.4 0.4 0.4], 'LineStyle', '--', ...
            'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', 5, 'MarkerFaceColor', [0.4 0.4 0.4]);
    end
    
    % Formatting Non-Spatial Axis
    set(axTuning, 'XScale', 'log', 'Box', 'off', 'TickDir', 'out', 'YAxisLocation', 'right');
    grid on;
    xticks(axTuning, [0.5, 1, 5, 10, 20, 40, 60, 100]);
    xticklabels(axTuning, {'0','1','5','10','20','40','60','100'});
    xlim(axTuning, [0.4, 110]);
    xlabel('Physical Wheel Speed (cm/s)', 'FontWeight', 'bold');
    ylabel(axTuning, sprintf('%s Response Value', signalField), 'Interpreter', 'none');
    
    % Add unique Legend configuration for the side pane
    if hasGray && hasDark
        legend(axTuning, {'Gray Screen (Open Loop)', 'Darkness (Pure Motor)'}, 'Location', 'southoutside', 'Box', 'off');
    elseif hasGray
        legend(axTuning, {'Gray Screen'}, 'Location', 'southoutside', 'Box', 'off');
    end
    
    linkaxes([ax1, ax2], 'x');
end