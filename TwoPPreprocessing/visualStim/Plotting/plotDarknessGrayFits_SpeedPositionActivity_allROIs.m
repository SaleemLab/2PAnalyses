function plotDarknessGrayFits_SpeedPositionActivity_allROIs(response_darkness, response_gray, respVR, pdfFullPath, useField, applySmoothing, smoothSigma)
% Iterates through EVERY single ROI in your dataset, creating a dedicated 
% 3-panel horizontal dashboard page per bouton, and appends them sequentially
% into a single multi-page master PDF document using the modern export engine.
%
% COLUMNS (Per Page):
%   Left Panel:   Darkness 1D Tuning Profile (Black + Fit Line + Shuffle p-val)
%   Center Panel: Gray Screen 1D Tuning Profile (Gray + Fit Line + Shuffle p-val)
%   Right Panel:  VR 2D Heatmap & 1D Speed-Stratified Spatial Profiles
%
% VISUAL ASSIGNMENT NOTE:
%   Panel titles automatically apply a soft green background fill if the specific 
%   condition meets your 99% significance criteria (p <= 0.01).
    if nargin < 4 || isempty(pdfFullPath)
        error('Please specify a target output path string ending in .pdf');
    end
    if nargin < 5, useField = 'dFFNeuropilCorrected'; end
    if nargin < 6, applySmoothing = true; end
    if nargin < 7, smoothSigma = [1.1, 1.5]; end
    
    [targetDir, ~, ~] = fileparts(pdfFullPath);
    if ~isempty(targetDir) && ~exist(targetDir, 'dir')
        mkdir(targetDir);
    end
    
    if exist(pdfFullPath, 'file')
        delete(pdfFullPath);
    end
    totalROIs = size(respVR.speedPositionActivity.matrix, 3);
    fprintf('=== Multi-Page PDF Engine Initialized ===\nCompiling %d total boutons. Master file: %s\n\n', totalROIs, pdfFullPath);
    
    FIXED_LOW_CUTOFF  = 15; 
    FIXED_HIGH_CUTOFF = 27; 
    landmarkCentres   = [40, 80, 120, 160];
    gaussFun = @(params, xdata) params(1) + params(2) .* exp(-(((xdata - params(3)).^2) / (2 * (params(4).^2))));
    
    colorDarkLine = [0, 0, 0];          
    colorGrayLine = [0.55, 0.55, 0.55];  
    
    % Color specifications for the text box titles
    sigGreenBackground  = [0.85, 0.95, 0.85]; % Soft green background for p <= 0.01
    defaultWhiteBackground = [1, 1, 1];        % Plain white background
    
    for targetROI = 1:totalROIs
        fprintf('Appending page for Bouton %d/%d...\n', targetROI, totalROIs);
        
        figHandle = figure('Name', sprintf('ROI %d PDF Canvas', targetROI), ...
                           'Position', [50, 50, 1500, 420], 'Color', 'w', 'Visible', 'off');
        
        vOffset = 0.15;
        rowHeights = 0.70;
        
        % PANEL 1: DARKNESS 1D TUNING CURVE
        edgesD   = response_darkness.tuningCurve.speedBins;
        centersD = edgesD(1:end-1) + diff(edgesD)/2;
        y_statD  = response_darkness.tuningCurve.(useField).statMean(targetROI);
        y_statED = response_darkness.tuningCurve.(useField).statSEM(targetROI);
        y_moveD  = response_darkness.tuningCurve.(useField).moveMean(targetROI, :);
        y_moveED = response_darkness.tuningCurve.(useField).moveSEM(targetROI, :);
        paramsD  = response_darkness.tuningCurve.(useField).classification.fitParams(targetROI, :);
        r2D      = response_darkness.tuningCurve.(useField).classification.R2(targetROI);
        typeD    = response_darkness.tuningCurve.(useField).classification.tuningType{targetROI};
        
        pValD = NaN;
        if isfield(response_darkness.tuningCurve.(useField), 'pValMoving')
            pValD = response_darkness.tuningCurve.(useField).pValMoving(targetROI);
        end
        if isnan(pValD), pStrD = 'NaN'; else, pStrD = sprintf('%.4f', pValD); end
        
        validD = ~isnan(y_moveD); centersD = centersD(validD); y_moveD = y_moveD(validD); y_moveED = y_moveED(validD);
        
        % PANEL 2: GRAY SCREEN 1D TUNING CURVE
        edgesG   = response_gray.tuningCurve.speedBins;
        centersG = edgesG(1:end-1) + diff(edgesG)/2;
        y_statG  = response_gray.tuningCurve.(useField).statMean(targetROI);
        y_statEG = response_gray.tuningCurve.(useField).statSEM(targetROI);
        y_moveG  = response_gray.tuningCurve.(useField).moveMean(targetROI, :);
        y_moveEG = response_gray.tuningCurve.(useField).moveSEM(targetROI, :);
        paramsG  = response_gray.tuningCurve.(useField).classification.fitParams(targetROI, :);
        r2G      = response_gray.tuningCurve.(useField).classification.R2(targetROI);
        typeG    = response_gray.tuningCurve.(useField).classification.tuningType{targetROI};
        
        pValG = NaN;
        if isfield(response_gray.tuningCurve.(useField), 'pValMoving')
            pValG = response_gray.tuningCurve.(useField).pValMoving(targetROI);
        end
        if isnan(pValG), pStrG = 'NaN'; else, pStrG = sprintf('%.4f', pValG); end
        
        validG = ~isnan(y_moveG); centersG = centersG(validG); y_moveG = y_moveG(validG); y_moveEG = y_moveEG(validG);
        
        all1DY = [y_statD, y_moveD, y_statG, y_moveG];
        yLimits1D = [min(all1DY)*0.85, max(all1DY)*1.15];
        
        % --- RENDERING PANEL 1: DARKNESS ---
        axStatD = axes('Position', [0.04, vOffset, 0.02, rowHeights]); hold(axStatD, 'on');
        axMoveD = axes('Position', [0.065, vOffset, 0.18, rowHeights]); hold(axMoveD, 'on');
        
        errorbar(axStatD, 0, y_statD, y_statED, 'ok', 'MarkerSize', 5.5, 'LineWidth', 1.1, 'MarkerEdgeColor', colorDarkLine);
        xlim(axStatD, [-0.5, 0.5]); axStatD.XTick = 0; axStatD.XTickLabel = {'0'}; ylim(axStatD, yLimits1D);
        set(axStatD, 'TickDir', 'out', 'Box', 'off'); ylabel(axStatD, '\DeltaF/F (Neu)', 'FontWeight', 'bold');
        
        errorbar(axMoveD, centersD, y_moveD, y_moveED, 'ok', 'MarkerFaceColor', colorDarkLine, 'MarkerEdgeColor', colorDarkLine, 'MarkerSize', 4.5);
        fineXD = logspace(log10(min(centersD)), log10(max(centersD)), 200);
        plot(axMoveD, fineXD, gaussFun(paramsD, fineXD), 'Color', colorDarkLine, 'LineWidth', 2.2);
        set(axMoveD, 'XScale', 'log', 'YTickLabel', '', 'Box', 'off', 'TickDir', 'out'); axMoveD.YAxis.Visible = 'off';
        
        axMoveD.XTick = [2, 5, 10, 20, 30, 40]; 
        axMoveD.XTickLabel = {'2', '5', '10', '20', '30', '40'}; 
        ylim(axMoveD, yLimits1D); grid(axMoveD, 'on');
        xlabel(axMoveD, 'Speed (cm/s)');
        
        % Assign Darkness Title Text and update Background if Significant
        titleTextD = sprintf('\\bfDarkness Tuning \\rm(%s)\\newline\\bfR^2 = %.3f | p_{shuffle} = %s', upper(typeD), r2D, pStrD);
        tH_D = title(axMoveD, titleTextD, 'Interpreter', 'tex');
        if ~isnan(pValD) && pValD <= 0.01
            set(tH_D, 'BackgroundColor', sigGreenBackground, 'EdgeColor', [0.4, 0.7, 0.4], 'Margin', 3);
        else
            set(tH_D, 'BackgroundColor', defaultWhiteBackground, 'EdgeColor', 'none');
        end
        
        % --- RENDERING PANEL 2: GRAY SCREEN ---
        axStatG = axes('Position', [0.29, vOffset, 0.02, rowHeights]); hold(axStatG, 'on');
        axMoveG = axes('Position', [0.315, vOffset, 0.18, rowHeights]); hold(axMoveG, 'on');
        
        errorbar(axStatG, 0, y_statG, y_statEG, 'ok', 'MarkerSize', 5.5, 'LineWidth', 1.1, 'MarkerEdgeColor', colorGrayLine);
        xlim(axStatG, [-0.5, 0.5]); axStatG.XTick = 0; axStatG.XTickLabel = {'0'}; ylim(axStatG, yLimits1D);
        set(axStatG, 'TickDir', 'out', 'Box', 'off', 'YTickLabel', '');
        
        errorbar(axMoveG, centersG, y_moveG, y_moveEG, 'ok', 'MarkerFaceColor', colorGrayLine, 'MarkerEdgeColor', colorGrayLine, 'MarkerSize', 4.5);
        fineXG = logspace(log10(min(centersG)), log10(max(centersG)), 200);
        plot(axMoveG, fineXG, gaussFun(paramsG, fineXG), 'Color', colorGrayLine, 'LineWidth', 2.2);
        set(axMoveG, 'XScale', 'log', 'YTickLabel', '', 'Box', 'off', 'TickDir', 'out'); axMoveG.YAxis.Visible = 'off';
        
        axMoveG.XTick = [2, 5, 10, 20, 30, 40]; 
        axMoveG.XTickLabel = {'2', '5', '10', '20', '30', '40'}; 
        ylim(axMoveG, yLimits1D); grid(axMoveG, 'on');
        xlabel(axMoveG, 'Speed (cm/s)');
        
        % Assign Gray Screen Title Text and update Background if Significant
        titleTextG = sprintf('\\bfGray Screen Tuning \\rm(%s)\\newline\\bfR^2 = %.3f | p_{shuffle} = %s', upper(typeG), r2G, pStrG);
        tH_G = title(axMoveG, titleTextG, 'Interpreter', 'tex');
        if ~isnan(pValG) && pValG <= 0.01
            set(tH_G, 'BackgroundColor', sigGreenBackground, 'EdgeColor', [0.4, 0.7, 0.4], 'Margin', 3);
        else
            set(tH_G, 'BackgroundColor', defaultWhiteBackground, 'EdgeColor', 'none');
        end
        
        % PANEL 3: VR POSITION-SPEED ACTIVITY HEATMAP & PROFILES
        tuningSurface = respVR.speedPositionActivity.matrix(:, :, targetROI);
        speedCenters  = respVR.speedPositionActivity.speedBinCenters;
        speedEdges    = respVR.speedPositionActivity.speedEdges;
        numPosBins    = size(tuningSurface, 2);
        x_pos         = 1:numPosBins;
        
        if applySmoothing
            mask = ~isnan(tuningSurface); dataZeroed = tuningSurface; dataZeroed(~mask) = 0;
            bData = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
            bMask = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
            tuningSurface = bData ./ bMask; tuningSurface(isnan(tuningSurface)) = 0;
        else
            tuningSurface(isnan(tuningSurface)) = 0;
        end
        
        lowRows  = find(speedCenters < FIXED_LOW_CUTOFF);
        medRows  = find(speedCenters >= FIXED_LOW_CUTOFF & speedCenters <= FIXED_HIGH_CUTOFF);
        highRows = find(speedCenters > FIXED_HIGH_CUTOFF);
        
        [lowL, lowE]   = getProfileStats(tuningSurface, lowRows);
        [medL, medE]   = getProfileStats(tuningSurface, medRows);
        [highL, highE] = getProfileStats(tuningSurface, highRows);
        
        if ~isempty(lowRows),  lowThresh  = speedEdges(max(lowRows) + 1);  else, lowThresh  = min(speedEdges); end
        if ~isempty(medRows),  highThresh = speedEdges(max(medRows) + 1);  else, highThresh = lowThresh; end
        minY = min(speedCenters); maxY = max(speedCenters);
        
        axProfVR = axes('Position', [0.55, vOffset + (rowHeights*0.58), 0.35, rowHeights*0.40]); hold(axProfVR, 'on');
        for cLines = landmarkCentres, xline(axProfVR, cLines, ':', 'Color', [0.6 0.6 0.6]); end
        renderShadedError(axProfVR, x_pos, lowL, lowE, [0, 0.75, 1]);
        renderShadedError(axProfVR, x_pos, medL, medE, [0.2, 0.2, 0.2]);
        renderShadedError(axProfVR, x_pos, highL, highE, [1, 0, 0.6]);
        plot(axProfVR, x_pos, lowL, 'Color', [0, 0.75, 1], 'LineWidth', 1.5);
        plot(axProfVR, x_pos, medL, 'Color', [0.2, 0.2, 0.2], 'LineWidth', 1.5);
        plot(axProfVR, x_pos, highL, 'Color', [1, 0, 0.6], 'LineWidth', 1.5);
        set(axProfVR, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
        ylabel(axProfVR, '\DeltaF/F');
        title(axProfVR, sprintf('\\bfBouton %d Pipeline Diagnostic Tracking (VR)', targetROI), 'FontSize', 11);
        
        axMapVR = axes('Position', [0.55, vOffset, 0.35, rowHeights*0.52]);
        imagesc(axMapVR, 1:numPosBins, speedCenters, tuningSurface);
        set(axMapVR, 'YDir', 'normal'); colormap(axMapVR, parula);
        set(axMapVR, 'XScale', 'linear', 'YScale', 'log', 'YMinorTick', 'on', 'TickDir', 'out');
        yticks(axMapVR, [2, 5, 10, 20, 30, 40]); yticklabels(axMapVR, {'2', '5', '10', '20', '30', '40'});
        xlim(axMapVR, [1, numPosBins]); ylim(axMapVR, [minY, maxY]);
        xticks(axMapVR, [1 40 80 120 160 200]); xticklabels(axMapVR, {'1', '40', '80', '120', '160', '200'});
        xlabel(axMapVR, 'Corridor Position (cm)'); ylabel(axMapVR, 'cm/s'); box(axMapVR, 'off');
        
        activeD = tuningSurface(tuningSurface > 0); maxV = 1; if ~isempty(activeD), maxV = prctile(activeD, 98.5); end
        set(axMapVR, 'CLim', [0, maxV]);
        
        axKeyVR = axes('Position', [0.915, vOffset, 0.01, rowHeights*0.52]); hold(axKeyVR, 'on');
        patch([0 1 1 0], [minY minY lowThresh lowThresh], [0.68, 0.92, 0.98], 'EdgeColor', 'none');
        patch([0 1 1 0], [lowThresh lowThresh highThresh highThresh], [0.88, 0.88, 0.88], 'EdgeColor', 'none');
        patch([0 1 1 0], [highThresh highThresh maxY maxY], [0.96, 0.64, 0.76], 'EdgeColor', 'none');
        set(axKeyVR, 'YScale', 'log', 'YLim', [minY, maxY], 'XLim', [0, 1], 'XTick', [], 'Box', 'off');
        set(axKeyVR, 'YTick', [mean([minY, lowThresh]), mean([lowThresh, highThresh]), mean([highThresh, maxY])], ...
                     'YTickLabel', {'L', 'M', 'H'}, 'YAxisLocation', 'right', 'FontSize', 7);
        
        linkaxes([axProfVR, axMapVR], 'x');
        
        if targetROI == 1
            exportgraphics(figHandle, pdfFullPath, 'ContentType', 'vector');
        else
            exportgraphics(figHandle, pdfFullPath, 'ContentType', 'vector', 'Append', true);
        end
        close(figHandle);
    end
    fprintf('\n=== Inspection Compilation Completed! ===\nManual updated with green significance tiles at:\n%s\n', pdfFullPath);
end

function [mLine, semLine] = getProfileStats(surface, rowIndices)
    if isempty(rowIndices), mLine = nan(1, size(surface, 2)); semLine = nan(1, size(surface, 2)); return; end
    mLine = mean(surface(rowIndices, :), 1, 'omitnan');
    nRows = length(rowIndices);
    if nRows > 1, semLine = std(surface(rowIndices, :), 0, 1, 'omitnan') ./ sqrt(nRows); else, semLine = zeros(1, size(surface, 2)); end
end

function renderShadedError(ax, x, m, err, col)
    if all(isnan(m)) || all(err == 0), return; end
    x_p = [x, fliplr(x)]; y_p = [(m + err), fliplr(m - err)];
    invalid = isnan(y_p); x_p(invalid) = []; y_p(invalid) = [];
    patch(ax, x_p, y_p, col, 'EdgeColor', 'none', 'FaceAlpha', 0.11);
end