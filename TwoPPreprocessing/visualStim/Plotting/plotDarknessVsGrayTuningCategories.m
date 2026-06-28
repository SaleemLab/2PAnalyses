function plotDarknessVsGrayTuning(responseDark, responseGray, useField)
% Finds the highest-R2 ROIs in Darkness for all 4 template categories
% and plots them with their corresponding Gray Screen profiles superimposed.
% Handles sessions with mismatched speed bin counts cleanly.
%
% USAGE:
%   plotDarknessVsGrayTuning(responseDark, responseGray, 'dFFNeuropilCorrected');

    if nargin < 3, useField = 'dFFNeuropilCorrected'; end
    
    % Core Gaussian Template Function
    gaussFun = @(params, xdata) params(1) + params(2) .* exp(-(((xdata - params(3)).^2) / (2 * (params(4).^2))));

    % 1. Isolate the Best Representative ROIs based on Darkness R2
    classDark   = responseDark.tuningCurve.(useField).classification;
    allTypes    = classDark.tuningType;
    r2Dark      = classDark.R2;
    
    categories = {'lowpass', 'highpass', 'bandpass', 'trough_inverted'};
    plotLabels = {'Low-Pass Archetype', ...
                  'High-Pass Archetype', ...
                  'Band-Pass Archetype', ...
                  'Trough / Inverted Archetype'};
    
    targetROIs = nan(1, 4);
    for c = 1:4
        idxList = find(strcmp(allTypes, categories{c}));
        if ~isempty(idxList)
            [~, bestMatch] = max(r2Dark(idxList));
            targetROIs(c) = idxList(bestMatch);
        end
    end

    % 2. Initialize Figure Window
    fig = figure('Name', 'Speed Tuning: Darkness vs. Gray Screen Superimposed', ...
                 'Position', [30, 150, 1600, 420], 'Color', 'w');
             
    colWidths = 0.22; colGap = 0.018; leftStart = 0.03;

    for c = 1:4
        r = targetROIs(c);
        xPos = leftStart + (c-1) * (colWidths + colGap);
        
        pPanel = uipanel(fig, 'Position', [xPos, 0.02, colWidths, 0.95], ...
                             'BackgroundColor', 'w', 'BorderType', 'none');
        
        if isempty(r) || isnan(r)
            axEmpty = axes(pPanel, 'Position', [0 0 1 1]);
            text(axEmpty, 0.5, 0.5, 'Category Empty', 'HorizontalAlignment', 'center');
            axis off; continue;
        end
        
        % Subsplit panel: 2 parts Stationary, 8 parts Locomotion
        tloInner = tiledlayout(pPanel, 1, 10, 'TileSpacing', 'none', 'Padding', 'none');
        
        %% --- EXTRACT DATA & SPEED AXIS: DARKNESS ---
        edgesD       = responseDark.tuningCurve.speedBins;
        centersD     = edgesD(1:end-1) + diff(edgesD)/2;
        
        y_statD      = responseDark.tuningCurve.(useField).statMean(r);
        y_statErrD   = responseDark.tuningCurve.(useField).statSEM(r);
        y_moveD      = responseDark.tuningCurve.(useField).moveMean(r, :);
        y_moveErrD   = responseDark.tuningCurve.(useField).moveSEM(r, :);
        paramsD      = responseDark.tuningCurve.(useField).classification.fitParams(r, :);
        r2ValD       = responseDark.tuningCurve.(useField).classification.R2(r);
        
        %% --- EXTRACT DATA & SPEED AXIS: GRAY SCREEN ---
        edgesG       = responseGray.tuningCurve.speedBins;
        centersG     = edgesG(1:end-1) + diff(edgesG)/2;
        
        y_statG      = responseGray.tuningCurve.(useField).statMean(r);
        y_statErrG   = responseGray.tuningCurve.(useField).statSEM(r);
        y_moveG      = responseGray.tuningCurve.(useField).moveMean(r, :);
        y_moveErrG   = responseGray.tuningCurve.(useField).moveSEM(r, :);
        paramsG      = responseGray.tuningCurve.(useField).classification.fitParams(r, :);
        r2ValG       = responseGray.tuningCurve.(useField).classification.R2(r);
        
        % Clean drop of any behavioral velocity NaNs per session
        validD = ~isnan(y_moveD); centersD = centersD(validD); y_moveD = y_moveD(validD); y_moveErrD = y_moveErrD(validD);
        validG = ~isnan(y_moveG); centersG = centersG(validG); y_moveG = y_moveG(validG); y_moveErrG = y_moveErrG(validG);
        
        % Global Y scaling boundary calculator
        allY = [y_statD+y_statErrD, y_statD-y_statErrD, y_moveD+y_moveErrD, y_moveD-y_moveErrD, ...
                y_statG+y_statErrG, y_statG-y_statErrG, y_moveG+y_moveErrG, y_moveG-y_moveErrG];
        yMin = min(allY); yMax = max(allY); yPad = (yMax - yMin) * 0.15; if yPad <= 1e-5, yPad = 1; end
        yLimits = [yMin - yPad, yMax + yPad];

        % Color Schemes
        if c == 1 || c == 4
            colorDark = [0.85, 0.33, 0.1];  colorGray = [0.93, 0.69, 0.13]; 
        else
            colorDark = [0, 0.45, 0.74];    colorGray = [0.30, 0.75, 0.93]; 
        end

        %% --- AXIS 1: DISCONNECTED STATIONARY SLOTS (X=0) ---
        axStat = nexttile(tloInner, 1, [1, 2]); hold(axStat, 'on');
        
        errorbar(axStat, -0.12, y_statD, y_statErrD, 'ok', 'MarkerSize', 6, ...
            'LineWidth', 1.2, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', colorDark);
        
        errorbar(axStat, 0.12, y_statG, y_statErrG, 'ok', 'MarkerSize', 6, ...
            'LineWidth', 1.2, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', colorGray);
        
        xlim(axStat, [-0.5, 0.5]); ylim(axStat, yLimits);
        axStat.XTick = 0; axStat.XTickLabel = {'0'};
        set(axStat, 'TickDir', 'out', 'LineWidth', 1.1, 'FontSize', 10);
        box(axStat, 'off'); grid(axStat, 'off');
        if c > 1, ylabel(axStat, ''); else ylabel(axStat, '\DeltaF/F (Neu)', 'FontSize', 12); end

        %% --- AXIS 2: LOCOMOTION DOMAIN (Log Scaling) ---
        axMove = nexttile(tloInner, 3, [1, 8]); hold(axMove, 'on');
        
        % FIXED: Plots using session-specific centers vectors now to eliminate size conflicts
        hD = errorbar(axMove, centersD, y_moveD, y_moveErrD, 'ok', 'MarkerFaceColor', colorDark, ...
            'MarkerEdgeColor', colorDark, 'MarkerSize', 5, 'LineWidth', 1.1);
        
        hG = errorbar(axMove, centersG, y_moveG, y_moveErrG, 'ok', 'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', colorGray, 'MarkerSize', 5, 'LineWidth', 1.1);
        
        % Generate high-density line overlays using continuous independent ranges
        fineAxisD = logspace(log10(min(centersD)), log10(max(centersD)), 200);
        y_fitD = gaussFun(paramsD, fineAxisD);
        plot(axMove, fineAxisD, y_fitD, 'Color', colorDark, 'LineWidth', 2.2);
        
        fineAxisG = logspace(log10(min(centersG)), log10(max(centersG)), 200);
        y_fitG = gaussFun(paramsG, fineAxisG);
        plot(axMove, fineAxisG, y_fitG, '--', 'Color', colorGray, 'LineWidth', 2.2);
        
        % Axis Scale Setup
        set(axMove, 'XScale', 'log');
        globalMin = min([centersD, centersG]); globalMax = max([centersD, centersG]);
        xlim(axMove, [globalMin*0.85, globalMax*1.15]); ylim(axMove, yLimits);
        
        axMove.XTick = [2, 5, 10, 20, 40]; axMove.XTickLabel = {'2', '5', '10', '20', '40'};
        set(axMove, 'TickDir', 'out', 'LineWidth', 1.1, 'FontSize', 10, 'YTickLabel', '');
        axMove.YAxis.Visible = 'off'; box(axMove, 'off');
        grid(axMove, 'on'); axMove.XGrid = 'on'; axMove.YGrid = 'off';
        
        % Header Information
        titleString = sprintf('\\bf%s (ROI %d)\\newline\\rm\\fontsize{9}Dark R^2 = %.2f | Gray R^2 = %.2f', ...
            plotLabels{c}, r, r2ValD, r2ValG);
        title(axMove, titleString, 'FontWeight', 'normal', 'Interpreter', 'tex', 'HorizontalAlignment', 'center');
        xlabel(axMove, 'Speed (cm/s)', 'FontSize', 11);
        
        if c == 1
            legend([hD, hG], {'Darkness', 'Gray Screen'}, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
        end
    end
end