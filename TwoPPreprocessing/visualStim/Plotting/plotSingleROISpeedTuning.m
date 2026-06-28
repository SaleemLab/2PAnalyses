function figHandle = plotSingleROISpeedTuning(response, roiIds, useField)
    % roiIds:  exactly 4 ROI indices e.g. [48, 129, 70, 233]
    % this version plots the gaussian template fits and 
    if nargin < 3, useField = 'dFFNeuropilCorrected'; end

    edges         = response.tuningCurve.speedBins;
    movingCenters = edges(1:end-1) + diff(edges)/2;

    statMean = response.tuningCurve.(useField).statMean;
    statSEM  = response.tuningCurve.(useField).statSEM;
    moveMean = response.tuningCurve.(useField).moveMean;
    moveSEM  = response.tuningCurve.(useField).moveSEM;

    dataStruct = response.tuningCurve.(useField);
    if isfield(dataStruct, 'classification')
        classStruct = dataStruct.classification;
        isLog2      = false;
    end

    r2Values      = classStruct.R2;
    fitParams     = classStruct.fitParams;
    allTypes      = classStruct.tuningType;

    if isfield(dataStruct, 'pValMoving')
        pValuesMoving = dataStruct.pValMoving;
    else
        pValuesMoving = zeros(size(r2Values));
    end

    gaussFun = @(params, xdata) params(1) + params(2) .* ...
        exp(-(((xdata - params(3)).^2) / (2 * (params(4).^2))));

    fineSpeedAxis_linear = logspace(log10(min(movingCenters)), log10(max(movingCenters)), 200);
    if isLog2
        fineSpeedAxis_fit = log2(fineSpeedAxis_linear);
    else
        fineSpeedAxis_fit = fineSpeedAxis_linear;
    end

    % Same size as plotAllSpeedTuningCategories
    figW_cm = 18;
    figH_cm = 5;
    figHandle = figure('Name', sprintf('ROI %s Speed Tuning', mat2str(roiIds)), ...
        'Units', 'centimeters', 'Position', [2, 2, figW_cm, figH_cm], ...
        'Color', 'w', 'PaperPositionMode', 'auto');

    colWidth    = 0.20;
    colGap      = 0.02;
    leftStart   = 0.04;
    rowHeight   = 0.68;
    bottomStart = 0.20;

    for c = 1:4
        r    = roiIds(c);
        xPos = leftStart + (c-1) * (colWidth + colGap);

        axStat = axes(figHandle, 'Position', [xPos,                 bottomStart, colWidth*0.25, rowHeight]); 
        axMove = axes(figHandle, 'Position', [xPos + colWidth*0.28, bottomStart, colWidth*0.72, rowHeight]); 

        y_stat    = statMean(r);
        y_statErr = statSEM(r);
        y_move    = moveMean(r, :);
        y_moveErr = moveSEM(r, :);
        params    = fitParams(r, :);
        currentR2 = r2Values(r);
        currentP  = pValuesMoving(r);
        roiType   = allTypes{r};

        allY    = [y_stat+y_statErr, y_stat-y_statErr, y_move+y_moveErr, y_move-y_moveErr];
        yLimits = [min(allY)*0.85, max(allY)*1.15];

        % Stationary
        hold(axStat, 'on');
        errorbar(axStat, 0, y_stat, y_statErr, 'ok', 'MarkerSize', 2.5, ...
            'LineWidth', 0.5, 'MarkerFaceColor', 'none');
        xlim(axStat, [-0.5, 0.5]); ylim(axStat, yLimits);
        axStat.XTick = 0; axStat.XTickLabel = {'0'};
        set(axStat, 'TickDir', 'out', 'LineWidth', 0.5, 'FontSize', 7);
        box(axStat, 'off');
        if c == 1
            ylabel(axStat, '\DeltaF/F', 'FontSize', 8);
        end

        % Moving
        hold(axMove, 'on');
        errorbar(axMove, movingCenters, y_move, y_moveErr, 'ok', ...
            'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerSize', 2.5, 'LineWidth', 0.5);

        if all(isfinite(params))
            y_fit = gaussFun(params, fineSpeedAxis_fit);
            lineColor = [0, 0.45, 0.74];
            if ismember(roiType, {'lowpass', 'trough_inverted'})
                lineColor = [0.85, 0.33, 0.1];
            end
            plot(axMove, fineSpeedAxis_linear, y_fit, 'Color', lineColor, 'LineWidth', 1.5);
        end

        set(axMove, 'XScale', 'log');
        xlim(axMove, [min(movingCenters)*0.85, max(movingCenters)*1.15]);
        ylim(axMove, yLimits);
        axMove.XTick      = [2, 5, 10, 20, 40];
        axMove.XTickLabel = {'2', '5', '10', '20', '40'};
        set(axMove, 'TickDir', 'out', 'LineWidth', 0.5, 'FontSize', 7, 'YTickLabel', '');
        axMove.YAxis.Visible = 'off';
        box(axMove, 'off');
        grid(axMove, 'on');

        titleString = sprintf('%s\\newline\\fontsize{7}ROI %d   (R^2 = %.3f | p = %.4f)', ...
            roiType, r, currentR2, currentP);
        title(axMove, titleString, 'FontWeight', 'normal', 'Interpreter', 'tex', ...
            'FontSize', 7, 'HorizontalAlignment', 'center');
        xlabel(axMove, 'Speed (cm/s)', 'FontSize', 7);
    end

    exportgraphics(figHandle, 'speedTuning_selectedROIs.pdf', 'ContentType', 'vector');
end