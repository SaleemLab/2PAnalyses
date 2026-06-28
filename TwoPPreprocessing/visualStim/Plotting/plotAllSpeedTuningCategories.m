% function figHandle = plotAllSpeedTuningCategories(response, useField)
% % Plots representative ROIs for all 4 competitive template classes with 
% % a logarithmic X-axis for moving speeds and a disconnected, unfilled 
% % marker at X=0 to represent the isolated stationary firing rate.
% %
% % USAGE:
% %   plotAllSpeedTuningCategories(response, 'dFFNeuropilCorrected');
% 
%     if nargin < 2, useField = 'dFFNeuropilCorrected'; end
% 
%     p_thresh = 0.01;
% 
%     % 1. Extract structural fields
%     edges         = response.tuningCurve.speedBins;
%     movingCenters = edges(1:end-1) + diff(edges)/2;
% 
%     statMean = response.tuningCurve.(useField).statMean;
%     statSEM  = response.tuningCurve.(useField).statSEM;
%     moveMean = response.tuningCurve.(useField).moveMean;
%     moveSEM  = response.tuningCurve.(useField).moveSEM;
% 
%     % Load log2 classification if available, otherwise fall back to standard
%     dataStruct = response.tuningCurve.(useField);
%     if isfield(dataStruct, 'classification')
%         classStruct = dataStruct.classification;
%         isLog2      = false;
% %     else
% %         classStruct = dataStruct.classification;
% %         isLog2      = false;
%     end
% 
%     allTypes  = classStruct.tuningType;
%     r2Values  = classStruct.R2;
%     fitParams = classStruct.fitParams;
% 
%     if isfield(dataStruct, 'pValMoving')
%         pValuesMoving = dataStruct.pValMoving;
%     else
%         warning('pValMoving not found. Bypassing shuffle check.');
%         pValuesMoving = zeros(size(r2Values));
%     end
% 
%     % 2. Find best representative ROI per category
%     categories = {'lowpass', 'highpass', 'bandpass', 'trough_inverted'};
%     plotLabels = {'Low-Pass', 'High-Pass', 'Band-Pass', 'Trough / Inverted'};
% 
%     targetROIs = nan(1, 4);
%     for c = 1:4
%         idxList = find(strcmp(allTypes, categories{c}) & (pValuesMoving <= p_thresh));
%         if ~isempty(idxList)
%             [~, bestMatch] = max(r2Values(idxList));
%             targetROIs(c)  = idxList(bestMatch);
%         else
%             idxListFallback = find(strcmp(allTypes, categories{c}));
%             if ~isempty(idxListFallback)
%                 [~, bestMatch] = max(r2Values(idxListFallback));
%                 targetROIs(c)  = idxListFallback(bestMatch);
%                 warning('Category %s: no cells pass p<=%.3f. Plotting unvalidated best fit.', categories{c}, p_thresh);
%             end
%         end
%     end
% 
%     gaussFun = @(params, xdata) params(1) + params(2) .* ...
%         exp(-(((xdata - params(3)).^2) / (2 * (params(4).^2))));
% 
%     % 3. Dense x axis for fit curve
%     % Always plot in linear cm/s; evaluate Gaussian in log2 if isLog2
%     fineSpeedAxis_linear = logspace(log10(min(movingCenters)), log10(max(movingCenters)), 200);
%     if isLog2
%         fineSpeedAxis_fit = log2(fineSpeedAxis_linear);  % log2 domain — matches fitParams
%     else
%         % fineSpeedAxis_fit = fineSpeedAxis_linear;       % linear domain (not used)
%         fineSpeedAxis_fit = fineSpeedAxis_linear;
%     end
% 
%     % 4. Figure
%     figHandle = figure('Name', 'Speed Tuning Representatives', ...
%         'Position', [50, 150, 1500, 380], 'Color', 'w');
% 
%     colWidths = 0.22;
%     colGap    = 0.02;
%     leftStart = 0.03;
% 
%     for c = 1:4
%         r    = targetROIs(c);
%         xPos = leftStart + (c-1) * (colWidths + colGap);
% 
%         pPanel = uipanel(figHandle, 'Position', [xPos, 0.02, colWidths, 0.95], ...
%             'BackgroundColor', 'w', 'BorderType', 'none');
% 
%         if isempty(r) || isnan(r)
%             axEmpty = axes(pPanel, 'Position', [0 0 1 1]); 
%             text(axEmpty, 0.5, 0.5, 'Category Empty', 'HorizontalAlignment', 'center', 'FontSize', 12);
%             title(axEmpty, plotLabels{c}); box off; axis off;
%             continue;
%         end
% 
%         tloInner = tiledlayout(pPanel, 1, 10, 'TileSpacing', 'none', 'Padding', 'none');
% 
%         y_stat    = statMean(r);
%         y_statErr = statSEM(r);
%         y_move    = moveMean(r, :);
%         y_moveErr = moveSEM(r, :);
%         params    = fitParams(r, :);   % in log2 space if isLog2
%         currentR2 = r2Values(r);
%         currentP  = pValuesMoving(r);
% 
%         allY    = [y_stat+y_statErr, y_stat-y_statErr, y_move+y_moveErr, y_move-y_moveErr];
%         yLimits = [min(allY)*0.85, max(allY)*1.15];
% 
%         %% --- STATIONARY SLOT ---
%         axStat = nexttile(tloInner, 1, [1, 2]);
%         hold(axStat, 'on');
%         errorbar(axStat, 0, y_stat, y_statErr, 'ok', 'MarkerSize', 6.5, ...
%             'LineWidth', 1.2, 'MarkerFaceColor', 'none', 'DisplayName', 'Resting');
%         xlim(axStat, [-0.5, 0.5]); ylim(axStat, yLimits);
%         axStat.XTick = 0; axStat.XTickLabel = {'0'};
%         set(axStat, 'TickDir', 'out', 'LineWidth', 1.1, 'FontSize', 10); box(axStat, 'off');
%         if c > 1
%             ylabel(axStat, '');
%         else
%             ylabel(axStat, '\DeltaF/F (Neu)', 'FontSize', 12);
%         end
% 
%         %% --- LOCOMOTION DOMAIN ---
%         axMove = nexttile(tloInner, 3, [1, 8]);
%         hold(axMove, 'on');
%         errorbar(axMove, movingCenters, y_move, y_moveErr, 'ok', ...
%             'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerSize', 5.5, ...
%             'LineWidth', 1.1, 'DisplayName', 'Moving Data');
% 
%         % Evaluate fit in the correct domain, plot against linear x
%         if all(isfinite(params))
%             y_fit = gaussFun(params, fineSpeedAxis_fit);
%             % y_fit = gaussFun(params, fineSpeedAxis_linear); % ← linear version (not used)
% 
%             if c == 1 || c == 4
%                 lineColor = [0.85, 0.33, 0.1];
%             else
%                 lineColor = [0, 0.45, 0.74];
%             end
%             plot(axMove, fineSpeedAxis_linear, y_fit, 'Color', lineColor, ...
%                 'LineWidth', 2.5, 'DisplayName', 'Fit');
%         end
% 
%         set(axMove, 'XScale', 'log');
%         xlim(axMove, [min(movingCenters)*0.85, max(movingCenters)*1.15]);
%         ylim(axMove, yLimits);
%         axMove.XTick      = [2, 5, 10, 20, 40];
%         axMove.XTickLabel = {'2', '5', '10', '20', '40'};
%         set(axMove, 'TickDir', 'out', 'LineWidth', 1.1, 'FontSize', 10, 'YTickLabel', '');
%         axMove.YAxis.Visible = 'off'; box(axMove, 'off'); grid(axMove, 'on');
% 
%         titleString = sprintf('%s\\newline\\bf\\fontsize{10}ROI %d   \\rm(R^2 = %.3f | p = %.4f)', ...
%             plotLabels{c}, r, currentR2, currentP);
%         title(axMove, titleString, 'FontWeight', 'normal', 'Interpreter', 'tex', ...
%             'HorizontalAlignment', 'center');
%         xlabel(axMove, 'Speed (cm/s)', 'FontSize', 11);
% 
%         if c == 1
%             legend(axMove, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
%         end
%     end
% end

function figHandle = plotAllSpeedTuningCategories(response, useField)
    if nargin < 2, useField = 'dFFNeuropilCorrected'; end

    p_thresh = 0.01;

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

    allTypes      = classStruct.tuningType;
    r2Values      = classStruct.R2;
    fitParams     = classStruct.fitParams;

    if isfield(dataStruct, 'pValMoving')
        pValuesMoving = dataStruct.pValMoving;
    else
        warning('pValMoving not found. Bypassing shuffle check.');
        pValuesMoving = zeros(size(r2Values));
    end

    categories = {'lowpass', 'highpass', 'bandpass', 'trough_inverted'};
    plotLabels = {'Low-Pass', 'High-Pass', 'Band-Pass', 'Trough / Inverted'};

    targetROIs = nan(1, 4);
    for c = 1:4
        idxList = find(strcmp(allTypes, categories{c}) & (pValuesMoving <= p_thresh));
        if ~isempty(idxList)
            [~, bestMatch] = max(r2Values(idxList));
            targetROIs(c)  = idxList(bestMatch);
        else
            idxListFallback = find(strcmp(allTypes, categories{c}));
            if ~isempty(idxListFallback)
                [~, bestMatch] = max(r2Values(idxListFallback));
                targetROIs(c)  = idxListFallback(bestMatch);
                warning('Category %s: no cells pass p<=%.3f. Plotting unvalidated best fit.', categories{c}, p_thresh);
            end
        end
    end

    gaussFun = @(params, xdata) params(1) + params(2) .* ...
        exp(-(((xdata - params(3)).^2) / (2 * (params(4).^2))));

    fineSpeedAxis_linear = logspace(log10(min(movingCenters)), log10(max(movingCenters)), 200);
    if isLog2
        fineSpeedAxis_fit = log2(fineSpeedAxis_linear);
    else
        fineSpeedAxis_fit = fineSpeedAxis_linear;
    end

    % Figure sized in cm at final intended size
    figW_cm = 18;
    figH_cm = 5;
    figHandle = figure('Name', 'Speed Tuning Representatives', ...
        'Units', 'centimeters', 'Position', [2, 2, figW_cm, figH_cm], ...
        'Color', 'w', 'PaperPositionMode', 'auto');

    colWidth    = 0.20;
    colGap      = 0.02;
    leftStart   = 0.04;
    rowHeight   = 0.68;
    bottomStart = 0.20;

    for c = 1:4
        r    = targetROIs(c);
        xPos = leftStart + (c-1) * (colWidth + colGap);

        if isnan(r)
            axEmpty = axes(figHandle, 'Position', [xPos, bottomStart, colWidth, rowHeight]); 
            text(axEmpty, 0.5, 0.5, 'Category Empty', 'HorizontalAlignment', 'center', 'FontSize', 8);
            title(axEmpty, plotLabels{c}); box off; axis off;
            continue;
        end

        axStat = axes(figHandle, 'Position', [xPos,                   bottomStart, colWidth*0.25, rowHeight]);
        axMove = axes(figHandle, 'Position', [xPos + colWidth*0.28,   bottomStart, colWidth*0.72, rowHeight]);

        y_stat    = statMean(r);
        y_statErr = statSEM(r);
        y_move    = moveMean(r, :);
        y_moveErr = moveSEM(r, :);
        params    = fitParams(r, :);
        currentR2 = r2Values(r);
        currentP  = pValuesMoving(r);

        allY    = [y_stat+y_statErr, y_stat-y_statErr, y_move+y_moveErr, y_move-y_moveErr];
        yLimits = [min(allY)*0.85, max(allY)*1.15];

        % Stationary
        hold(axStat, 'on');
        errorbar(axStat, 0, y_stat, y_statErr, 'ok', 'MarkerSize', 2, ...
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
            'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerSize', 2, ...
            'LineWidth', 0.5, 'DisplayName', 'Moving Data');

        if all(isfinite(params))
            y_fit = gaussFun(params, fineSpeedAxis_fit);
            if c == 1 || c == 4
                lineColor = [0.85, 0.33, 0.1];
            else
                lineColor = [0, 0.45, 0.74];
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
            plotLabels{c}, r, currentR2, currentP);
        title(axMove, titleString, 'FontWeight', 'normal', 'Interpreter', 'tex', ...
            'FontSize', 7, 'HorizontalAlignment', 'center');
        xlabel(axMove, 'Speed (cm/s)', 'FontSize', 7);
    end

    exportgraphics(figHandle, 'speedTuning_categories.pdf', 'ContentType', 'vector');
end