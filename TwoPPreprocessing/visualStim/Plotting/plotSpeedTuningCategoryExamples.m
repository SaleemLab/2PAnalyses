function fig = plotSpeedTuningCategoryExamples(response, targetStruct, useField, categoryName, plotOpts)
% PLOTSPEEDTUNINGCATEGORYEXAMPLES
%
% Example:
%   plotOpts = struct();
%   plotOpts.maxPlots = 12;
%   plotOpts.onlySignificantMoving = true;
%   plotOpts.sortBy = 'R2';
%   plotSpeedTuningCategoryExamples(response, 'tuningCurve', 'dFFNeuropilCorrected', 'bandpass', plotOpts)

    if nargin < 5 || isempty(plotOpts)
        plotOpts = struct();
    end
    if ~isfield(plotOpts, 'maxPlots'),              plotOpts.maxPlots = 12; end
    if ~isfield(plotOpts, 'nCols'),                 plotOpts.nCols = 4; end
    if ~isfield(plotOpts, 'onlySignificantMoving'), plotOpts.onlySignificantMoving = true; end
    if ~isfield(plotOpts, 'sortBy'),                plotOpts.sortBy = 'R2'; end
    if ~isfield(plotOpts, 'showSEM'),               plotOpts.showSEM = true; end
    if ~isfield(plotOpts, 'roiIdx'),                plotOpts.roiIdx = []; end
    if ~isfield(plotOpts, 'figureName'),            plotOpts.figureName = ''; end
    if ~isfield(plotOpts, 'classficationType'),    plotOpts.classficationType = 'classification'; end 

    dataStruct = response.(targetStruct).(useField);

    % Load log2 classification if available, otherwise fall back to standard
%     if isfield(dataStruct, 'classification') && contains(plotOpts.classficationType, 'classification')
%         cls = dataStruct.classification;
%         isLog2 = true;
%     else
    cls = dataStruct.classification;
    isLog2 = false;
%     end

    edges = response.(targetStruct).speedBins;
    x     = (edges(1:end-1) + diff(edges)/2)';   % linear cm/s for plotting

    yMean    = dataStruct.moveMean;
    ySEM     = dataStruct.moveSEM;
    statMean = dataStruct.statMean;
    statSEM  = dataStruct.statSEM;

    % Gaussian function — params always in the fitted domain
    gaussFun = @(params, xdata) ...
        params(1) + params(2) .* exp(-(((xdata - params(3)).^2) ./ (2 * (params(4).^2))));

    idx = find(strcmp(cls.tuningType, categoryName));
    if plotOpts.onlySignificantMoving && isfield(dataStruct, 'isSignificantMoving_999')
        idx = idx(dataStruct.isSignificantMoving_999(idx));
    end
    if ~isempty(plotOpts.roiIdx)
        idx = intersect(idx, plotOpts.roiIdx(:)');
    end
    if isempty(idx)
        warning('No ROIs found for category: %s', categoryName);
        fig = [];
        return;
    end

    switch lower(plotOpts.sortBy)
        case 'r2'
            [~, ord] = sort(cls.R2(idx), 'descend');
        case 'preferredspeed'
            [~, ord] = sort(cls.preferredSpeed(idx), 'ascend');
        otherwise
            ord = 1:numel(idx);
    end
    idx    = idx(ord);
    idx    = idx(1:min(plotOpts.maxPlots, numel(idx)));
    nPlots = numel(idx);
    nCols  = plotOpts.nCols;
    nRows  = ceil(nPlots / nCols);

    if isempty(plotOpts.figureName)
        figName = sprintf('%s | %s | %s', targetStruct, useField, categoryName);
    else
        figName = plotOpts.figureName;
    end

    fig = figure('Color', 'w', 'Name', figName, 'Position', [50 50 1500 900]);

    %
    % Always plot in linear cm/s on screen.
    % For log2 fits: evaluate Gaussian in log2 space, then plot against
    % the corresponding linear speeds.
    xDense_linear = logspace(log10(min(x)), log10(max(x)), 200);   % linear cm/s
    if isLog2
        xDense_fit = log2(xDense_linear);   % log2 domain — matches fitParams
    else
        % xDense_fit = xDense_linear;        % linear domain (commented out — using log2 classifier)
        xDense_fit = xDense_linear;
    end

    for i = 1:nPlots
        roi = idx(i);

        y_stat    = statMean(roi);
        y_statErr = statSEM(roi);
        y_move    = yMean(roi, :);
        y_moveErr = ySEM(roi, :);
        p         = cls.fitParams(roi, :);   % in log2 space if isLog2

        allY    = [y_stat + y_statErr, y_stat - y_statErr, y_move + y_moveErr, y_move - y_moveErr];
        yLimits = [min(allY)*0.85, max(allY)*1.15];
        if yLimits(1) == yLimits(2), yLimits = [yLimits(1)-0.1, yLimits(1)+0.1]; end

        colIdx = mod(i-1, nCols);
        rowIdx = floor((i-1) / nCols);

        panelLeft   = colIdx / nCols;
        panelBottom = 1 - ((rowIdx + 1) / nRows);
        panelWidth  = 1 / nCols;
        panelHeight = 1 / nRows;

        subPanel = uipanel(fig, 'Units', 'normalized', 'Position', ...
            [panelLeft, panelBottom, panelWidth, panelHeight], ...
            'BackgroundColor', 'w', 'BorderType', 'none');

        tloInner = tiledlayout(subPanel, 1, 10, 'TileSpacing', 'none', 'Padding', 'none');

        %% --- PANEL 1: STATIONARY ---
        axStat = nexttile(tloInner, 1, [1, 2]);
        hold(axStat, 'on');
        if plotOpts.showSEM && ~isnan(y_statErr)
            errorbar(axStat, 0, y_stat, y_statErr, 'ok', 'MarkerSize', 5, ...
                'LineWidth', 1, 'MarkerFaceColor', 'none');
        else
            plot(axStat, 0, y_stat, 'ok', 'MarkerSize', 5);
        end
        xlim(axStat, [-0.5, 0.5]); ylim(axStat, yLimits);
        axStat.XTick = 0; axStat.XTickLabel = {'0'};
        set(axStat, 'TickDir', 'out', 'LineWidth', 1.0, 'FontSize', 8); box(axStat, 'off');
        if mod(i-1, nCols) == 0
            ylabel(axStat, 'Response');
        end

        %% --- PANEL 2: LOCOMOTION ---
        axMove = nexttile(tloInner, 3, [1, 8]);
        hold(axMove, 'on');

        valid = ~isnan(y_move);
        xv    = x(valid);
        yv    = y_move(valid);
        ev    = y_moveErr(valid);

        if plotOpts.showSEM && ~isempty(ev)
            errorbar(axMove, xv, yv, ev, 'ok', ...
                'MarkerFaceColor', [0.15 0.15 0.15], ...
                'MarkerSize', 4, 'LineWidth', 1);
        else
            plot(axMove, xv, yv, 'ok', ...
                'MarkerFaceColor', [0.15 0.15 0.15], ...
                'MarkerSize', 4);
        end

        % Plot fitted curve
        % For log2: evaluate Gaussian in log2 space, display against linear x
        % For linear: evaluate Gaussian in linear space (commented out)
        if all(isfinite(p))
            yFit = gaussFun(p, xDense_fit);   % xDense_fit is log2 if isLog2
            % yFit = gaussFun(p, xDense_linear); % ← linear version (not used)

            if strcmp(categoryName, 'lowpass') || strcmp(categoryName, 'trough_inverted')
                lineColor = [0.85, 0.33, 0.1];
            else
                lineColor = [0, 0.45, 0.74];
            end
            plot(axMove, xDense_linear, yFit, '-', 'Color', lineColor, 'LineWidth', 2);
        end

        % X axis always in linear cm/s with log scale display
        set(axMove, 'XScale', 'log');
        xlim(axMove, [min(x)*0.85, max(x)*1.15]);
        ylim(axMove, yLimits);
        axMove.XTick      = [2, 5, 10, 20, 40];
        axMove.XTickLabel = {'2', '5', '10', '20', '40'};
        set(axMove, 'TickDir', 'out', 'LineWidth', 1.0, 'FontSize', 8, 'YTickLabel', '');
        axMove.YAxis.Visible = 'off'; box(axMove, 'off'); grid(axMove, 'on');

        titleString = sprintf('ROI %d | R2=%.2f | pref=%.1f', ...
            roi, cls.R2(roi), cls.preferredSpeed(roi));
        title(axMove, titleString, 'FontSize', 8, 'FontWeight', 'normal');

        if ceil(i / nCols) == nRows
            xlabel(axMove, 'Speed (cm s^{-1})');
        end
    end

    sgtitle(fig, sprintf('%s examples (%d shown)', categoryName, nPlots), 'FontWeight', 'bold');
end
