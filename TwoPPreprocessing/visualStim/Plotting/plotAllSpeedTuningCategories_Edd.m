function figHandle = plotAllSpeedTuningCategories_Edd(response, useField, alphaShuffle, customOffset, axisMode)
% Plots the best-fitting example of each tuning shape, using Edd's
% fitGaussianTemplates_tuning (same classification as the population script).
%
%   useField     - e.g. 'dFFNeuropilCorrected'
%   alphaShuffle - selection threshold on pValMoving (default 0.01)
%   customOffset - Edd's offset (default 0.5, as for visual speed)
%   axisMode     - x-axis for the tuning curves (default 'linear'):
%                    'linear' : running speed in cm/s, linear axis
%                    'log'    : running speed in cm/s, log axis (as in Saleem et al., 2013)
%                    'rank'   : bin rank, labelled with bin-centre speeds
%                  This only affects the display. The templates are always
%                  fitted in bin-rank units and mapped to cm/s for plotting.

    if nargin < 2 || isempty(useField),     useField     = 'dFFNeuropilCorrected'; end
    if nargin < 3 || isempty(alphaShuffle), alphaShuffle = 0.01;     end
    if nargin < 4 || isempty(customOffset), customOffset = 0.5;      end
    if nargin < 5 || isempty(axisMode),     axisMode     = 'rank'; end

    tc = response.tuningCurve;
    E  = classifySessionEdd(tc, useField, alphaShuffle, customOffset);

    statMean = tc.(useField).statMean;
    statSEM  = tc.(useField).statSEM;
    moveMean = tc.(useField).moveMean;
    moveSEM  = tc.(useField).moveSEM;
    pVals    = tc.(useField).pValMoving;

    centres = E.centres;
    nB      = numel(centres);
    rDense  = linspace(1, nB, 200);                 % fit is in bin-rank units
    xDense  = interp1(1:nB, centres, rDense);       % mapped to cm/s for display
    xMax    = ceil(max(centres) / 10) * 10;

    plotLabels = {'Low-Pass', 'High-Pass', 'Band-Pass', 'Trough / Inverted'};   % Edd codes 1-4
    lineColors = {[0.85 0.33 0.10], [0 0.45 0.74], [0 0.45 0.74], [0.85 0.33 0.10]};

    % Best-fitting speed-tuned example of each shape
    targetK = nan(1, 4);
    for c = 1:4
        idx = find(E.code == c);
        if ~isempty(idx)
            [~, best] = max(E.R2(idx));
            targetK(c) = idx(best);
        end
    end

    figHandle = figure('Name', sprintf('Speed Tuning Representatives (Edd, %s axis)', axisMode), ...
        'Units', 'centimeters', 'Position', [2, 2, 18, 5], ...
        'Color', 'w', 'PaperPositionMode', 'auto');

    colWidth = 0.20; colGap = 0.02; leftStart = 0.04;
    rowHeight = 0.68; bottomStart = 0.20;

    for c = 1:4
        xPos = leftStart + (c - 1) * (colWidth + colGap);
        k    = targetK(c);

        if isnan(k)
            axEmpty = axes(figHandle, 'Position', [xPos, bottomStart, colWidth, rowHeight]);
            text(axEmpty, 0.5, 0.5, 'Category empty', 'HorizontalAlignment', 'center', 'FontSize', 8);
            title(axEmpty, plotLabels{c}); axis(axEmpty, 'off');
            continue;
        end

        r = E.sigIdx(k);
        p = E.params(k, :);

        axStat = axes(figHandle, 'Position', [xPos,                 bottomStart, colWidth * 0.25, rowHeight]);
        axMove = axes(figHandle, 'Position', [xPos + colWidth*0.28, bottomStart, colWidth * 0.72, rowHeight]);

        y_move = moveMean(r, :); y_moveErr = moveSEM(r, :);
        allY   = [statMean(r) + statSEM(r), statMean(r) - statSEM(r), y_move + y_moveErr, y_move - y_moveErr];
        yLim   = [min(allY), max(allY)] + [-0.1 0.1] * range(allY);

        % Stationary (speed <= 1 cm/s)
        hold(axStat, 'on');
        errorbar(axStat, 0, statMean(r), statSEM(r), 'ok', 'MarkerSize', 2, 'LineWidth', 0.5);
        xlim(axStat, [-0.5 0.5]); ylim(axStat, yLim);
        axStat.XTick = 0; axStat.XTickLabel = {'0'};
        set(axStat, 'TickDir', 'out', 'LineWidth', 0.5, 'FontSize', 7); box(axStat, 'off');
        if c == 1, ylabel(axStat, '\DeltaF/F', 'FontSize', 8); end

        % Moving, with Edd's fitted template
        hold(axMove, 'on');
        yFit = p(1) + p(2) * exp(-((rDense - p(3)).^2) / (2 * p(4)^2));

        switch axisMode
            case 'linear'   % running speed (cm/s), linear axis
                errorbar(axMove, centres, y_move, y_moveErr, 'ok', ...
                    'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerSize', 2, 'LineWidth', 0.5);
                plot(axMove, xDense, yFit, 'Color', lineColors{c}, 'LineWidth', 1.5);
                xlim(axMove, [0, xMax]);
                axMove.XTick = 0:10:xMax;

            case 'log'      % running speed (cm/s), log axis, as in Saleem et al. (2013)
                errorbar(axMove, centres, y_move, y_moveErr, 'ok', ...
                    'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerSize', 2, 'LineWidth', 0.5);
                plot(axMove, xDense, yFit, 'Color', lineColors{c}, 'LineWidth', 1.5);
                set(axMove, 'XScale', 'log');
                axMove.XTick = [2 5 10 20 40]; axMove.XTickLabel = {'2', '5', '10', '20', '40'};
                xlim(axMove, [min(centres) * 0.85, max(centres) * 1.15]);

            case 'rank'     % bin rank, labelled with bin-centre speeds
                errorbar(axMove, 1:nB, y_move, y_moveErr, 'ok', ...
                    'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerSize', 2, 'LineWidth', 0.5);
                plot(axMove, rDense, yFit, 'Color', lineColors{c}, 'LineWidth', 1.5);
                xlim(axMove, [0.5, nB + 0.5]);
                tickIdx = unique(round(linspace(1, nB, 4)));
                axMove.XTick      = tickIdx;
                axMove.XTickLabel = arrayfun(@(v) sprintf('%.0f', v), centres(tickIdx), 'UniformOutput', false);

            otherwise
                error('axisMode must be ''linear'', ''log'' or ''rank''.');
        end

        ylim(axMove, yLim);
        set(axMove, 'TickDir', 'out', 'LineWidth', 0.5, 'FontSize', 7, 'YTickLabel', '');
        axMove.YAxis.Visible = 'off'; box(axMove, 'off');

        % p-values from 1,000 shuffles cannot resolve below 0.001
        if pVals(r) < 0.001
            pStr = 'p < 0.001';
        else
            pStr = sprintf('p = %.3f', pVals(r));
        end
        title(axMove, sprintf('%s\\newline\\fontsize{7}ROI %d   (R^2 = %.3f | %s)', ...
            plotLabels{c}, r, E.R2(k), pStr), ...
            'FontWeight', 'normal', 'Interpreter', 'tex', 'FontSize', 7);
        xlabel(axMove, 'Speed (cm/s)', 'FontSize', 7);
    end
end

%% ---------------------------------------------------------------------
function E = classifySessionEdd(tc, useField, alphaShuffle, customOffset)
% Selects speed-tuned boutons (circular-shift test on moving frames) and
% classifies their tuning shape with Edd's fitGaussianTemplates_tuning.
    Y       = tc.(useField).moveMean;
    pvals   = tc.(useField).pValMoving(:);
    edges   = tc.speedBins(:)';
    E.centres = edges(1:end-1) + diff(edges) / 2;
    E.nROIs   = size(Y, 1);
    E.sigIdx  = find(pvals <= alphaShuffle);

    n  = numel(E.sigIdx);
    nB = numel(E.centres);
    E.code = nan(n, 1); E.pref = nan(n, 1); E.R2 = nan(n, 1); E.params = nan(n, 4);

    for k = 1:n
        y = Y(E.sigIdx(k), :);
        if any(isnan(y)), continue; end
        [p, c, r2] = fitGaussianTemplates_tuning(num2cell(y(:)), customOffset, false);
        E.code(k) = c; E.R2(k) = r2; E.params(k, :) = p;
        if c == 3 || c == 4
            E.pref(k) = interp1(1:nB, E.centres, min(max(p(3), 1), nB));
        end
    end
end