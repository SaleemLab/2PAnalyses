function plotRFTraceGrid(fig, panelPos, boutonData, uAz, uEl_plot, timeVector)
% plotRFTraceGrid: Draws a small-multiples grid of raw temporal traces, one
% subplot per Az/El position, spatially laid out to match the orientation of
% the corresponding heatmap (plotRFHeatmapOnly). Meant to be placed alongside
% that function so you get one clean spatial map + one panel with the
% underlying temporal traces, instead of traces overlaid on the heatmap.
%
% All subplots share a common y-axis (based on this bouton's own min/max
% across all positions + blank) so relative amplitude across positions is
% visually honest. The subplot matching boutonData.centerAz/centerEl is
% outlined in yellow.
%
% INPUTS:
%   fig         - parent figure handle
%   panelPos    - [x y w h] normalized figure position for the whole grid
%   boutonData  - single element of an RFMapping struct array, must contain:
%                   .meanTemporalResponse (nTpts x nEl x nAz)
%                   .meanBlankResponse (nTpts x 1)
%                   .centerAz, .centerEl, .isResponsive
%   uAz         - vector of unique azimuth positions
%   uEl_plot    - vector of unique elevation positions
%   timeVector  - shared time vector (nTpts x 1)
%
% USAGE:
%   figA = figure('Color','w','Position',[50 50 1800 900]);
%   axHeat = axes('Position', [0.05 0.55 0.2 0.35]);
%   plotRFHeatmapOnly(axHeat, allRFMapping(iBouton), uAz, uEl_plot);
%   plotRFTraceGrid(figA, [0.05 0.05 0.2 0.42], allRFMapping(iBouton), ...
%       uAz, uEl_plot, timeVector);

    nAz = length(uAz);
    nEl = length(uEl_plot);

    tempStack  = boutonData.meanTemporalResponse; % nTpts x nEl x nAz
    blankTrace = boutonData.meanBlankResponse;

    yMax = max(tempStack(:), [], 'omitnan');
    yMin = min([tempStack(:); blankTrace(:)], [], 'omitnan');
    if isnan(yMin) || isnan(yMax) || yMin == yMax
        yMin = -1; yMax = 1; % fallback so ylim never errors on degenerate data
    end

    cellW = panelPos(3) / nAz;
    cellH = panelPos(4) / nEl;

    for r = 1:nEl
        for c = 1:nAz
            % row 1 corresponds to uEl_plot(1); flip so higher elevation
            % values sit higher on the figure, matching imagesc/YDir normal
            xPos = panelPos(1) + (c-1)*cellW;
            yPos = panelPos(2) + (nEl - r)*cellH;

            axSub = axes('Parent', fig, 'Position', ...
                [xPos yPos cellW*0.95 cellH*0.95]);
            hold(axSub, 'on');

            plot(axSub, timeVector, blankTrace, 'Color', [0.6 0.6 0.6], ...
                'LineWidth', 0.6, 'LineStyle', ':');
            plot(axSub, timeVector, squeeze(tempStack(:, r, c)), ...
                'Color', [0.85 0.2 0.15], 'LineWidth', 1);

            ylim(axSub, [yMin yMax]);
            xlim(axSub, [min(timeVector) max(timeVector)]);
            axis(axSub, 'off');

            % highlight the cell matching this bouton's peak/center position
            if boutonData.isResponsive && ...
               abs(uAz(c) - boutonData.centerAz) < 1e-6 && ...
               abs(uEl_plot(r) - boutonData.centerEl) < 1e-6
                set(axSub, 'Box', 'on', 'Visible', 'on', 'XTick', [], 'YTick', [], ...
                    'XColor', [1 0.85 0], 'YColor', [1 0.85 0], 'LineWidth', 2);
            end
        end
    end
end
