function plotRFTintedGrid(tl, boutonData, uAz, uEl_plot, timeVector)
% plotRFTintedGrid: Small-multiples grid (one subplot per Az/El position,
% matching physical screen layout). Each tile's background is tinted by
% that position's peak response amplitude (relative to the bouton's own
% max), so spatial selectivity is visible at a glance while traces stay
% fully legible (no overlap, unlike plotRFHeatmapWithTraces).
%
% INPUTS:
%   tl          - target tiledlayout handle (nEl x nAz), e.g.
%                   tl = tiledlayout(figure, nEl, nAz, 'TileSpacing', 'tight');
%   boutonData  - single element of an RFMapping struct array (from
%                 analyseRFMapping.m), must contain:
%                   .meanGridResponse (nEl x nAz)
%                   .meanTemporalResponse (nTpts x nEl x nAz)
%                   .meanBlankResponse (nTpts x 1)
%                   .centerAz, .centerEl, .isResponsive
%   uAz         - vector of unique azimuth positions
%   uEl_plot    - vector of unique elevation positions
%   timeVector  - shared time vector (nTpts x 1)
%
% USAGE:
%   fig = figure;
%   tl = tiledlayout(fig, nEl, nAz, 'TileSpacing', 'tight', 'Padding', 'compact');
%   plotRFTintedGrid(tl, allRFMapping(iBouton), uAz, uEl_plot, timeVector);

    nAz = length(uAz);
    nEl = length(uEl_plot);

    tempStack  = boutonData.meanTemporalResponse;   % nTpts x nEl x nAz
    blankTrace = boutonData.meanBlankResponse;
    grid       = boutonData.meanGridResponse;
    gridMax    = max(grid(:), [], 'omitnan');
    if gridMax <= 0 || isnan(gridMax), gridMax = eps; end

    yMax = max(tempStack(:), [], 'omitnan');
    yMin = min([tempStack(:); blankTrace(:)], [], 'omitnan');
    if isnan(yMax) || yMax == yMin, yMax = 1; yMin = -1; end

    cmap = bone(256);  % light-to-dark colormap for tile backgrounds

    for r = 1:nEl
        for c = 1:nAz
            ax = nexttile(tl, (r-1)*nAz + c);

            % tile background color from this position's amplitude
            normAmp   = max(0, min(1, grid(r, c) / gridMax));
            tileColor = cmap(round(normAmp * 255) + 1, :);
            set(ax, 'Color', tileColor);
            hold(ax, 'on');

            plot(ax, timeVector, squeeze(tempStack(:, r, c)), ...
                'Color', [0.85 0.15 0.1], 'LineWidth', 1.1);
            plot(ax, timeVector, blankTrace, ...
                'Color', [1 1 1], 'LineWidth', 0.8, 'LineStyle', ':');

            if boutonData.isResponsive && r == find(uEl_plot == boutonData.centerEl, 1) && ...
                    c == find(uAz == boutonData.centerAz, 1)
                set(ax, 'XColor', [1 0.85 0], 'YColor', [1 0.85 0], 'LineWidth', 2.5, 'Box', 'on');
            end

            ylim(ax, [yMin yMax]);
            xlim(ax, [timeVector(1) timeVector(end)]);
            set(ax, 'XTick', [], 'YTick', []);
        end
    end
end
