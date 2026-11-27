function sk_plot_cell_cards(source, cellIdx, name)
% One card per cell, on a shared track axis: the measured map (odd and even base trials),
% GLM1's full and vision + speed predictions, the position kernel G (with +-1 SE where the
% fit gives SEs), and the session's mean running speed. Up to 16 cells per page.
% CELLIDX: CellIdx values (V1_GLMData), or a number N to take the N significant cells with
% the deepest kernels beyond 10 cm (the onset transient would otherwise dominate). Writes figures/cards_<source>_<name>_<page>.png.
if nargin < 2, cellIdx = 16; end
if nargin < 3, name = 'top'; end
C = sk_config(); ks = sk_kernelset(source); x = ks.x; tex = [40 80 120 160];
if isscalar(cellIdx)
    d = range(ks.G(:, x > 10), 2); d(~ks.sig) = -inf;         % depth beyond the onset bins
    [~, o] = sort(d, 'descend'); rows = o(1:min(cellIdx, nnz(ks.sig)));
else
    [tf, rows] = ismember(cellIdx, ks.cells.CellIdx); rows = rows(tf);
end
seSc = 1 ./ mean(ks.comp.mapFull, 2, 'omitnan');             % SE on the G scale (identity link)
for pg = 1:ceil(numel(rows) / 16)
    rr = rows((pg-1)*16 + 1:min(pg*16, numel(rows)));
    fig = figure('Position', [20 20 1800 1000], 'Color', 'w'); tl = tiledlayout(4, 4, 'TileSpacing', 'compact');
    for r = rr(:)'
        t = tiledlayout(tl, 2, 1, 'TileSpacing', 'none'); t.Layout.Tile = find(rr == r);
        ax = nexttile(t); hold(ax, 'on');
        for L = tex, patch(ax, [L-3 L+3 L+3 L-3], [0 0 1 1] * max(ks.comp.mapObs(r,:)) * 1.3, [.95 .88 .5], 'EdgeColor', 'none', 'FaceAlpha', .6); end
        plot(ax, x, ks.comp.mapObsA(r,:), 'Color', [.6 .6 .6]); plot(ax, x, ks.comp.mapObsB(r,:), 'Color', [.6 .6 .6]);
        plot(ax, x, ks.comp.mapObs(r,:), 'k', 'LineWidth', 1.3);
        plot(ax, x, ks.comp.mapVS(r,:), 'Color', [.85 .25 .15], 'LineWidth', 1.1);
        plot(ax, x, ks.comp.mapFull(r,:), '--', 'Color', [.2 .45 .8], 'LineWidth', 1.1);
        set(ax, 'XTick', [], 'XLim', [0 200]); axis(ax, 'tight'); xlim(ax, [0 200]); box(ax, 'off');
        sg = ''; if ks.sig(r), sg = '*'; end
        title(ax, sprintf('%s%s  (map r_{odd,even} %.2f)', strrep(ks.cells.CellName{r}, '_', ' '), sg, ks.comp.relMap(r)), 'FontSize', 8, 'FontWeight', 'normal');
        ax = nexttile(t); hold(ax, 'on');
        g = ks.G(r,:); se = ks.SE(r,:) * seSc(r);
        if any(~isnan(se)), fill(ax, [x fliplr(x)], [g - se fliplr(g + se)], [.8 .8 .8], 'EdgeColor', 'none'); end
        plot(ax, x, g, 'k', 'LineWidth', 1.4); yline(ax, 0, ':');
        for L = tex, xline(ax, L, ':', 'Color', [.8 .6 .1]); end
        yyaxis(ax, 'right'); plot(ax, x, ks.comp.speedMap(r,:), 'Color', [.3 .6 .3 .6]); ax.YAxis(2).Color = [.3 .6 .3];
        yyaxis(ax, 'left'); xlim(ax, [0 200]); box(ax, 'off');
    end
    title(tl, sprintf(['%s. Top: measured map (black; grey odd/even trials), GLM1 vision + speed (red), GLM1 full (blue dashed). ' ...
        'Bottom: position kernel as fractional change (black), mean speed (green). * position significant'], strrep(source, '_', ' ')), 'FontSize', 9);
    exportgraphics(fig, fullfile(C.figRoot, sprintf('cards_%s_%s_%d.png', source, name, pg)), 'Resolution', 110);
    close(fig);
end
end
