function sk_plot_kernelset(source, order)
% Population view of one kernel set (significant cells): the kernel next to the measured
% map, GLM1's vision + speed prediction and GLM1's position prediction, all z-scored per
% cell and in the same row order. ORDER: 'peak' (main kernel peak, default), 'session',
% or a vector of cluster labels (one per significant cell, e.g. from sk_kernel_clusters).
% Texture positions are marked. Writes figures/kernelset_<source>.png.
if nargin < 2, order = 'peak'; end
C = sk_config(); ks = sk_kernelset(source); m = ks.sig; x = ks.x; tex = [40 80 120 160];
zs = @(M) (M - mean(M, 2, 'omitnan')) ./ std(M, 0, 2, 'omitnan');
Z = zs(ks.G(m,:)); [~, mp] = max(Z, [], 2); sess = ks.cells.sessionid(m);
if ischar(order) && strcmp(order, 'session'), [~, o] = sortrows([sess mp]); lab = sess; tag = 'session';
elseif ischar(order), [~, o] = sort(mp); lab = []; tag = 'peak';
else, [~, o] = sortrows([order(:) mp]); lab = order(:); tag = 'cluster'; end
panels = {Z, 'position kernel'; zs(ks.comp.mapObs(m,:)), 'measured map (base trials)'; ...
    zs(ks.comp.mapVS(m,:)), 'GLM1 vision + speed prediction'; zs(ks.comp.mapPos(m,:)), 'GLM1 position prediction'};
fig = figure('Position', [30 30 1700 900], 'Color', 'w'); tl = tiledlayout(1, 4, 'TileSpacing', 'compact');
for p = 1:4
    nexttile; imagesc(x, 1:nnz(m), panels{p,1}(o,:), [-2 3]); hold on;
    for L = tex, xline(L, 'w:', 'LineWidth', 1); end
    if ~isempty(lab), b = find(diff(lab(o)) ~= 0); yline(b + 0.5, 'w-', 'LineWidth', 1); end
    xlabel('position (cm)'); title(panels{p,2}); if p == 1, ylabel(sprintf('cells, sorted by %s', tag)); else, set(gca, 'YTick', []); end
end
colormap(parula); cb = colorbar; cb.Label.String = 'z-score within cell';
title(tl, sprintf('%s: %d cells with significant position (dotted: textures)', strrep(source, '_', ' '), nnz(m)));
exportgraphics(fig, fullfile(C.figRoot, sprintf('kernelset_%s_%s.png', source, tag)), 'Resolution', 120);
end
