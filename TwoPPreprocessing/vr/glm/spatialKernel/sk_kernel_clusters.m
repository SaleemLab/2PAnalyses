function [lab, Cn] = sk_kernel_clusters(source, k)
% Recurring kernel shapes in one kernel set (significant cells). Kernels (G beyond 10 cm,
% z-scored per cell; the onset transient is measured separately) are reduced by PCA and grouped by k-means on the top 10 components. Reports:
%   - variance explained by the first components
%   - each cluster's mean kernel, with the mean measured map, GLM1 vision + speed
%     prediction and speed profile of its cells
%   - how clusters are spread over sessions (a shape found in one session only is suspect)
%   - stability: clustering two random halves of the cells separately (and, with more than
%     one mouse, each mouse group separately) and matching centroids by correlation
% Returns cluster labels for the significant cells (in ks order) and the centroids.
% Writes figures/clusters_<source>.png, clusters_<source>.txt and the sorted image.
if nargin < 2, k = 6; end
C = sk_config(); rng(1); ks = sk_kernelset(source); m = find(ks.sig); x = ks.x; tex = [40 80 120 160];
zs = @(M) (M - mean(M, 2, 'omitnan')) ./ std(M, 0, 2, 'omitnan');
sel = x > 10; x = x(sel);                  % onset bins left out: the onset is its own feature
Z = zs(ks.G(m, sel)); nPC = 10;
[coef, sc, ~, ~, ve] = pca(Z);
km = @(S) kmeans(S, k, 'Replicates', 20, 'MaxIter', 500);
lab = km(sc(:, 1:nPC));
Cn = cell2mat(arrayfun(@(j) mean(Z(lab == j,:), 1), (1:k)', 'UniformOutput', false));
% order clusters by size
[~, o] = sort(accumarray(lab, 1), 'descend'); r(o) = 1:k; lab = r(lab)'; Cn = Cn(o,:);

% stability
cen = @(idx) centroids(Z(idx,:), coef(:, 1:nPC), k);
match = @(A, B) matchCorr(A, B);
h = randperm(numel(m)) <= numel(m)/2; rHalf = match(cen(h), cen(~h));
mice = ks.cells.mouseID(m); um = unique(mice);
if numel(um) > 1
    g = ismember(mice, um(1:2:end)); rMouse = match(cen(g), cen(~g));
else, rMouse = nan(1, k); end
% how well do the full-data centroids match each half (cluster by cluster)
rFullHalf = match(Cn, cen(h));

fid = fopen(fullfile(C.outRoot, sprintf('clusters_%s.txt', source)), 'w');
out = @(varargin) cellfun(@(f) fprintf(f, varargin{:}), {1, fid}, 'UniformOutput', false);
out('%s: %d significant cells, k = %d, PCA variance explained by PC1-5: %s%%\n', source, numel(m), k, mat2str(round(ve(1:5)', 1)));
out('Matched centroid correlation, random halves: %s\n', mat2str(round(rHalf, 2)));
out('Matched centroid correlation, mouse groups (%s vs rest): %s\n', strjoin(um(1:2:end), '+'), mat2str(round(rMouse, 2)));
sess = ks.cells.sessionid(m); us = unique(sess);
out('%-8s %5s %8s  %s\n', 'cluster', 'n', 'mainPk', ['fraction of each session''s cells: ' sprintf('%d ', us)]);
for j = 1:k
    [~, p] = max(Cn(j,:)); fr = arrayfun(@(s) mean(lab(sess == s) == j), us);
    out('%-8d %5d %8d  %s\n', j, nnz(lab == j), x(p), sprintf('%.2f ', fr));
end
fclose(fid);

% figure
nc = ceil(k/2) + 1;
fig = figure('Position', [20 20 1800 760], 'Color', 'w'); tl = tiledlayout(2, nc, 'TileSpacing', 'compact');
V = zs(ks.comp.mapVS(m, sel)); O = zs(ks.comp.mapObs(m, sel)); S = zs(ks.comp.speedMap(m, sel));
for j = 1:k
    nexttile; hold on; s = lab == j;
    for L = tex, patch([L-3 L+3 L+3 L-3], [-1.5 -1.5 4 4], [.95 .88 .5], 'EdgeColor', 'none', 'FaceAlpha', .6); end
    q = prctile(Z(s,:), [25 75]); fill([x fliplr(x)], [q(1,:) fliplr(q(2,:))], [.85 .85 .85], 'EdgeColor', 'none');
    plot(x, mean(O(s,:)), 'Color', [.2 .45 .8]); plot(x, mean(V(s,:)), 'Color', [.85 .25 .15]); plot(x, mean(S(s,:)), 'Color', [.3 .6 .3]);
    plot(x, Cn(j,:), 'k', 'LineWidth', 1.6); xlim([10 200]); ylim([-1.5 4]); box off;
    title(sprintf('cluster %d: %d cells (half-split r = %.2f)', j, nnz(s), rFullHalf(j)));
    if j == 1, legend({'','','','','IQR','measured map','GLM1 vision+speed','speed','kernel'}, 'Box', 'off', 'Location', 'northeast'); end
end
nexttile; bar(ve(1:10)); xlabel('PC'); ylabel('% variance'); box off; title('PCA of z-scored kernels');
nexttile; hold on; plot(x, coef(:, 1:3), 'LineWidth', 1.2); for L = tex, xline(L, ':'); end
legend({'PC1','PC2','PC3'}, 'Box', 'off'); xlabel('position (cm)'); box off; title('First components');
title(tl, sprintf('%s: recurring kernel shapes (z-scored, k-means on %d PCs)', strrep(source, '_', ' '), nPC));
exportgraphics(fig, fullfile(C.figRoot, sprintf('clusters_%s.png', source)), 'Resolution', 120);
sk_plot_kernelset(source, lab);
end

function Cn = centroids(Zs, coef, k)
lab = kmeans(Zs * coef, k, 'Replicates', 20, 'MaxIter', 500);
Cn = cell2mat(arrayfun(@(j) mean(Zs(lab == j,:), 1), (1:k)', 'UniformOutput', false));
end

function r = matchCorr(A, B)
% best one-to-one matching of centroids (rows) by correlation; r(i) for each row of A
R = corr(A', B'); M = matchpairs(-R, 10); r = nan(1, size(A, 1));
for i = 1:size(M, 1), r(M(i,1)) = R(M(i,1), M(i,2)); end
end
