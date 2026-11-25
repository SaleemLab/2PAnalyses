function S = sk_kernel_properties()
% Do the earlier kernel observations (coverage, multi-peaks, co-location with
% the vision+speed peaks, ~25 cm peak spacing, onset edge peaks) survive in the
% refitted kernels, in particular the Poisson + vision x speed (VSIP) kernels?
%
% Kernel sets compared (spatial cells):
%   orig   original GLM1 position kernels (Gaussian, lasso) from V1_GLMData.mat
%   gauss  Gaussian ridge refit, VSP
%   pois   Poisson + vision x speed refit, VSIP (log multiplicative gain)
%   poisSig  the pois kernels of cells where position still significantly improves
%            held-out likelihood (p_pois_posI < 0.05, paired t over the 10 folds)
% Writes kernel_properties_summary.txt and a figure.

C = sk_config();
K = load(fullfile(C.outRoot, 'refit_kernels.mat'));
R = readtable(fullfile(C.outRoot, 'refit_spatialcells.csv'));
assert(isequal(R.CellIdx, K.CellIdx), 'CSV and kernel file are not in the same order');
f = 'D:\DATA\Sonali\GLM\ForAS\V1_GLMData.mat';
Korig = h5read(f, '/Kernels/PosKernel');            % [allCells x 100], row = CellIdx
VSm   = h5read(f, '/Maps/VS/mean');                 % [allCells x 6 x 100]
T = readtable(C.tablePath);
assert(isequal(T.CellIdx, (1:height(T))'), 'table rows are not in CellIdx order');
x = C.binCentres;

sets.orig  = Korig(K.CellIdx, :);
sets.gauss = K.gauss_VSP;
sets.pois  = K.pois_VSIP;
sig = R.p_pois_posI < 0.05 & R.LLHrelI_pois > 1;
sets.poisSig = K.pois_VSIP(sig, :);
VSbase = squeeze(VSm(K.CellIdx, 1, :));
vsOf.orig = VSbase; vsOf.gauss = VSbase; vsOf.pois = VSbase; vsOf.poisSig = VSbase(sig,:);
names = fieldnames(sets);

rng(3);
fid = fopen(fullfile(C.outRoot, 'kernel_properties_summary.txt'), 'w');
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('Spatial cells refitted: %d; position significant under Poisson + VxS: %d (%.0f%%)\n\n', ...
    height(R), nnz(sig), 100*mean(sig));
out('%-8s %6s %10s %10s %10s %12s %12s %10s %10s\n', 'set', 'n', 'multi@.3', 'multi@.5', 'multi@.7', ...
    'peak@0-2cm', 'coloc obs', 'coloc null', '25cm obs/null');
for i = 1:numel(names)
    Kn = sets.(names{i}); V = vsOf.(names{i});
    S.(names{i}) = kernelStats(Kn, V, x);
    s = S.(names{i});
    out('%-8s %6d %10.2f %10.2f %10.2f %12.2f %12.2f %10.2f %6.2f/%.2f\n', names{i}, size(Kn,1), ...
        s.multi(1), s.multi(2), s.multi(3), s.edge, s.coloc, s.colocNull, s.sp25, s.sp25null);
end
r = arrayfun(@(c) corr(sets.orig(c,:)', sets.pois(c,:)', 'rows','complete'), (1:size(sets.orig,1))');
rs = r(sig);
out('\nCorrelation original vs Poisson+VxS kernel: median %.2f (all), %.2f (significant cells)\n', ...
    median(r,'omitnan'), median(rs,'omitnan'));
out('Coverage (main-peak position, 20 cm bins, fraction of cells):\n');
edges = 0:20:200;
for i = 1:numel(names)
    out('   %-8s %s\n', names{i}, mat2str(histcounts(x(S.(names{i}).mainPeak), edges, 'Normalization','probability'), 2));
end
fclose(fid);

% ---- figure
fig = figure('Position',[40 40 1600 800],'Color','w');
tl = tiledlayout(fig, 2, 4, 'TileSpacing','compact');
cols = lines(numel(names));
nexttile; hold on;
for i = 1:numel(names)
    histogram(x(S.(names{i}).mainPeak), 0:10:200, 'Normalization','probability', 'DisplayStyle','stairs', ...
        'EdgeColor', cols(i,:), 'LineWidth', 1.3);
end
xlabel('main peak position (cm)'); ylabel('fraction of cells'); box off; title('Coverage');
legend(names, 'Box','off');
nexttile; hold on;
M = cell2mat(cellfun(@(n) S.(n).multi(:)', names, 'UniformOutput', false));
bar(M'); set(gca,'XTick',1:3,'XTickLabel',{'prom .3','prom .5','prom .7'}); ylabel('fraction multi-peaked'); box off;
title('Multi-peaked kernels'); legend(names,'Box','off','Location','best');
nexttile; hold on;
for i = 1:numel(names)
    plot(i, S.(names{i}).coloc, 'o', 'Color', cols(i,:), 'MarkerFaceColor', cols(i,:));
    plot(i, S.(names{i}).colocNull, 'x', 'Color', cols(i,:), 'MarkerSize', 9, 'LineWidth', 1.5);
end
set(gca,'XTick',1:numel(names),'XTickLabel',names); xlim([.5 numel(names)+.5]); ylim([0 1]);
ylabel('kernel peaks within 4 cm of a VS peak'); box off; title('Co-location (o obs, x shifted null)');
nexttile; hold on;
for i = 1:numel(names)
    histogram(S.(names{i}).dists, 0:2:150, 'Normalization','pdf', 'DisplayStyle','stairs', ...
        'EdgeColor', cols(i,:), 'LineWidth', 1.2);
end
xline(25:25:150, 'k:'); xlabel('distance between peaks (cm)'); ylabel('density'); box off; title('Peak spacing');
nexttile; histogram(r, -1:0.05:1, 'FaceColor', [.6 .6 .6]); hold on;
histogram(rs, -1:0.05:1, 'FaceColor', [.8 .2 .2]); box off;
xlabel('corr(original, Poisson+VxS kernel)'); ylabel('cells'); legend({'all','significant'},'Box','off');
title('Kernel shape agreement');
% example significant cells
[~, ord] = sort(R.LLHrelI_pois .* sig, 'descend');
for k = 1:3
    nexttile; c = ord(k);
    % additive kernels scaled to their max |value| (original and refit units differ)
    nrm = @(v) v ./ max(abs(v), [], 'omitnan');
    yyaxis left; hold on;
    h1 = plot(x, nrm(sets.orig(c,:)), '-', 'Color', [0 0 0]);
    h2 = plot(x, nrm(K.gauss_VSP(c,:)), '-', 'Color', [0 .45 .9]);
    ylabel('additive kernel (scaled to max)'); set(gca, 'YColor', [0 0 0]);
    yyaxis right; h3 = plot(x, K.pois_VSIP(c,:), '-', 'Color', [.85 .1 .1]);
    ylabel('log gain (Poisson+VxS)'); set(gca, 'YColor', [.85 .1 .1]);
    xlabel('position (cm)'); box off;
    title(sprintf('cell %d: LLHrel orig %.2f, pois+VxS %.2f', R.CellIdx(c), R.LLHrel_table(c), R.LLHrelI_pois(c)));
    if k == 1, legend([h1 h2 h3], {'original','Gaussian refit','Poisson+VxS'},'Box','off','Location','best'); end
end
title(tl, 'Do the position-kernel properties survive a multiplicative model?');
exportgraphics(fig, fullfile(C.figRoot, 'kernel_properties.png'), 'Resolution', 150);
end

% ------------------------------------------------------------------------
function s = kernelStats(Kn, V, x)
proms = [0.3 0.5 0.7];
n = size(Kn,1); nP = zeros(n, numel(proms)); mainPeak = zeros(n,1);
pk = cell(n,1); dists = []; colocHit = 0; colocN = 0; nullHit = 0; nullN = 0;
for c = 1:n
    k = Kn(c,:); if all(isnan(k)), mainPeak(c) = NaN; continue; end
    k = fillmissing(k, 'linear', 'EndValues', 'nearest');
    rng_ = max(k) - min(k); if rng_ == 0, mainPeak(c) = NaN; continue; end
    [~, mainPeak(c)] = max(k);
    for j = 1:numel(proms)
        nP(c,j) = numel(findpeaks([min(k) k min(k)], 'MinPeakProminence', proms(j)*rng_));
    end
    [~, loc] = findpeaks([min(k) k min(k)], 'MinPeakProminence', 0.5*rng_); loc = loc - 1;
    pk{c} = loc;
    if numel(loc) > 1, d = abs(loc(:) - loc(:)'); dists = [dists; 2*d(triu(true(numel(loc)),1))]; end %#ok<AGROW>
    v = V(c,:); v = fillmissing(v, 'linear', 'EndValues', 'nearest');
    [~, vl] = findpeaks([min(v) v min(v)], 'MinPeakProminence', 0.2*(max(v)-min(v))); vl = vl - 1;
    if isempty(vl), continue; end
    colocHit = colocHit + sum(min(abs(loc(:) - vl(:)'), [], 2) <= 2); colocN = colocN + numel(loc);
    for r = 1:20      % null: kernel peaks shifted by a random offset (wrapping)
        ls = mod(loc - 1 + randi(99), 100) + 1;
        nullHit = nullHit + sum(min(abs(ls(:) - vl(:)'), [], 2) <= 2); nullN = nullN + numel(ls);
    end
end
v = ~isnan(mainPeak);
s.mainPeak = mainPeak(v);
s.multi = mean(nP(v,:) > 1, 1);
s.edge = mean(mainPeak(v) <= 1);                    % main peak in the first 2 cm bin
s.coloc = colocHit / max(colocN,1); s.colocNull = nullHit / max(nullN,1);
s.dists = dists;
near25 = @(d) mean(abs(d - 25*round(d/25)) <= 3 & d >= 20);
s.sp25 = near25(dists);
s.sp25null = near25(2*abs(randi(100, 1e5, 1) - randi(100, 1e5, 1)));   % uniform peak positions
end
