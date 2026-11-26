function sk_kernel_shapes()
% Position-kernel shapes: additive (Poisson VSP) vs multiplicative (Poisson VSIP), same 892 cells, earlier refit
C = sk_config(); od = C.outRoot; rng(1);
R = readtable(fullfile(od, 'refit_spatialcells.csv'));
K = load(fullfile(od, 'refit_kernels.mat')); assert(isequal(K.CellIdx, R.CellIdx));
Korig = h5read('D:\DATA\Sonali\GLM\ForAS\V1_GLMData.mat', '/Kernels/PosKernel'); Korig = Korig(K.CellIdx,:);
x = K.binCentres;                                   % 1:2:199 cm
texLM = [40 80 120 160];                            % wall textures pass the animal (37-43, 77-83, 117-123, 157-163 cm)
sigM = R.p_pois_posI < 0.05 & R.LLHrelI_pois > 1;
sets = {'GLM1 (additive, Gaussian)', Korig; 'Poisson additive (VSP)', K.pois_VSP; 'Poisson multiplicative (VSIP)', K.pois_VSIP};
ns = size(sets,1); n = height(R);

% per-cell metrics
F = struct();
for s = 1:ns
    Kk = fillmissing(sets{s,2}, 'linear', 2, 'EndValues', 'nearest');
    f.depth = max(Kk,[],2) - min(Kk,[],2);
    med = median(Kk, 2); sd = std(Kk, 0, 2);
    [mx, mp] = max(Kk, [], 2); f.mainPeak = x(mp)';
    f.peakZ = (mx - med) ./ sd;                                       % main peak height in SDs
    f.fwhm = nan(n,1); f.nPk = nan(n,1); pk = cell(n,1);
    for c = 1:n
        k = Kk(c,:); h = med(c) + (mx(c) - med(c))/2;
        l = mp(c); while l > 1 && k(l-1) > h, l = l - 1; end
        r = mp(c); while r < numel(k) && k(r+1) > h, r = r + 1; end
        f.fwhm(c) = 2*(r - l + 1);
        rg = range(k); if rg == 0, pk{c} = []; f.nPk(c) = 0; continue; end
        [~, loc] = findpeaks([min(k) k min(k)], 'MinPeakProminence', 0.5*rg); pk{c} = x(loc - 1);
        f.nPk(c) = numel(loc);
    end
    g = exp(Kk - max(Kk,[],2)); f.sparse = 1 - (mean(g,2).^2 ./ mean(g.^2,2));    % 0 flat, ->1 one sharp peak
    f.pk = pk; F(s).f = f; %#ok<SAGROW>
end
rAM = arrayfun(@(c) corr(K.pois_VSP(c,:)', K.pois_VSIP(c,:)', 'rows','complete'), (1:n)');
shift = abs(F(3).f.mainPeak - F(2).f.mainPeak);

% landmark and periodicity measures on pooled peaks
near = @(p, L, w) mean(min(abs(p(:) - L(:)'), [], 2) <= w);
uni = near(x, texLM, 4);                                            % expected fraction if peaks were uniform
vs = @(p, P) abs(mean(exp(2i*pi*p(:)/P)));                          % vector strength at period P
fid = fopen(fullfile(od, 'kernel_shapes_summary.txt'), 'w');
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('Position-kernel shape, %d GLM1 spatial cells (sig = multiplicative position significant, n = %d)\n', n, nnz(sigM));
out('Texture landmarks at %s cm; background elements repeat every 24 cm\n\n', mat2str(texLM));
hdr = '%-31s %7s %7s %7s %7s %7s %7s %8s %8s %8s %8s\n';
out(hdr, 'set', 'depth', 'peakZ', 'FWHMcm', 'nPeaks', 'multi%', 'sparse', 'pk<10cm', 'nearTex', 'VS24cm', 'VS40cm');
for sub = 1:2
    if sub == 1, m = true(n,1); out('All cells\n'); else, m = sigM; out('Multiplicative-significant cells\n'); end
    for s = 1:ns
        f = F(s).f; p = [f.pk{m}]; p = p(p > 10);                    % ignore onset-edge peaks for landmark tests
        out('%-31s %7.2f %7.2f %7.1f %7.2f %7.0f %7.2f %8.2f %8.2f %8.3f %8.3f\n', sets{s,1}, median(f.depth(m)), median(f.peakZ(m)), ...
            median(f.fwhm(m)), mean(f.nPk(m)), 100*mean(f.nPk(m) > 1), median(f.sparse(m)), mean(f.mainPeak(m) < 10), ...
            near(p, texLM, 4), vs(p, 24), vs(p, 40));
    end
end
out('\nUniform expectation: nearTex %.2f; vector strength null (n peaks ~%d) ~%.3f\n', uni, numel([F(3).f.pk{:}]), sqrt(pi/(4*numel([F(3).f.pk{:}]))));
out('Additive vs multiplicative kernel, same cell: median corr %.2f (sig cells %.2f); main peak moves >10 cm in %.0f%% (sig %.0f%%)\n', ...
    median(rAM,'omitnan'), median(rAM(sigM),'omitnan'), 100*mean(shift > 10), 100*mean(shift(sigM) > 10));
pw = @(a,b) signrank(a, b);
out('Paired additive vs multiplicative (Poisson, all cells), signrank p: depth %.1g, peakZ %.1g, FWHM %.1g, nPeaks %.1g, sparse %.1g\n', ...
    pw(F(2).f.depth, F(3).f.depth), pw(F(2).f.peakZ, F(3).f.peakZ), pw(F(2).f.fwhm, F(3).f.fwhm), pw(F(2).f.nPk, F(3).f.nPk), pw(F(2).f.sparse, F(3).f.sparse));
% where did the gain go: mean kernel difference along the track (additive - multiplicative), mean-centred per cell
dK = (K.pois_VSP - mean(K.pois_VSP,2,'omitnan')) - (K.pois_VSIP - mean(K.pois_VSIP,2,'omitnan'));
fclose(fid);

% figure
cs = [.4 .4 .4; .2 .45 .8; .85 .25 .15];
fig = figure('Position',[30 30 1750 900],'Color','w'); tl = tiledlayout(2,4,'TileSpacing','compact');
mets = {'depth','modulation depth (log gain)'; 'peakZ','main peak height (SD)'; 'fwhm','main peak FWHM (cm)'; 'mainPeak','main peak position (cm)'};
for i = [2 3 4]
    nexttile; hold on; v = cellfun(@(f) f.(mets{i,1}), {F.f}, 'UniformOutput', false);
    for s = 1:ns
        q = v{s}; e = prctile(q, [1 99]); ed = linspace(min(cellfun(@(z) prctile(z,1), v)), max(cellfun(@(z) prctile(z,99), v)), 40);
        histogram(q, ed, 'DisplayStyle','stairs','EdgeColor',cs(s,:),'LineWidth',1.4,'Normalization','probability');
    end
    xlabel(mets{i,2}); ylabel('fraction of cells'); box off; title(mets{i,2});
    if i == 2, legend(sets(:,1),'Box','off'); end
end
nexttile; hold on; ct = zeros(ns, 5);
for s = 1:ns, ct(s,:) = histcounts(min(F(s).f.nPk, 5), 0.5:1:5.5, 'Normalization','probability'); end
b = bar(1:5, ct', 'grouped'); for s = 1:ns, b(s).FaceColor = cs(s,:); end
set(gca,'XTickLabel',{'1','2','3','4','5+'}); xlabel('peaks per kernel (prominence > 50% of range)'); ylabel('fraction of cells'); box off; title('Multi-peakedness');
% pooled peak density with landmarks
nexttile([1 2]); hold on; ed = 0:4:200;
for L = texLM, patch([L-3 L+3 L+3 L-3], [0 0 .06 .06], [.95 .85 .3], 'EdgeColor','none', 'FaceAlpha', .5); end
for b24 = 0:24:192, xline(b24 + 1, ':', 'Color', [.6 .6 .6]); end
for s = 1:ns
    p = [F(s).f.pk{:}]; histogram(p, ed, 'DisplayStyle','stairs','EdgeColor',cs(s,:),'LineWidth',1.5,'Normalization','probability');
end
xlabel('position (cm)'); ylabel('fraction of all kernel peaks'); box off; xlim([0 200]); ylim([0 .06]);
title('Where kernel peaks fall (yellow: texture passes the animal; dotted: 24 cm background period)');
nexttile; hold on;
mA = median(K.pois_VSP - median(K.pois_VSP,2), 1, 'omitnan'); mM = median(K.pois_VSIP - median(K.pois_VSIP,2), 1, 'omitnan');
for L = texLM, patch([L-3 L+3 L+3 L-3], [-.05 -.05 .15 .15], [.95 .85 .3], 'EdgeColor','none', 'FaceAlpha', .5); end
plot(x, mA, 'Color', cs(2,:), 'LineWidth', 1.5); plot(x, mM, 'Color', cs(3,:), 'LineWidth', 1.5);
plot(x, median(dK,1,'omitnan'), 'k', 'LineWidth', 1);
xlabel('position (cm)'); ylabel('median log gain re cell median'); box off; ylim([-.05 .15]);
legend({'','','','','additive','multiplicative','add - mult'},'Box','off'); title('Population kernel');
title(tl, sprintf('Position-kernel shape: additive vs multiplicative (Poisson ridge refit, %d GLM1 spatial cells)', n));
exportgraphics(fig, fullfile(C.figRoot, 'kernel_shapes_add_vs_mult.png'), 'Resolution', 140);
end
