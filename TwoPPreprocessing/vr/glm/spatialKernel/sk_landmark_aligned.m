function A = sk_landmark_aligned(sets, nShuf)
% Are position kernels systematically higher just before or after the textures, or at a
% fixed phase of the 24 cm background period? For each kernel set (significant cells):
%   - each kernel (G) is z-scored and cut in a +-20 cm window around each texture
%     (40, 80, 120, 160 cm); the mean over cells is compared with a null in which each
%     kernel is shifted circularly by a random amount over 10-200 cm (onset bins kept out)
%   - the same windows of GLM1's vision + speed prediction, measured map and mean speed
%     show where the visual drive and behaviour change around each texture
%   - the kernels folded modulo 24 cm (phase of the background period), away from the
%     textures, against the same null
%   - signed distances of all kernel peaks to the nearest texture against the null
% Writes figures/landmark_aligned.png and landmark_aligned_summary.txt.
if nargin < 1, sets = {'glm1', 'refit1_vsp', 'refit1_vsip', 'refit2_en_vsip'}; end
if nargin < 2, nShuf = 500; end
C = sk_config(); rng(1); tex = [40 80 120 160]; w = 10;           % window +-10 bins = +-20 cm
lagX = 2*(-w:w); ns = numel(sets); A = struct();
fig = figure('Position', [20 20 1900 330*ns], 'Color', 'w'); tl = tiledlayout(ns, 5, 'TileSpacing', 'compact');
fid = fopen(fullfile(C.outRoot, 'landmark_aligned_summary.txt'), 'w');
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('Texture-aligned kernels (z-scored G), significant cells; null: circular shifts over 10-200 cm, %d shuffles\n', nShuf);
for s = 1:ns
    ks = sk_kernelset(sets{s}); x = ks.x; m = ks.sig;
    zs = @(M) (M - mean(M, 2, 'omitnan')) ./ std(M, 0, 2, 'omitnan');
    Z = zs(ks.G(m,:)); V = zs(ks.comp.mapVS(m,:)); O = zs(ks.comp.mapObs(m,:)); S = zs(ks.comp.speedMap(m,:));
    ci = arrayfun(@(L) find(x == L - 1) + (-w:w), tex, 'UniformOutput', false);   % bin with centre L-1 (39, 79, ...)
    win = @(M) cell2mat(cellfun(@(c) mean(M(:, c), 1, 'omitnan'), ci', 'UniformOutput', false));   % 4 x 21
    real = win(Z); rv = win(V); ro = win(O); rs = win(S);
    % null: shift each cell's kernel over bins with x > 10
    sel = find(x > 10); nl = zeros(nShuf, 4, numel(lagX)); fold = @(M) foldp(M, x); nf = zeros(nShuf, 12);
    pkReal = peakDist(Z, x, tex); pkNull = [];
    for k = 1:nShuf
        Zs = Z; sh = randi([5 numel(sel) - 5], size(Z, 1), 1);
        for c = 1:size(Z, 1), Zs(c, sel) = circshift(Z(c, sel), sh(c)); end
        nl(k,:,:) = win(Zs); nf(k,:) = fold(Zs);
        if k <= 50, pkNull = [pkNull; peakDist(Zs, x, tex)]; end %#ok<AGROW>
    end
    lo = squeeze(prctile(nl, 2.5, 1)); hi = squeeze(prctile(nl, 97.5, 1));
    % mean over textures 2-4 (the first texture comes right after the onset region)
    r234 = mean(real(2:4,:), 1); n234 = squeeze(mean(nl(:, 2:4, :), 2));
    lo234 = prctile(n234, 2.5); hi234 = prctile(n234, 97.5);
    pre = lagX >= -8 & lagX < 0; post = lagX > 0 & lagX <= 8;
    dReal = mean(r234(pre)) - mean(r234(post)); dNull = mean(n234(:, pre), 2) - mean(n234(:, post), 2);
    pDiff = mean(abs(dNull) >= abs(dReal));
    fr = fold(Z); pF = mean(range(nf, 2) >= range(fr));
    out('\n%s (%d cells)\n', sets{s}, nnz(m));
    out('  textures 2-4: mean z 0-8 cm before %.3f, 0-8 cm after %.3f; before - after %.3f (null 95%%: %.3f to %.3f, p = %.3f)\n', ...
        mean(r234(pre)), mean(r234(post)), dReal, prctile(dNull, 2.5), prctile(dNull, 97.5), pDiff);
    for t = 1:4
        [mx, j] = max(real(t,:)); out('  texture at %3d cm: max %.3f at %+d cm (null upper 95%% %.3f); bins above null: %s\n', tex(t), mx, lagX(j), hi(t, j), ...
            mat2str(lagX(real(t,:) > hi(t,:))));
    end
    out('  24 cm phase fold: range %.3f (null p = %.3f), peak at phase %d cm\n', range(fr), pF, 2*(find(fr == max(fr), 1) - 1));
    out('  peaks within 8 cm before / after a texture: %.1f%% / %.1f%% of peaks (null %.1f%% / %.1f%%)\n', ...
        100*mean(pkReal >= -8 & pkReal < 0), 100*mean(pkReal > 0 & pkReal <= 8), 100*mean(pkNull >= -8 & pkNull < 0), 100*mean(pkNull > 0 & pkNull <= 8));
    A.(sets{s}) = struct('real', real, 'lo', lo, 'hi', hi, 'vis', rv, 'obs', ro, 'speed', rs, 'fold', fr, 'pDiff', pDiff, 'pFold', pF);

    % figure row
    cols = lines(4);
    nexttile; hold on;
    for t = 1:4, plot(lagX, real(t,:), 'Color', cols(t,:), 'LineWidth', 1.3); end
    xline(0, 'k:'); xlabel('cm from texture'); ylabel(sprintf('%s\nmean kernel (z)', strrep(sets{s}, '_', ' ')));
    if s == 1, legend(compose('%d cm', tex), 'Box', 'off', 'Location', 'best'); end; box off; title('Kernel around each texture');
    nexttile; hold on; fill([lagX fliplr(lagX)], [lo234 fliplr(hi234)], [.85 .85 .85], 'EdgeColor', 'none');
    plot(lagX, r234, 'k', 'LineWidth', 1.6); xline(0, 'k:'); xlabel('cm from texture'); box off;
    title(sprintf('Textures 2-4 vs shuffle (before-after p = %.2f)', pDiff));
    nexttile; hold on; plot(lagX, mean(rv(2:4,:)), 'Color', [.85 .25 .15], 'LineWidth', 1.4);
    plot(lagX, mean(ro(2:4,:)), 'Color', [.2 .45 .8], 'LineWidth', 1.4); plot(lagX, mean(rs(2:4,:)), 'Color', [.3 .6 .3], 'LineWidth', 1.4);
    plot(lagX, r234, 'k', 'LineWidth', 1.4); xline(0, 'k:'); xlabel('cm from texture'); box off; title('Textures 2-4: other signals (z)');
    if s == 1, legend({'GLM1 vision+speed', 'measured map', 'running speed', 'kernel'}, 'Box', 'off', 'Location', 'best'); end
    nexttile; hold on; ph = 0:2:22; fill([ph fliplr(ph)], [prctile(nf, 2.5) fliplr(prctile(nf, 97.5))], [.85 .85 .85], 'EdgeColor', 'none');
    plot(ph, fr, 'k-o', 'LineWidth', 1.3, 'MarkerSize', 3); xlabel('position mod 24 (cm)'); box off; title(sprintf('24 cm phase fold (p = %.2f)', pF));
    nexttile; hold on; ed = -20:2:20;
    histogram(pkNull, ed, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', [.6 .6 .6], 'LineWidth', 1.3);
    histogram(pkReal, ed, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1.5);
    xline(0, 'k:'); xlabel('peak - nearest texture (cm)'); ylabel('fraction of peaks'); box off; title('Kernel peaks near textures');
    if s == 1, legend({'shuffled', 'real'}, 'Box', 'off'); end
end
fclose(fid);
title(tl, 'Position kernels aligned to textures (significant cells; textures pass the animal at 37-43, 77-83, 117-123, 157-163 cm)');
exportgraphics(fig, fullfile(C.figRoot, 'landmark_aligned.png'), 'Resolution', 120);
end

function f = foldp(M, x)
% mean over cells and over bins with x > 10 at each 2 cm phase of the 24 cm period,
% leaving out bins within 10 cm of a texture (textures sit at fixed phases of 24 cm)
ph = mod(x - 1, 24); sel = x > 10 & min(abs(x(:) - [39 79 119 159]), [], 2)' > 10; f = zeros(1, 12);
for p = 0:11, f(p+1) = mean(M(:, sel & ph == 2*p), 'all', 'omitnan'); end
end

function d = peakDist(Z, x, tex)
% signed distance (cm) of every peak (prominence > 50% of range, x > 10) to the nearest texture
d = [];
for c = 1:size(Z, 1)
    g = Z(c,:); rg = range(g); if ~(rg > 0), continue; end
    [~, lp] = findpeaks([min(g) g min(g)], 'MinPeakProminence', 0.5*rg); p = x(lp - 1); p = p(p > 10);
    dd = p(:) - tex; [~, j] = min(abs(dd), [], 2); d = [d; dd(sub2ind(size(dd), (1:numel(p))', j))]; %#ok<AGROW>
end
end
