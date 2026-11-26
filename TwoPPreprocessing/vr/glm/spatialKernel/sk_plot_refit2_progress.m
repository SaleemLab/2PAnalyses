function sk_plot_refit2_progress()
% Results so far from the full refit2 run: counts per session and sorted kernel images
C = sk_config(); od = fullfile(C.outRoot, 'refit2');
T = readtable(C.tablePath);
Korig = h5read('D:\DATA\Sonali\GLM\ForAS\V1_GLMData.mat', '/Kernels/PosKernel');
d = dir(fullfile(od, 'refit2_M*.csv')); ns = numel(d);
nrm = @(K) (K - min(K,[],2)) ./ max(range(K,2), eps);
cm = interp1([0 .5 1], [1 1 1; .55 .75 .9; .03 .19 .42], linspace(0,1,256));
fig = figure('Position',[40 40 1700 430*ns],'Color','w'); tl = tiledlayout(ns,4,'TileSpacing','compact');
fprintf('%-18s %5s | %-28s | %-28s\n', 'session', 'cells', 'EN: t / LR/phi / both', 'ridge: t / LR/phi / both');
tot = zeros(1,7);
for s = 1:ns
    N = readtable(fullfile(od, d(s).name)); L = load(fullfile(od, strrep(d(s).name,'.csv','.mat'))); Ex = L.Ex;
    tag = erase(erase(d(s).name, 'refit2_'), '.csv'); n = height(N);
    [~, it] = ismember([N.sessionid N.cellNum], [T.sessionid T.cellNum], 'rows'); cidx = T.CellIdx(it);
    c = zeros(1,6); k = 0;
    for p = {'en','ridge'}
        sT = N.(['p_' p{1} '_posI']) < 0.05 & N.(['LLHrelI_' p{1}]) > 1; sQ = N.(['pLRq_' p{1} '_posI']) < 0.05;
        c(k+(1:3)) = [nnz(sT) nnz(sQ) nnz(sT&sQ)]; k = k + 3;
    end
    tot = tot + [n c];
    fprintf('%-18s %5d | %4d (%2.0f%%) %4d (%2.0f%%) %4d | %4d (%2.0f%%) %4d (%2.0f%%) %4d   phi median %.2f\n', tag, n, ...
        c(1), 100*c(1)/n, c(2), 100*c(2)/n, c(3), c(4), 100*c(4)/n, c(5), 100*c(5)/n, c(6), median(N.phi_en_posI));
    x = Ex.binCentres; Ko = Korig(cidx,:); Kr = Ex.kernels.ridge_VSIP; Ke = Ex.kernels.en_VSIP;
    [~, po] = max(Ko,[],2); [~, oA] = sort(po);
    sQ = N.pLRq_en_posI < 0.05; [~, pe] = max(Ke,[],2); iS = find(sQ); [~, o] = sort(pe(iS)); oS = iS(o);
    P = {Ko(oA,:), 'Original GLM1'; Kr(oA,:), 'Ridge, +VxS'; Ke(oA,:), 'Elastic net, +VxS'; ...
         Ke(oS,:), sprintf('EN, LR/phi significant (n = %d), own sort', numel(oS))};
    for j = 1:4
        nexttile; K = nrm(P{j,1}); K(isnan(K)) = 0;
        imagesc(x, 1:size(K,1), K); colormap(cm); clim([0 1]); box off; set(gca,'TickDir','out');
        xlabel('Position (cm)'); title(P{j,2}, 'FontWeight','normal');
        if j == 1, ylabel({strrep(tag,'_',' '), sprintf('%d cells, sorted by original peak', n)}); end
    end
end
n = tot(1);
fprintf('%-18s %5d | %4d (%2.0f%%) %4d (%2.0f%%) %4d | %4d (%2.0f%%) %4d (%2.0f%%) %4d\n', 'total', n, ...
    tot(2), 100*tot(2)/n, tot(3), 100*tot(3)/n, tot(4), tot(5), 100*tot(5)/n, tot(6), 100*tot(6)/n, tot(7));
title(tl, sprintf('Refit v2, sessions finished so far (%d cells); each row scaled to its own min-max', n));
exportgraphics(fig, fullfile(C.figRoot, 'refit2_sofar_kernel_images.png'), 'Resolution', 130);
end
