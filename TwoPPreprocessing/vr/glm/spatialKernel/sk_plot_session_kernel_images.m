function sk_plot_session_kernel_images(tag)
% Sorted kernel images (cells as rows) for one session (tag, default
% 'M25126_20260313'): original GLM1, earlier Poisson ridge refit, and refit2
% ridge and EN (+VxS). The session must be in both refit and refit2.
if nargin < 1, tag = 'M25126_20260313'; end
C = sk_config(); od = C.outRoot;
N = readtable(fullfile(od, 'refit2', ['refit2_' tag '.csv']));
L = load(fullfile(od, 'refit2', ['refit2_' tag '.mat'])); Ex = L.Ex;
OK = load(fullfile(od, 'refit_kernels.mat'));
T = readtable(C.tablePath);
[~, it] = ismember([N.sessionid N.cellNum], [T.sessionid T.cellNum], 'rows'); cidx = T.CellIdx(it);
[~, ik] = ismember(cidx, OK.CellIdx);
Korig = h5read('D:\DATA\Sonali\GLM\ForAS\V1_GLMData.mat', '/Kernels/PosKernel');
x = Ex.binCentres;
sets = {'Original GLM1 (lasso, Gaussian)', Korig(cidx,:); 'Poisson ridge, +VxS (earlier refit)', OK.pois_VSIP(ik,:); ...
        'Ridge, +VxS (this run)', Ex.kernels.ridge_VSIP; 'Elastic net, +VxS (this run)', Ex.kernels.en_VSIP};
sQ = N.pLRq_en_posI < 0.05;
nrm = @(K) (K - min(K,[],2)) ./ max(range(K,2), eps);
[~, pOrig] = max(sets{1,2}, [], 2); [~, oAll] = sort(pOrig);
[~, pEN] = max(sets{4,2}, [], 2); iS = find(sQ); [~, o] = sort(pEN(iS)); oSig = iS(o);
cm = interp1([0 .5 1], [1 1 1; .55 .75 .9; .03 .19 .42], linspace(0,1,256));
fig = figure('Position',[40 40 1700 900],'Color','w'); tl = tiledlayout(2,4,'TileSpacing','compact');
rows = {oAll, sprintf('All spatial cells (n = %d), sorted by original GLM1 peak', numel(oAll)); ...
        oSig, ['EN position significant, LR/phi ' sprintf('(n = %d), sorted by EN peak', numel(oSig))]};
for r = 1:2
    for s = 1:4
        nexttile; K = nrm(sets{s,2}(rows{r,1},:)); K(isnan(K)) = 0;
        imagesc(x, 1:size(K,1), K); colormap(cm); clim([0 1]); box off; set(gca,'TickDir','out');
        xlabel('Position (cm)'); title(sets{s,1}, 'FontWeight','normal');
        if s == 1, ylabel({rows{r,2}, 'cells'}); end
    end
end
title(tl, ['Position kernels, ' strrep(tag,'_',' ') ' (each row scaled to its own min-max)']);
exportgraphics(fig, fullfile(C.figRoot, ['refit2_kernel_images_' tag '.png']), 'Resolution', 150);
end
