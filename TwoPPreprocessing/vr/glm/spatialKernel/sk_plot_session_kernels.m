function sk_plot_session_kernels(tag)
% Example kernels for one session (tag, default 'M25126_20260313') from the refit2
% run with the likelihood-ratio test: summary panels plus the top EN cells that
% pass both tests and that pass only the dispersion-corrected LR test. The
% session must be in both refit and refit2.
if nargin < 1, tag = 'M25126_20260313'; end
C = sk_config(); od = C.outRoot;
N = readtable(fullfile(od, 'refit2', ['refit2_' tag '.csv']));
L = load(fullfile(od, 'refit2', ['refit2_' tag '.mat'])); Ex = L.Ex;
O = readtable(fullfile(od, 'refit', ['refit_' tag '.csv']));
OK = load(fullfile(od, 'refit_kernels.mat'));
T = readtable(C.tablePath);
[~, io] = ismember([N.sessionid N.cellNum], [O.sessionid O.cellNum], 'rows'); O = O(io,:);
[~, it] = ismember([N.sessionid N.cellNum], [T.sessionid T.cellNum], 'rows'); cidx = T.CellIdx(it);
[~, ik] = ismember(cidx, OK.CellIdx);
x = Ex.binCentres;
sets = {'old pois VSIP', OK.pois_VSIP(ik,:); 'new ridge VSIP', Ex.kernels.ridge_VSIP; 'new EN VSIP', Ex.kernels.en_VSIP};
sT = N.p_en_posI < 0.05 & N.LLHrelI_en > 1; sQ = N.pLRq_en_posI < 0.05;
grp = 1 + sQ + sT;            % 1 neither, 2 LR/phi only, 3 t & LR/phi (sT without sQ is 1 cell -> grp 2)
fprintf('EN: t&LR/phi %d, LR/phi only %d, t only %d, neither %d\n', nnz(sT&sQ), nnz(sQ&~sT), nnz(sT&~sQ), nnz(~sT&~sQ));

fig = figure('Position',[40 40 1500 1030],'Color','w'); tl = tiledlayout(3,4,'TileSpacing','compact');
cols = [.6 .6 .6; .95 .6 .1; .85 .1 .1];
nexttile; hold on;
for g = 1:3, s = (sQ + sT == g-1) | (g==3 & sT & sQ); end
cat3 = zeros(height(N),1); cat3(sQ & ~sT) = 1; cat3(sT & ~sQ) = 1; cat3(sT & sQ) = 2;
for g = 0:2, s = cat3 == g; scatter(O.LLHrelI_pois(s), N.LLHrelI_en(s), 14, cols(g+1,:), 'filled'); end
refline(1,0); xlabel('old LLHrelI (pois ridge)'); ylabel('new LLHrelI (EN)'); box off; title('Position gain, +VxS');
legend({'neither','one test','t-test & LR/\phi'},'Box','off','Location','northwest');
nexttile; hold on;
for g = 0:2, s = cat3 == g; scatter(N.LLHrelI_ridge(s), N.LLHrelI_en(s), 14, cols(g+1,:), 'filled'); end
refline(1,0); xlabel('ridge'); ylabel('EN'); box off; title('new: ridge vs EN');
nexttile; hold on;
for s = 1:3, Kk = sets{s,2}; [~, mp] = max(Kk,[],2); histogram(x(mp), 0:10:200, 'Normalization','probability','DisplayStyle','stairs','LineWidth',1.3); end
legend(sets(:,1),'Box','off'); xlabel('main peak (cm)'); box off; title('Coverage');
nexttile; hold on;
B = Ex.onOffBasis; nb = size(B,2); lag = Ex.onOffLagS;
W = Ex.onoff.en_VSIP; on = W(:,1:nb)*B'; off = W(:,nb+1:end)*B';
plot(lag, on', 'Color', [0 .45 .9 .15]); plot(lag, median(on), 'b', 'LineWidth', 2);
plot(lag, off', 'Color', [.85 .1 .1 .15]); plot(lag, median(off), 'r', 'LineWidth', 2);
xlabel('time since onset (blue) / offset (red) (s)'); ylabel('log gain'); box off; title('On/off kernels (EN, VSIP)');
rows = {sT & sQ, 't-test & LR/\phi'; sQ & ~sT, 'LR/\phi only'};
for r = 1:2
    sc = N.LLHrelI_en; sc(~rows{r,1}) = -Inf; [~, ord] = sort(sc, 'descend');
    for k = 1:4
        nexttile; c = ord(k); hold on;
        h = [plot(x, sets{1,2}(c,:), 'k'), plot(x, sets{2,2}(c,:), 'Color', [0 .45 .9]), plot(x, sets{3,2}(c,:), 'Color', [.85 .1 .1])];
        xlabel('position (cm)'); ylabel('log gain'); box off;
        fu = 10*mean(Ex.foldLL(c,:,4,1) > Ex.foldLL(c,:,3,1));
        title({sprintf('%s: cell %d', rows{r,2}, cidx(c)), [sprintf('LLHrelI %.2f, ', N.LLHrelI_en(c)) '\phi ' sprintf('%.1f, %d/10 folds better', N.phi_en_posI(c), round(fu))]}, 'FontWeight','normal');
        if r == 1 && k == 1, legend(h, {'old ridge','new ridge','new EN'},'Box','off'); end
    end
end
title(tl, ['Refit v2 with likelihood-ratio test (measured \phi), ' strrep(tag,'_',' ') ': top EN cells per test']);
exportgraphics(fig, fullfile(C.figRoot, ['refit2_LR_kernels_' tag '.png']), 'Resolution', 150);
end
