function sk_add_vs_mult()
% Additive (VS -> VSP) vs multiplicative (VSI -> VSIP) position models: recovered cells and gain from position
C = sk_config(); od = C.outRoot;
R = readtable(fullfile(od, 'refit_spatialcells.csv'));          % Poisson ridge, all 4 models, 892 cells
K = load(fullfile(od, 'refit_kernels.mat')); assert(isequal(K.CellIdx, R.CellIdx));
E = readtable(fullfile(od, 'refit2_full', 'refit2_M25126_20260313.csv'));   % EN, all models, session 33
sid = unique(R.sessionid);
sname = cell(numel(sid),1);
for i = 1:numel(sid)
    j = find([C.sessions{:,3}] == sid(i)); sname{i} = sprintf('%s %d', C.sessions{j,1}, C.sessions{j,2});
end
% phi from the current run where available (VSIP, ridge); the LR test uses 100 dof
phi = nan(height(R),1);
for f = dir(fullfile(od, 'refit2', 'refit2_M*.csv'))'
    N = readtable(fullfile(f.folder, f.name)); [tf, ia] = ismember([R.sessionid R.cellNum], [N.sessionid N.cellNum], 'rows');
    phi(tf) = N.phi_ridge_posI(ia(tf));
end
dF = @(m1,m2) 2*sum(K.foldLL(:,:,m2,2) - K.foldLL(:,:,m1,2), 2);      % dist 2 = Poisson
lrA = chi2cdf(max(dF(1,2),0)./phi, 100, 'upper') < 0.05;
lrM = chi2cdf(max(dF(3,4),0)./phi, 100, 'upper') < 0.05;
tA = R.p_pois_pos < 0.05 & R.LLHrel_pois > 1; tM = R.p_pois_posI < 0.05 & R.LLHrelI_pois > 1;
hasPhi = ~isnan(phi);
fid = fopen(fullfile(od, 'add_vs_mult_summary.txt'), 'w');
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('Poisson ridge refit, GLM1 spatial cells. Additive: VSP vs VS. Multiplicative: VSIP vs VSI.\n');
out('%-18s %5s | %-22s | %-22s | %s\n', 'session', 'cells', 't-test add / mult', 'LR/phi add / mult', 'median LLHrel add / mult');
for i = 0:numel(sid)
    if i == 0, s = true(height(R),1); nm = 'all'; else, s = R.sessionid == sid(i); nm = sname{i}; end
    sp = s & hasPhi;
    if any(sp), lrs = sprintf('%4d (%2.0f%%) %4d (%2.0f%%)', nnz(lrA&sp), 100*mean(lrA(sp)), nnz(lrM&sp), 100*mean(lrM(sp)));
    else, lrs = sprintf('%22s', 'n/a (no phi yet)'); end
    out('%-18s %5d | %4d (%2.0f%%) %4d (%2.0f%%) | %s | %.3f / %.3f\n', nm, nnz(s), nnz(tA&s), 100*mean(tA(s)), nnz(tM&s), 100*mean(tM(s)), ...
        lrs, median(R.LLHrel_pois(s)), median(R.LLHrelI_pois(s)));
end
out('\nt-test categories (all cells): both %d, additive only %d, multiplicative only %d, neither %d\n', nnz(tA&tM), nnz(tA&~tM), nnz(~tA&tM), nnz(~tA&~tM));
out('LR/phi categories (%d cells with phi): both %d, additive only %d, multiplicative only %d, neither %d\n', nnz(hasPhi), ...
    nnz(lrA&lrM&hasPhi), nnz(lrA&~lrM&hasPhi), nnz(~lrA&lrM&hasPhi), nnz(~lrA&~lrM&hasPhi));
gA = R.LLHi_pois_VSP - R.LLHi_pois_VS; gM = R.LLHi_pois_VSIP - R.LLHi_pois_VSI;
out('Position gain (bits/spike): median additive %.4f, multiplicative %.4f; per-cell mult/add ratio median %.2f\n', median(gA), median(gM), median(gM(gA>0)./gA(gA>0)));
out('Base model: median LLHi VS %.3f, VSI %.3f (+%.0f%%); full: VSP %.3f, VSIP %.3f\n', median(R.LLHi_pois_VS), median(R.LLHi_pois_VSI), ...
    100*(median(R.LLHi_pois_VSI./R.LLHi_pois_VS)-1), median(R.LLHi_pois_VSP), median(R.LLHi_pois_VSIP));
eA = E.p_en_pos < 0.05 & E.LLHrel_en > 1; eM = E.p_en_posI < 0.05 & E.LLHrelI_en > 1;
out('\nElastic net, M25126 20260313 (%d cells): t-test additive %d, multiplicative %d; median LLHrel %.3f / %.3f\n', height(E), nnz(eA), nnz(eM), median(E.LLHrel_en), median(E.LLHrelI_en));
fclose(fid);

% figure
cA = [.2 .45 .8]; cM = [.85 .25 .15]; cc = [.7 .7 .7; cA; cM; .45 .2 .55];
fig = figure('Position',[40 40 1700 820],'Color','w'); tl = tiledlayout(2,4,'TileSpacing','compact');
nexttile; fr = zeros(numel(sid)+1, 2);
for i = 1:numel(sid), s = R.sessionid == sid(i); fr(i,:) = [mean(tA(s)) mean(tM(s))]; end
fr(end,:) = [mean(tA) mean(tM)];
b = bar(fr, 'grouped'); b(1).FaceColor = cA; b(2).FaceColor = cM;
set(gca, 'XTick', 1:numel(sid)+1, 'XTickLabel', [sname; {'all'}]); xtickangle(20); ylim([0 1]);
ylabel('fraction of GLM1 spatial cells'); legend({'additive','multiplicative'},'Box','off'); box off; title('Recovered (fold t-test)');

nexttile; fr2 = nan(numel(sid)+1, 2);
for i = 1:numel(sid), s = R.sessionid == sid(i) & hasPhi; if any(s), fr2(i,:) = [mean(lrA(s)) mean(lrM(s))]; end; end
fr2(end,:) = [mean(lrA(hasPhi)) mean(lrM(hasPhi))];
b = bar(fr2, 'grouped'); b(1).FaceColor = cA; b(2).FaceColor = cM;
set(gca, 'XTick', 1:numel(sid)+1, 'XTickLabel', [sname; {'all with \phi'}]); xtickangle(20); ylim([0 1]); box off;
title(['Recovered (LR/\phi, ' sprintf('%d cells with ', nnz(hasPhi)) '\phi)']);

nexttile; ct = [nnz(tA&tM) nnz(tA&~tM) nnz(~tA&tM) nnz(~tA&~tM)];
b = bar(ct, 'FaceColor', 'flat'); b.CData = [cc(4,:); cA; cM; cc(1,:)];
set(gca, 'XTickLabel', {'both','add only','mult only','neither'}); ylabel('cells'); box off; title('Overlap (t-test, all sessions)');
text(1:4, ct, string(ct), 'HorizontalAlignment','center','VerticalAlignment','bottom');

nexttile; hold on; M = [R.LLHi_pois_VS R.LLHi_pois_VSP R.LLHi_pois_VSI R.LLHi_pois_VSIP];
md = median(M); ci = bootci(1000, @median, M);
bar(1:4, md, 'FaceColor', 'flat', 'CData', [cA*.5+.5; cA; cM*.5+.5; cM]);
errorbar(1:4, md, md-ci(1,:), ci(2,:)-md, 'k.', 'LineStyle','none');
set(gca, 'XTick', 1:4, 'XTickLabel', {'VS','VSP','VSI','VSIP'}); ylabel('median held-out info (bits/spike)'); box off; title('Model performance');

nexttile; hold on; cat = 1 + tA + 2*tM;
for k = 1:4, s = cat == k; scatter(R.LLHrel_pois(s), R.LLHrelI_pois(s), 8, cc(k,:), 'filled', 'MarkerFaceAlpha', .6); end
lim = [0.9 1.6]; plot(lim, lim, 'k:'); xlim(lim); ylim(lim); axis square;
xlabel('LLHrel additive (VSP/VS)'); ylabel('LLHrel multiplicative (VSIP/VSI)'); box off; title('Relative gain from position');
legend({'neither','add only','mult only','both'},'Box','off','Location','northwest');

nexttile; hold on; scatter(gA, gM, 8, 'k', 'filled', 'MarkerFaceAlpha', .3);
lim = [-0.005 0.06]; plot(lim, lim, 'k:'); xlim(lim); ylim(lim); axis square;
xlabel('\DeltaLLH position, additive (bits/spike)'); ylabel('\DeltaLLH position, multiplicative'); box off; title('Absolute gain from position');

nexttile; hold on; ed = linspace(-0.01, 0.06, 50);
histogram(gA, ed, 'DisplayStyle','stairs','EdgeColor',cA,'LineWidth',1.5); histogram(gM, ed, 'DisplayStyle','stairs','EdgeColor',cM,'LineWidth',1.5);
xline(0,'k:'); xlabel('\DeltaLLH position (bits/spike)'); ylabel('cells'); legend({'additive','multiplicative'},'Box','off'); box off; title('Gain from position');

nexttile; hold on; cat = 1 + eA + 2*eM;
for k = 1:4, s = cat == k; scatter(E.LLHrel_en(s), E.LLHrelI_en(s), 14, cc(k,:), 'filled'); end
lim = [0.9 1.5]; plot(lim, lim, 'k:'); xlim(lim); ylim(lim); axis square;
xlabel('LLHrel additive'); ylabel('LLHrel multiplicative'); box off;
title(sprintf('Elastic net, M25126 20260313: %d add / %d mult', nnz(eA), nnz(eM)));
title(tl, sprintf('Additive vs multiplicative position models, Poisson ridge refit (%d GLM1 spatial cells)', height(R)));
exportgraphics(fig, fullfile(C.figRoot, 'additive_vs_multiplicative.png'), 'Resolution', 140);
end
