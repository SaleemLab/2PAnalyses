function sk_summarise_refit(csvPath)
% Summarise the refits (steps 1-2): does position still add cross-validated
% likelihood once the vision+speed model is Poisson (multiplicative) and/or
% has an explicit vision x speed term? Writes refit_summary.txt and a figure.
C = sk_config();
if nargin < 1
    csvPath = fullfile(C.outRoot, 'refit_spatialcells.csv');
    if ~exist(csvPath,'file')   % partial run: combine whatever sessions exist
        d = dir(fullfile(C.outRoot,'refit','refit_M*.csv'));
        R = cellfun(@(f) readtable(fullfile(C.outRoot,'refit',f)), {d.name}, 'UniformOutput', false);
        R = vertcat(R{:});
        T = readtable(C.tablePath);
        [~, ia] = ismember([R.sessionid R.cellNum], [T.sessionid T.cellNum], 'rows');
        R.CellIdx = T.CellIdx(ia); R.LLHrel_table = T.LLHrel(ia);
        R.LLHi_vis_table = T.LLHi_vis(ia); R.LLHi_pos_table = T.LLHi_pos(ia);
    end
end
if ~exist('R','var'), R = readtable(csvPath); end

% position gain in bits/spike, and relative to vision+speed
dG  = R.LLHi_gauss_VSP - R.LLHi_gauss_VS;   dGI = R.LLHi_gauss_VSIP - R.LLHi_gauss_VSI;
dP  = R.LLHi_pois_VSP  - R.LLHi_pois_VS;    dPI = R.LLHi_pois_VSIP  - R.LLHi_pois_VSI;
dT  = R.LLHi_pos_table - R.LLHi_vis_table;
rel = {R.LLHrel_table, R.LLHrel_gauss, R.LLHrelI_gauss, R.LLHrel_pois, R.LLHrelI_pois};
dd  = {dT, dG, dGI, dP, dPI};
names = {'table (lasso, Gaussian)','Gaussian','Gaussian + V\timesS','Poisson','Poisson + V\timesS'};
plain = {'table','gauss','gauss+VxS','pois','pois+VxS'};

fid = fopen(fullfile(C.outRoot,'refit_summary.txt'),'w');
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('Refit of %d spatial cells from %d sessions\n', height(R), numel(unique(R.sessionid)));
out('corr with table: LLHi_vis %.3f, LLHi_pos %.3f\n', corr(R.LLHi_vis_table, R.LLHi_gauss_VS, 'rows','complete'), ...
    corr(R.LLHi_pos_table, R.LLHi_gauss_VSP, 'rows','complete'));
out('%-12s %14s %22s %26s\n', 'model', 'median LLHrel', 'median dLLH (bits/spk)', 'frac cells position helps');
for i = 1:5
    out('%-12s %14.3f %22.4f %26.2f\n', plain{i}, median(rel{i},'omitnan'), median(dd{i},'omitnan'), mean(dd{i} > 0,'omitnan'));
end
% how much of the position gain survives each change, per cell (relative to Gaussian)
k = dG > 0;
out('Among cells where position helps in the Gaussian refit (n=%d):\n', nnz(k));
out('   median fraction of position gain kept: +V x S %.2f | Poisson %.2f | Poisson+V x S %.2f\n', ...
    median(dGI(k)./dG(k)), median(dP(k)./dG(k)), median(dPI(k)./dG(k)));
out('Lambda (Poisson, chosen on VSP): %s\n', mat2str(arrayfun(@(v) nnz(R.lam_pois_VSP==v), [1e2 1e3 1e4])));
fclose(fid);

% figure
f = figure('Position',[50 50 1500 420],'Color','w');
tl = tiledlayout(f,1,4,'TileSpacing','compact');
nexttile; hold on;
for i = 1:5
    v = rel{i}; v = v(isfinite(v)); ci = bootci(500, @median, v);
    errorbar(i, median(v), median(v)-ci(1), ci(2)-median(v), 'ko', 'MarkerFaceColor', 'k');
end
yline(1,'b-'); set(gca,'XTick',1:5,'XTickLabel',names,'XTickLabelRotation',30); xlim([.5 5.5]);
ylabel('median LLH_{rel} (VSP / VS)'); box off; title('Position gain relative to vision+speed');
nexttile; hold on;
for i = 1:5, plot(i, mean(dd{i}>0,'omitnan'), 'ko', 'MarkerFaceColor','k'); end
yline(0.5,'k:'); set(gca,'XTick',1:5,'XTickLabel',names,'XTickLabelRotation',30); xlim([.5 5.5]); ylim([0 1]);
ylabel('fraction of cells with \DeltaLLH > 0'); box off; title('Cells where position improves CV fit');
nexttile; lim = prctile([dG; dP], [1 99]);
scatter(dG, dP, 6, 'k', 'filled', 'MarkerFaceAlpha', .25); hold on; plot(lim, lim, 'r-'); xline(0,'b:'); yline(0,'b:');
xlim(lim); ylim(lim); axis square; box off; xlabel('\DeltaLLH position, Gaussian'); ylabel('\DeltaLLH position, Poisson');
title('Per cell (bits/spike)');
nexttile; scatter(dT, dG, 6, 'k', 'filled', 'MarkerFaceAlpha', .25); hold on;
lim2 = prctile([dT; dG], [1 99]); plot(lim2, lim2, 'r-'); xline(0,'b:'); yline(0,'b:'); xlim(lim2); ylim(lim2); axis square; box off;
xlabel('\DeltaLLH position, table (lasso)'); ylabel('\DeltaLLH position, Gaussian ridge refit'); title('Original vs refit');
title(tl, sprintf('Steps 1-2: position gain under Gaussian/Poisson, with/without vision x speed (n = %d spatial cells)', height(R)));
exportgraphics(f, fullfile(C.figRoot, 'refit_position_gain.png'), 'Resolution', 150);
end
