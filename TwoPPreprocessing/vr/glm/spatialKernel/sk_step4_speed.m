function R = sk_step4_speed(mode)
% Step 4: does the residual "position kernel" scale with running speed?
%
% Base trials only. At every 2 cm bin, trials are split at the median running
% speed in that bin (fast / slow). For each cell:
%   r_fast = D_fast - V_fast,  r_slow = D_slow - V_slow   (V = vision+speed GLM prediction,
%                                                          which already includes speed additively)
%   K      = mean residual over all base trials
%   field  = bins where K >= 50% of its maximum
%   dObs   = mean_field(r_fast - r_slow) / mean_field(K)
% Additive position predicts dObs = 0. A gain on vision+speed (D = g.*V, g from
% base) predicts r_fast - r_slow = (g-1).*(V_fast - V_slow), i.e.
%   dGain  = mean_field((g-1).*(V_fast - V_slow)) / mean_field(K).
% A random (speed-independent) split gives the null distribution of dObs.
% mode: 'real' | 'simAdd' | 'simGain' (see sk_load_session).

if nargin < 1, mode = 'real'; end
rng(2);
C = sk_config();
T = readtable(C.tablePath);
nNull = 20;

out = {};
for iS = 1:size(C.sessions,1)
    S = sk_load_session(iS, T, mode);
    b = find(S.cond == 1);
    D = S.D(:,b,:); V = S.V(:,b,:); sp = S.speed(b,:);
    nC = size(D,1); nB = size(D,3);
    Db = squeeze(mean(D,2,'omitnan')); Vb = squeeze(mean(V,2,'omitnan'));
    [~, ~, gp] = sk_gain(Db, Vb);
    K = Db - Vb;
    field = K >= 0.5*max(K,[],2) & ~isnan(K);
    Kf = sum(K.*field,2,'omitnan')./sum(field,2);

    % median split per bin
    med = median(sp,1,'omitnan');
    H = sp > med; Lo = sp <= med; H(isnan(sp)) = false; Lo(isnan(sp)) = false;
    [dObs, dGain, vFast, vSlow] = splitStats(D, V, H, Lo, gp, field, Kf, sp);
    dNull = zeros(nC, nNull);
    for k = 1:nNull
        Hr = false(size(sp)); Lr = false(size(sp));
        for x = 1:nB
            v = find(~isnan(sp(:,x))); v = v(randperm(numel(v)));
            nh = nnz(H(:,x)); Hr(v(1:nh),x) = true; Lr(v(nh+1:end),x) = true;
        end
        dNull(:,k) = splitStats(D, V, Hr, Lr, gp, field, Kf, sp);
    end
    z = (dObs - mean(dNull,2)) ./ std(dNull,0,2);
    tb = table(S.T.CellIdx, repmat(S.sid,nC,1), S.T.cellNum, dObs, dGain, mean(dNull,2), std(dNull,0,2), z, ...
        repmat(vFast,nC,1), repmat(vSlow,nC,1), Kf, sum(field,2), S.T.spatialcells, S.T.goodcells, S.T.LLHrel, ...
        'VariableNames', {'CellIdx','sessionid','cellNum','dObs','dGainPred','dNullMean','dNullSD','zObs', ...
        'speedFast','speedSlow','kernelFieldMean','nFieldBins','spatialcells','goodcells','LLHrel'});
    out{end+1} = tb; %#ok<AGROW>
    s = S.T.spatialcells == 1;
    fprintf('%s: speed fast/slow %.1f/%.1f cm/s, median dObs %.3f, dGain %.3f (spatial)\n', ...
        S.tag, vFast, vSlow, median(dObs(s),'omitnan'), median(dGain(s),'omitnan'));
end
R = sortrows(vertcat(out{:}), 'CellIdx');
sfx = ''; if ~strcmp(mode,'real'), sfx = ['_' mode]; end
writetable(R, fullfile(C.outRoot, ['step4_speed_percell' sfx '.csv']));
save(fullfile(C.outRoot, ['step4_speed' sfx '.mat']), 'R');
plotStep4(R, C, sfx);
end

function [dObs, dGain, vFast, vSlow] = splitStats(D, V, H, Lo, gp, field, Kf, sp)
H3 = reshape(H, [1 size(H)]); L3 = reshape(Lo, [1 size(Lo)]);
mn = @(M, W) squeeze(sum(M.*W,2,'omitnan') ./ sum(W.*~isnan(M),2));
Dh = mn(D, H3); Dl = mn(D, L3); Vh = mn(V, H3); Vl = mn(V, L3);
dr = (Dh - Vh) - (Dl - Vl);
dObs = sum(dr.*field,2,'omitnan')./sum(field,2) ./ Kf;
gd = (gp(Vh) - Vh) - (gp(Vl) - Vl);      % gain-predicted residual difference
dGain = sum(gd.*field,2,'omitnan')./sum(field,2) ./ Kf;
vFast = mean(sp(H),'omitnan'); vSlow = mean(sp(Lo),'omitnan');
end

function plotStep4(R, C, sfx)
s = R.spatialcells == 1 & isfinite(R.dObs) & isfinite(R.dGainPred);
Rs = R(s,:);
f = figure('Position',[50 50 1400 420],'Color','w');
tl = tiledlayout(f,1,4,'TileSpacing','compact');
lim = [-1.5 1.5];
nexttile; histogram(max(min(Rs.dObs,2),-2), -2:0.1:2, 'FaceColor',[.4 .4 .4], 'Normalization','pdf'); hold on;
nullSamp = Rs.dNullMean + Rs.dNullSD.*randn(height(Rs),1);   % one draw per cell from its null
histogram(max(min(nullSamp,2),-2), -2:0.1:2, 'DisplayStyle','stairs','EdgeColor',[.5 .5 1], 'Normalization','pdf','LineWidth',1.2);
xline(0,'b-'); xline(median(Rs.dObs),'k--',sprintf('median %.2f',median(Rs.dObs)));
xlabel('(r_{fast} - r_{slow}) / K   in field'); ylabel('spatial cells'); legend({'speed split','random split'},'Box','off'); box off;
title('Kernel change with speed');
nexttile; c = Rs.dObs; p = Rs.dGainPred; k = abs(c)<3 & abs(p)<3;
scatter(p(k), c(k), 6, 'k', 'filled', 'MarkerFaceAlpha', .2); hold on;
plot(lim, lim, 'r-'); plot(lim, [0 0], 'b-'); xlim(lim); ylim(lim);
sl = robustfit(p(k), c(k)); plot(lim, sl(1)+sl(2)*lim, 'k--');
xlabel('gain-model prediction'); ylabel('observed'); box off; axis square;
title(sprintf('slope %.2f (0 additive, 1 gain)', sl(2)));
nexttile; q = discretize(Rs.LLHrel, quantile(Rs.LLHrel,[0 .25 .5 .75 1]));
mo = splitapply(@median, Rs.dObs, q); mp = splitapply(@median, Rs.dGainPred, q); mn = splitapply(@median, Rs.dNullMean, q);
plot(1:4, mo, 'ko-', 'MarkerFaceColor','k'); hold on; plot(1:4, mp, 'ro-'); plot(1:4, mn, 'o-', 'Color', [.6 .6 .6]);
yline(0,'b-'); xlim([.5 4.5]); set(gca,'XTick',1:4); xlabel('LLH_{rel} quartile'); ylabel('median'); box off;
legend({'observed','gain prediction','random split'},'Box','off','Location','best'); title('By position-model gain');
nexttile; histogram(max(min(Rs.zObs,6),-6), -6:0.25:6, 'FaceColor',[.4 .4 .4]); xline(0,'b-');
xlabel('z vs random split'); ylabel('spatial cells'); box off;
title(sprintf('%.0f%% cells z>2, %.0f%% z<-2', 100*mean(Rs.zObs>2), 100*mean(Rs.zObs<-2)));
title(tl, ['Step 4: does the position kernel scale with running speed? (base trials) ' strrep(sfx,'_',' ')]);
exportgraphics(f, fullfile(C.figRoot, ['step4_speed' sfx '.png']), 'Resolution', 150);
end
