function R = sk_step3_swapomit(mode)
% mode: 'real' (default), or 'simAdd' / 'simGain' to calibrate the estimator
% on simulated data with known ground truth (see sk_load_session).
if nargin < 1, mode = 'real'; end
rng(1);
% Step 3: does the position kernel stay put in swap/omit trials (additive
% position signal) or follow the visual drive (gain on vision+speed)?
%
% For each cell, from BASE trials only:
%   K(x) = Dbase(x) - Vbase(x)                       additive position kernel
%   g(x) = (Dbase(x)+e) / (Vbase(x)+e)               multiplicative gain (sk_gain)
% where V is the stored GLM vision+speed prediction (clipped at 0) and e = 20% of mean rate.
% Then, in each manipulated condition m (V_m already reflects the changed cues):
%   additive prediction  Padd  = V_m + K
%   gain prediction      Pgain = (V_m+e).*g - e
% Population beta = sum(num)/sum(den) over cells (precision weighted), with a
% bootstrap over cells; per-cell beta is reported too but is heavy-tailed.
% evaluated only on bins where the two predictions differ (|Pgain-Padd| >= 25% of its max).
%   beta = projection of (D_m - Padd) onto (Pgain - Padd): 0 = additive, 1 = gain
% Single trials: residual r_t = D_t - V_t is compared with K and with the
% gain-predicted residual G_t = (V_t+e).*g - e - V_t; we record whether the
% residual peak lies closer to the K peak or to the G_t peak (trials where those
% two peaks are >= 10 cm apart).

C = sk_config();
T = readtable(C.tablePath);
if ~exist(C.figRoot,'dir'), mkdir(C.figRoot); end
minTrials = 3; frac = 0.25; minSep = 5;   % bins (10 cm)

out = {};
ex = struct('tag',{},'cell',{},'cond',{},'Db',{},'Vb',{},'Dm',{},'Vm',{},'Padd',{},'Pgain',{});
for iS = 1:size(C.sessions,1)
    S = sk_load_session(iS, T, mode);
    nC = size(S.D,1);
    b = S.cond == 1;
    Db = squeeze(mean(S.D(:,b,:),2,'omitnan'));
    Vb = squeeze(mean(S.V(:,b,:),2,'omitnan'));
    K  = Db - Vb;
    [~, ~, gp] = sk_gain(Db, Vb);

    num = zeros(nC,1); den = zeros(nC,1);
    sseA = zeros(nC,1); sseG = zeros(nC,1); nEval = zeros(nC,1);
    betaC = nan(nC,6);
    nCloserG = zeros(nC,1); nCloserA = zeros(nC,1);
    for m = 2:6
        tm = find(S.cond == m);
        if numel(tm) < minTrials, continue; end
        Dm = squeeze(mean(S.D(:,tm,:),2,'omitnan'));
        Vm = squeeze(mean(S.V(:,tm,:),2,'omitnan'));
        Padd  = Vm + K;
        Pgain = gp(Vm);
        dP = Pgain - Padd;
        msk = abs(dP) >= frac*max(abs(dP),[],2) & ~isnan(Dm) & ~isnan(dP);
        y = (Dm - Padd); y(~msk) = 0; dPm = dP; dPm(~msk) = 0;
        nm = sum(y.*dPm,2); dm = sum(dPm.^2,2);
        betaC(:,m) = nm ./ dm;
        num = num + nm; den = den + dm;
        rA = Dm - Padd; rG = Dm - Pgain; rA(~msk)=0; rG(~msk)=0;
        sseA = sseA + sum(rA.^2,2); sseG = sseG + sum(rG.^2,2); nEval = nEval + sum(msk,2);

        % keep a few examples (spatial cells with big kernels)
        if numel(ex) < 40
            cand = find(S.T.spatialcells==1 & dm > prctile(dm,90));
            for c = cand(1:min(2,end))'
                ex(end+1) = struct('tag',S.tag,'cell',c,'cond',m,'Db',Db(c,:),'Vb',Vb(c,:), ...
                    'Dm',Dm(c,:),'Vm',Vm(c,:),'Padd',Padd(c,:),'Pgain',Pgain(c,:)); %#ok<AGROW>
            end
        end

        % single-trial peak test
        for t = tm'
            Dt = squeeze(S.D(:,t,:)); Vt = squeeze(S.V(:,t,:));
            rt = Dt - Vt;
            Gt = gp(Vt) - Vt;
            valid = ~isnan(rt) & ~isnan(Gt) & ~isnan(K);
            rt(~valid) = -Inf; Kv = K; Kv(~valid) = -Inf; Gv = Gt; Gv(~valid) = -Inf;
            [~,pr] = max(rt,[],2); [~,pk] = max(Kv,[],2); [~,pg] = max(Gv,[],2);
            use = abs(pk - pg) >= minSep & any(valid,2);
            dA = abs(pr - pk); dG = abs(pr - pg);
            nCloserG = nCloserG + (use & dG < dA);
            nCloserA = nCloserA + (use & dA < dG);
        end
    end
    beta = num ./ den;
    tb = table(S.T.CellIdx, repmat(S.sid,nC,1), S.T.cellNum, beta, num, den, betaC(:,2), betaC(:,3), ...
        betaC(:,4), betaC(:,5), betaC(:,6), sseA./max(nEval,1), sseG./max(nEval,1), nEval, ...
        nCloserA, nCloserG, S.T.spatialcells, S.T.goodcells, S.T.LLHrel, ...
        'VariableNames', {'CellIdx','sessionid','cellNum','beta','betaNum','betaDen','beta_swap23','beta_swap34', ...
        'beta_omit2','beta_omit3','beta_omit4','mseAdd','mseGain','nEvalBins', ...
        'nTrialsPeakCloserAdd','nTrialsPeakCloserGain','spatialcells','goodcells','LLHrel'});
    out{end+1} = tb; %#ok<AGROW>
    fprintf('%s: median beta (spatial) %.2f\n', S.tag, median(beta(S.T.spatialcells==1),'omitnan'));
end
R = sortrows(vertcat(out{:}), 'CellIdx');
sfx = ''; if ~strcmp(mode,'real'), sfx = ['_' mode]; end
writetable(R, fullfile(C.outRoot, ['step3_swapomit_percell' sfx '.csv']));
save(fullfile(C.outRoot, ['step3_swapomit' sfx '.mat']), 'R', 'ex');
sk_plot_step3(R, ex, C, sfx);
end

function sk_plot_step3(R, ex, C, sfx)
sp = R.spatialcells == 1 & ~isnan(R.beta) & isfinite(R.beta);
Rs = R(sp,:);
f = figure('Position',[50 50 1500 800],'Color','w');
tl = tiledlayout(f,2,4,'TileSpacing','compact');

nexttile; b = max(min(Rs.beta, 2.5), -1.5);
histogram(b, -1.5:0.1:2.5, 'FaceColor',[.4 .4 .4]); hold on;
xline(0,'b-','additive','LineWidth',1.5); xline(1,'r-','gain','LineWidth',1.5);
xline(median(Rs.beta),'k--', sprintf('median %.2f', median(Rs.beta)));
popB = sum(Rs.betaNum)/sum(Rs.betaDen);
ciB = bootci(1000, @(n,d) sum(n)/sum(d), Rs.betaNum, Rs.betaDen);
xlabel('\beta per cell (clipped for display)'); ylabel('spatial cells'); box off;
title(sprintf('n = %d; population \\beta = %.2f [%.2f %.2f]', height(Rs), popB, ciB));

nexttile; names = {'swap23','swap34','omit2','omit3','omit4'};
B = Rs{:, strcat('beta_', names)}; B(~isfinite(B)) = NaN;
med = median(B,'omitnan'); ci = bootci(500, @(x) median(x,'omitnan'), B);
bar(med,'FaceColor',[.6 .6 .6]); hold on; errorbar(1:5, med, med-ci(1,:), ci(2,:)-med, 'k.', 'LineWidth',1);
yline(0,'b-'); yline(1,'r-'); set(gca,'XTickLabel',names); ylabel('median \beta (95% CI)'); box off;
title('\beta by manipulation');

nexttile; q = discretize(Rs.LLHrel, quantile(Rs.LLHrel,[0 .25 .5 .75 1]));
mq = splitapply(@(x) median(x,'omitnan'), Rs.beta, q);
ciq = zeros(2,4); for i=1:4, ciq(:,i) = bootci(500,@median,Rs.beta(q==i)); end
errorbar(1:4, mq, mq-ciq(1,:)', ciq(2,:)'-mq, 'ko-','MarkerFaceColor','k'); hold on;
yline(0,'b-'); yline(1,'r-'); xlim([.5 4.5]); set(gca,'XTick',1:4);
xlabel('LLH_{rel} quartile'); ylabel('median \beta'); box off; title('\beta vs position-model gain');

nexttile; nA = Rs.nTrialsPeakCloserAdd; nG = Rs.nTrialsPeakCloserGain;
pq = zeros(4,1); pci = zeros(4,2);
for i = 1:4
    a = sum(nA(q==i)); gg = sum(nG(q==i)); [pq(i), pci(i,:)] = binofit(gg, a+gg);
end
errorbar(1:4, pq, pq-pci(:,1), pci(:,2)-pq, 'ko-','MarkerFaceColor','k'); yline(0.5,'k:');
xlim([.5 4.5]); set(gca,'XTick',1:4); ylim([0 1]);
xlabel('LLH_{rel} quartile'); ylabel('P(trial peak closer to gain peak)'); box off;
tot = sum(nA)+sum(nG);
title(sprintf('Single-trial peaks: %.2f gain (n=%d trials\\timescells)', sum(nG)/tot, tot));

% examples
for k = 1:min(4, numel(ex))
    nexttile; E = ex(k); x = C.binCentres;
    plot(x,E.Db,'Color',[.6 .6 .6]); hold on; plot(x,E.Dm,'k','LineWidth',1.5);
    plot(x,E.Vm,'g'); plot(x,E.Padd,'b'); plot(x,E.Pgain,'r');
    box off; xlim([0 200]); xlabel('position (cm)'); ylabel('Hz');
    title(sprintf('%s c%d %s', strrep(E.tag,'_',' '), E.cell, C.condNames{E.cond}));
    if k == 1, legend({'data base','data manip','VS manip','VS+K','VS\timesg'},'Box','off','Location','best'); end
end
title(tl, ['Step 3: additive position kernel vs gain, tested on swap/omit trials (K, g from base trials only) ' strrep(sfx,'_',' ')]);
exportgraphics(f, fullfile(C.figRoot,['step3_swapomit' sfx '.png']), 'Resolution', 150);
end
