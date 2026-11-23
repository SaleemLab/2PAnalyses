% Run steps 3 and 4 on real data and on the two simulated ground truths,
% and write a one-page text summary to <outRoot>\steps34_summary.txt.
C = sk_config();
fid = fopen(fullfile(C.outRoot,'steps34_summary.txt'),'w');
p = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
modes = {'real','simAdd','simGain'};
for i = 1:3
    R = sk_step3_swapomit(modes{i});
    s = R.spatialcells == 1 & isfinite(R.beta);
    popB = sum(R.betaNum(s))/sum(R.betaDen(s));
    ci = bootci(1000, @(n,d) sum(n)/sum(d), R.betaNum(s), R.betaDen(s));
    a = sum(R.nTrialsPeakCloserAdd(s)); g = sum(R.nTrialsPeakCloserGain(s));
    p('STEP3 %-8s n=%d  population beta %.2f [%.2f %.2f]  median cell beta %.2f  additive wins %.2f  trial peak nearer gain %.3f\n', ...
        modes{i}, nnz(s), popB, ci, median(R.beta(s)), mean(R.mseAdd(s) < R.mseGain(s)), g/(a+g));
    q = discretize(R.LLHrel(s), quantile(R.LLHrel(s), [0 .25 .5 .75 1]));
    nb = R.betaNum(s); db = R.betaDen(s);
    p('      population beta by LLHrel quartile: %s\n', mat2str(arrayfun(@(k) sum(nb(q==k))/sum(db(q==k)), 1:4), 2));
    sids = unique(R.sessionid(s));
    p('      population beta by session: %s\n', mat2str(arrayfun(@(k) sum(nb(R.sessionid(s)==k))/sum(db(R.sessionid(s)==k)), sids'), 2));
end
for i = 1:3
    R = sk_step4_speed(modes{i});
    s = R.spatialcells == 1 & isfinite(R.dObs) & isfinite(R.dGainPred);
    k = s & abs(R.dObs) < 3 & abs(R.dGainPred) < 3;
    b = robustfit(R.dGainPred(k), R.dObs(k));
    p('STEP4 %-8s n=%d  median dObs %.3f  median gain-pred %.3f  median null %.3f  slope obs~gainpred %.2f  z>2 %.2f  z<-2 %.2f\n', ...
        modes{i}, nnz(s), median(R.dObs(s)), median(R.dGainPred(s)), median(R.dNullMean(s)), b(2), mean(R.zObs(s)>2), mean(R.zObs(s)<-2));
end
fclose(fid);
