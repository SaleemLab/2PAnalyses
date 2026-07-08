% DotFields_CheckTrialCountVsResponsiveFraction.m
%
% Tests whether trial count is driving the stationary-vs-running
% responsive-fraction gap, by checking whether SESSIONS with more
% trials (in a given state) also show a HIGHER responsive fraction in
% that same state. If trial count were the real driver, this
% correlation should be clearly positive. If it's a genuine
% locomotion-gain effect (not a counting artifact), this correlation
% should be weak/absent -- the effect should show up broadly regardless
% of how many trials a given session happened to have.
%
% NOTE: only 6 sessions -- this correlation will necessarily be noisy
% and low-powered on its own, but it's still informative as a sanity
% check, especially combined with the qualitative pattern already seen
% (stationary generally has MORE trials but FEWER responsive boutons,
% which already argues against a trial-count explanation).

uniqueSessions = unique({allDotUnits.sessionLabel}, 'stable');
nSessions = numel(uniqueSessions);

nTrials_stat = nan(nSessions, 1);
nTrials_run  = nan(nSessions, 1);
fracResp_stat = nan(nSessions, 1);
fracResp_run  = nan(nSessions, 1);

for iSess = 1:nSessions
    thisSession = uniqueSessions{iSess};
    sessMask = strcmp({allDotUnits.sessionLabel}, thisSession);
    b = find(sessMask, 1);

    nTrials_stat(iSess) = sum(cellfun(@numel, allDotUnits(b).alltraces(:,1)));
    nTrials_run(iSess)  = sum(cellfun(@numel, allDotUnits(b).alltraces(:,2)));

    fracResp_stat(iSess) = mean([allDotUnits(sessMask).isResponsive_stat]);
    fracResp_run(iSess)  = mean([allDotUnits(sessMask).isResponsive_run]);
end

fprintf('%-30s %12s %12s %14s %14s\n', 'Session', 'nTrials(stat)', 'nTrials(run)', 'fracResp(stat)', 'fracResp(run)');
for iSess = 1:nSessions
    fprintf('%-30s %12d %12d %14.3f %14.3f\n', uniqueSessions{iSess}, ...
        nTrials_stat(iSess), nTrials_run(iSess), fracResp_stat(iSess), fracResp_run(iSess));
end

[r_stat, p_stat] = corr(nTrials_stat, fracResp_stat);
[r_run, p_run]   = corr(nTrials_run, fracResp_run);

fprintf('\n=== Correlation: trial count vs responsive fraction (n=%d sessions) ===\n', nSessions);
fprintf('Stationary: r=%.3f, p=%.3f\n', r_stat, p_stat);
fprintf('Running:    r=%.3f, p=%.3f\n', r_run, p_run);
fprintf('\n(A strong positive r would suggest trial count is driving responsive rate.\n');
fprintf('A weak/near-zero r argues AGAINST a trial-count explanation.)\n');

%% ===================== plot =====================
figure('Position', [100 100 900 400]);

subplot(1,2,1);
scatter(nTrials_stat, fracResp_stat, 60, 'filled', 'MarkerFaceColor', 'k');
lsline;
xlabel('Number of stationary trials (this session)'); ylabel('Fraction responsive (stationary, ANOVA+2SD)');
title(sprintf('Stationary: r=%.2f, p=%.2f', r_stat, p_stat));

subplot(1,2,2);
scatter(nTrials_run, fracResp_run, 60, 'filled', 'MarkerFaceColor', 'r');
lsline;
xlabel('Number of running trials (this session)'); ylabel('Fraction responsive (running, ANOVA+2SD)');
title(sprintf('Running: r=%.2f, p=%.2f', r_run, p_run));

sgtitle('Does trial count predict responsive fraction, per session?');