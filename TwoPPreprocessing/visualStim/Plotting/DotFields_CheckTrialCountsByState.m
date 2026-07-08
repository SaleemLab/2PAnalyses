% DotFields_CheckTrialCountsByState.m
%
% Reports trial counts per behavioral state (stationary vs running),
% per session, summed across speed conditions (plus blank separately).
% Trial counts are the same for every bouton WITHIN a given session
% (since state classification happens at the trial/session level, not
% per bouton), so this picks ONE representative bouton per unique
% session rather than redundantly recomputing the same numbers for
% every bouton in that session.

uniqueSessions = unique({allDotUnits.sessionLabel}, 'stable');
nSessions = numel(uniqueSessions);

nStatTrials_bySpeed = nan(nSessions, 1);
nRunTrials_bySpeed  = nan(nSessions, 1);
nStatBlank = nan(nSessions, 1);
nRunBlank  = nan(nSessions, 1);

fprintf('%-30s %12s %12s %12s %12s\n', 'Session', 'Stat(speed)', 'Run(speed)', 'Stat(blank)', 'Run(blank)');
for iSess = 1:nSessions
    thisSession = uniqueSessions{iSess};
    b = find(strcmp({allDotUnits.sessionLabel}, thisSession), 1); % first bouton in this session

    nStatTrials_bySpeed(iSess) = sum(cellfun(@numel, allDotUnits(b).alltraces(:,1)));
    nRunTrials_bySpeed(iSess)  = sum(cellfun(@numel, allDotUnits(b).alltraces(:,2)));
    nStatBlank(iSess) = numel(allDotUnits(b).blankTrials{1});
    nRunBlank(iSess)  = numel(allDotUnits(b).blankTrials{2});

    fprintf('%-30s %12d %12d %12d %12d\n', thisSession, ...
        nStatTrials_bySpeed(iSess), nRunTrials_bySpeed(iSess), nStatBlank(iSess), nRunBlank(iSess));
end

fprintf('\n=== Summary across %d sessions ===\n', nSessions);
fprintf('Stationary (speed trials):  median=%.0f, min=%.0f, max=%.0f\n', ...
    median(nStatTrials_bySpeed), min(nStatTrials_bySpeed), max(nStatTrials_bySpeed));
fprintf('Running (speed trials):     median=%.0f, min=%.0f, max=%.0f\n', ...
    median(nRunTrials_bySpeed), min(nRunTrials_bySpeed), max(nRunTrials_bySpeed));
fprintf('Stationary (blank trials):  median=%.0f, min=%.0f, max=%.0f\n', ...
    median(nStatBlank), min(nStatBlank), max(nStatBlank));
fprintf('Running (blank trials):     median=%.0f, min=%.0f, max=%.0f\n', ...
    median(nRunBlank), min(nRunBlank), max(nRunBlank));

nLowRun  = sum(nRunTrials_bySpeed < 10);
nLowStat = sum(nStatTrials_bySpeed < 10);
fprintf('\nSessions with <10 total running speed-trials: %d / %d\n', nLowRun, nSessions);
fprintf('Sessions with <10 total stationary speed-trials: %d / %d\n', nLowStat, nSessions);

%% ===================== plots =====================
figure('Position', [100 100 900 400]);

subplot(1,2,1);
bar([nStatTrials_bySpeed, nRunTrials_bySpeed]);
legend({'Stationary', 'Running'}, 'Location', 'best');
xlabel('Session index'); ylabel('Total trials (summed across speeds)');
title('Speed-condition trial counts by state, per session');

subplot(1,2,2);
scatter(nStatTrials_bySpeed, nRunTrials_bySpeed, 40, 'filled');
hold on;
maxVal = max([nStatTrials_bySpeed; nRunTrials_bySpeed]);
plot([0 maxVal], [0 maxVal], 'k--');
xlabel('Stationary trials'); ylabel('Running trials');
title('Stationary vs running trial count (each point = one session)');
axis equal; xlim([0 maxVal]); ylim([0 maxVal]);

%% ===================== per-speed-condition breakdown =====================
% Totals can hide a thin individual speed condition (e.g. a "total" of
% 20 running trials spread across 6 speeds is only ~3 trials/speed on
% average -- right at or below what cross-val R^2 (kval=3) needs).

minTrialWarn = 5; % flag any speed x state cell below this count

fprintf('\n=== Per-speed-condition trial counts (flagging < %d) ===\n', minTrialWarn);
for iSess = 1:nSessions
    thisSession = uniqueSessions{iSess};
    b = find(strcmp({allDotUnits.sessionLabel}, thisSession), 1);
    nSpeedsHere = size(allDotUnits(b).alltraces, 1);

    fprintf('\n--- %s ---\n', thisSession);
    fprintf('%-10s %12s %12s\n', 'Speed#', 'Stationary', 'Running');
    for s = 1:nSpeedsHere
        nStat = numel(allDotUnits(b).alltraces{s,1});
        nRun  = numel(allDotUnits(b).alltraces{s,2});
        statFlag = ''; runFlag = '';
        if nStat < minTrialWarn, statFlag = ' <-- LOW'; end
        if nRun  < minTrialWarn, runFlag  = ' <-- LOW'; end
        fprintf('%-10d %12d%-8s %12d%-8s\n', s, nStat, statFlag, nRun, runFlag);
    end
end