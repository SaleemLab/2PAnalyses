% Session-level (and mouse-level) summary: stationary vs. locomotion
%
% Treats SESSION (n=6) -- and optionally MOUSE (n=3, more conservative)
% -- as the unit of replication, instead of pooling thousands of
% non-independent boutons. This avoids the pseudoreplication problem in
% the bouton-level panels (Figure4_5_DotFields_TuningCurves_behSplit.m),
% at the cost of much lower statistical power -- expected and honest for
% a dataset this size. Statistics here (sign test, signrank) should be
% read as "how consistent is the direction of the effect across
% sessions/mice", not as a strong significance claim.
%
% Requires: allDotUnits, statR2, runR2, isTunedStat, isTunedRun,
% validIdx_both already in workspace (run
% DotFields_IncTuningCurveAnalysis_compareStatesV2_2PData first).

%% Build per-session summary table
sessionLabels_all = {allDotUnits.sessionLabel}';
uniqueSessionsAll  = unique(sessionLabels_all, 'stable');
nSessions = numel(uniqueSessionsAll);

Session = strings(nSessions,1);
PctTunedStat = nan(nSessions,1); PctTunedRun = nan(nSessions,1);
MeanR2Stat   = nan(nSessions,1); MeanR2Run   = nan(nSessions,1);
PctPrefSpeedMatch = nan(nSessions,1); MeanAbsPrefSpeedDiff = nan(nSessions,1);
PctShapeMatch = nan(nSessions,1);
NValidPair = nan(nSessions,1); NDualTuned = nan(nSessions,1);

for s = 1:nSessions
    thisSession = uniqueSessionsAll{s};
    sessIdx = strcmp(sessionLabels_all, thisSession);

    % --- Prevalence: % tuned, among boutons with a valid R^2 in both states ---
    validPair = sessIdx(:) & ~isnan(statR2) & ~isnan(runR2);
    nValid = sum(validPair);
    NValidPair(s) = nValid;
    if nValid > 0
        PctTunedStat(s) = 100 * sum(isTunedStat(validPair)) / nValid;
        PctTunedRun(s)  = 100 * sum(isTunedRun(validPair))  / nValid;

        % --- R^2 magnitude, among boutons tuned in that state ---
        tunedStatIdx = validPair & isTunedStat;
        tunedRunIdx  = validPair & isTunedRun;
        if any(tunedStatIdx), MeanR2Stat(s) = mean(statR2(tunedStatIdx)); end
        if any(tunedRunIdx),  MeanR2Run(s)  = mean(runR2(tunedRunIdx));  end
    end

    % --- Preferred speed & shape agreement, dual-tuned boutons only ---
    dualIdxAll = find(sessIdx(:));
    dualIdx = intersect(dualIdxAll, validIdx_both);
    NDualTuned(s) = numel(dualIdx);
    if ~isempty(dualIdx)
        prefStat = cat(1, allDotUnits(dualIdx).prefSpeed_stat);
        prefRun  = cat(1, allDotUnits(dualIdx).prefSpeed_run);
        PctPrefSpeedMatch(s)     = 100 * mean(prefStat == prefRun);
        MeanAbsPrefSpeedDiff(s)  = mean(abs(prefStat - prefRun));

        shapeStat = cat(1, allDotUnits(dualIdx).gaussChar_stat);
        shapeRun  = cat(1, allDotUnits(dualIdx).gaussChar_run);
        validShape = ~isnan(shapeStat) & ~isnan(shapeRun);
        if any(validShape)
            PctShapeMatch(s) = 100 * mean(shapeStat(validShape) == shapeRun(validShape));
        end
    end

    Session(s) = thisSession;
end

sessionTable = table(Session, NValidPair, PctTunedStat, PctTunedRun, MeanR2Stat, MeanR2Run, ...
    NDualTuned, PctPrefSpeedMatch, MeanAbsPrefSpeedDiff, PctShapeMatch);
disp(sessionTable);

%% Paired dot-and-line plots + sign test / signrank, per question
% Sign test is the more honest test at n=6 -- signrank's own minimum
% possible p-value with 6 pairs is ~0.03 even with perfect agreement, so
% don't over-read small differences between the two.

metricPairs = {
    'PctTunedStat', 'PctTunedRun', 'Tuning prevalence (%)',     'Prevalence';
    'MeanR2Stat',   'MeanR2Run',   'Mean R^2 (tuned boutons)',  'R2 magnitude';
    };

for m = 1:size(metricPairs,1)
    colStat = metricPairs{m,1};
    colRun  = metricPairs{m,2};
    yLabelStr  = metricPairs{m,3};
    metricName = metricPairs{m,4};

    valStat = sessionTable.(colStat);
    valRun  = sessionTable.(colRun);
    keepIdx = ~isnan(valStat) & ~isnan(valRun);
    valStat = valStat(keepIdx); valRun = valRun(keepIdx);
    nSessM  = numel(valStat);

    pSign = signtest(valStat, valRun);
    pSR   = signrank(valStat, valRun);

    figure('Color', 'w', 'Name', sprintf('Session-level: %s', metricName), 'Position', [200, 200, 350, 420]);
    hold on;
    for i = 1:nSessM
        plot([1 2], [valStat(i) valRun(i)], '-o', 'Color', [0.6 0.6 0.6], ...
            'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerSize', 6);
    end
    xlim([0.5 2.5]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Stationary', 'Locomotion'});
    ylabel(yLabelStr);
    title(sprintf('%s, per session (n=%d)\nSign test p=%.3g | signrank p=%.3g', ...
        metricName, nSessM, pSign, pSR), 'FontSize', 10);
    box off; hold off;
end

%% Preferred speed & shape agreement, per session (single value, not paired)
% These aren't stat-vs-run pairs -- they're "how much do the two states
% agree" per session, so shown as a single dot per session against a
% chance-level reference line rather than a two-armed spaghetti plot.

figure('Color', 'w', 'Name', 'Agreement metrics, per session', 'Position', [200, 200, 700, 400]);

subplot(1,2,1);
scatter(1:nSessions, sessionTable.PctPrefSpeedMatch, 60, [0.2 0.4 0.7], 'filled');
yline(100/6, '--', 'Color', [0.6 0.6 0.6], 'Label', 'chance (1/6)');
xlim([0.5, nSessions+0.5]); ylim([0 100]);
set(gca, 'XTick', 1:nSessions, 'XTickLabel', sessionTable.Session, 'XTickLabelRotation', 45);
ylabel('% exact preferred-speed match'); title('Preferred speed agreement');
box off;

subplot(1,2,2);
scatter(1:nSessions, sessionTable.PctShapeMatch, 60, [0.7 0.3 0.3], 'filled');
yline(100/4, '--', 'Color', [0.6 0.6 0.6], 'Label', 'chance (1/4)');
xlim([0.5, nSessions+0.5]); ylim([0 100]);
set(gca, 'XTick', 1:nSessions, 'XTickLabel', sessionTable.Session, 'XTickLabelRotation', 45);
ylabel('% same tuning shape'); title('Tuning shape agreement');
box off;

%% Optional: same summary collapsed to MOUSE level (n=3), more conservative
mouseIDs_all = {allDotUnits.mouseID}';
uniqueMice = unique(mouseIDs_all, 'stable');
nMice = numel(uniqueMice);

MouseID = strings(nMice,1);
PctTunedStat_m = nan(nMice,1); PctTunedRun_m = nan(nMice,1);
MeanR2Stat_m   = nan(nMice,1); MeanR2Run_m   = nan(nMice,1);

for m = 1:nMice
    thisMouse = uniqueMice{m};
    mouseIdx = strcmp(mouseIDs_all, thisMouse);

    validPair = mouseIdx(:) & ~isnan(statR2) & ~isnan(runR2);
    if any(validPair)
        PctTunedStat_m(m) = 100 * sum(isTunedStat(validPair)) / sum(validPair);
        PctTunedRun_m(m)  = 100 * sum(isTunedRun(validPair))  / sum(validPair);
        tunedStatIdx = validPair & isTunedStat;
        tunedRunIdx  = validPair & isTunedRun;
        if any(tunedStatIdx), MeanR2Stat_m(m) = mean(statR2(tunedStatIdx)); end
        if any(tunedRunIdx),  MeanR2Run_m(m)  = mean(runR2(tunedRunIdx));  end
    end
    MouseID(m) = thisMouse;
end

mouseTable = table(MouseID, PctTunedStat_m, PctTunedRun_m, MeanR2Stat_m, MeanR2Run_m);
fprintf('\n--- Mouse-level summary (n=%d mice, most conservative) ---\n', nMice);
disp(mouseTable);
