% DotFields_TrialCountMatchedComparison.m
%
% Directly tests whether the running > stationary responsive-fraction
% gap is a trial-count artifact: for each bouton, downsamples BOTH
% states to the SAME (smaller) number of trials per speed condition
% (and per blank), then recomputes the SD-heuristic and ANOVA-protected
% t-test criteria on this trial-count-MATCHED data. If the gap persists
% after matching, trial count isn't the main driver. If it shrinks a
% lot, trial count was contributing more than the earlier correlational
% checks suggested.
%
% Uses REPEATED random subsampling (nRepeats) averaged together, since a
% single random subsample could be lucky/unlucky, especially at small n.
%
% Requires allDotUnits with alltraces, blankTrials (per state).

nRepeats = 20; % repeated random subsamples per bouton, averaged for stability
ALPHA = 0.05;
NSD = 2;
rng(1); % reproducible

nBoutonsTotal = numel(allDotUnits);

isResp_matched_stat  = false(nBoutonsTotal, nRepeats);
isResp_matched_run   = false(nBoutonsTotal, nRepeats);
isRespT_matched_stat = false(nBoutonsTotal, nRepeats);
isRespT_matched_run  = false(nBoutonsTotal, nRepeats);

for b = 1:nBoutonsTotal
    alltraces   = allDotUnits(b).alltraces;   % nSpeeds x 2
    blankTrials = allDotUnits(b).blankTrials; % 1 x 2
    nSpeeds = size(alltraces, 1);

    % per-speed matched trial count = min across the two states
    nMatchedPerSpeed = nan(nSpeeds, 1);
    for s = 1:nSpeeds
        n1 = sum(~isnan(alltraces{s,1}));
        n2 = sum(~isnan(alltraces{s,2}));
        nMatchedPerSpeed(s) = min(n1, n2);
    end
    nMatchedBlank = min(sum(~isnan(blankTrials{1})), sum(~isnan(blankTrials{2})));

    if any(nMatchedPerSpeed < 2) || nMatchedBlank < 2
        continue; % can't do a meaningful matched comparison for this bouton
    end

    for rep = 1:nRepeats
        for si = 1:2
            % --- build matched-trial-count data for this state ---
            matchedTraces = cell(nSpeeds, 1);
            for s = 1:nSpeeds
                vals = alltraces{s, si}(~isnan(alltraces{s, si}));
                sampIdx = randperm(numel(vals), nMatchedPerSpeed(s));
                matchedTraces{s} = vals(sampIdx);
            end
            blankVals = blankTrials{si}(~isnan(blankTrials{si}));
            matchedBlank = blankVals(randperm(numel(blankVals), nMatchedBlank));

            % --- ANOVA(speeds+blank) + 2SD threshold, on matched data ---
            y = []; grp = [];
            for s = 1:nSpeeds
                y = [y; matchedTraces{s}(:)];
                grp = [grp; repmat(s, numel(matchedTraces{s}), 1)];
            end
            y = [y; matchedBlank(:)];
            grp = [grp; repmat(nSpeeds+1, numel(matchedBlank), 1)];
            anovaP_m = anova1(y, grp, 'off');

            blankMean_m = mean(matchedBlank);
            blankSD_m   = std(matchedBlank);
            tuningVals = cellfun(@mean, matchedTraces);
            [~, prefIdx_m] = max(tuningVals);
            meanPref_m = mean(matchedTraces{prefIdx_m});

            isResponsive_m = (anovaP_m < ALPHA) && (meanPref_m > blankMean_m + NSD*blankSD_m);

            % --- ANOVA-protected t-test, on matched data ---
            [~, ttestP_m] = ttest2(matchedTraces{prefIdx_m}, matchedBlank, 'Vartype', 'unequal');
            isRespT_m = (anovaP_m < ALPHA) && (ttestP_m < ALPHA) && (meanPref_m > blankMean_m);

            if si == 1
                isResp_matched_stat(b, rep)  = isResponsive_m;
                isRespT_matched_stat(b, rep) = isRespT_m;
            else
                isResp_matched_run(b, rep)  = isResponsive_m;
                isRespT_matched_run(b, rep) = isRespT_m;
            end
        end
    end

    if mod(b, 200) == 0
        fprintf('Matched comparison: %d / %d boutons...\n', b, nBoutonsTotal);
    end
end

% average across repeats -> "probability of being responsive" per bouton, then threshold at 0.5 for a final call
fracResp_matched_stat  = mean(isResp_matched_stat, 2);
fracResp_matched_run   = mean(isResp_matched_run, 2);
fracRespT_matched_stat = mean(isRespT_matched_stat, 2);
fracRespT_matched_run  = mean(isRespT_matched_run, 2);

isResp_matched_stat_final  = fracResp_matched_stat  >= 0.5;
isResp_matched_run_final   = fracResp_matched_run   >= 0.5;
isRespT_matched_stat_final = fracRespT_matched_stat >= 0.5;
isRespT_matched_run_final  = fracRespT_matched_run  >= 0.5;

fprintf('\n=== TRIAL-COUNT-MATCHED comparison (n=%d boutons, %d repeats averaged) ===\n', nBoutonsTotal, nRepeats);
fprintf('%-40s %10s %10s\n', 'Method', 'Stationary', 'Running');
fprintf('%-40s %10d %10d\n', 'SD-heuristic (matched)', sum(isResp_matched_stat_final), sum(isResp_matched_run_final));
fprintf('%-40s %10d %10d\n', 'ANOVA-protected t-test (matched)', sum(isRespT_matched_stat_final), sum(isRespT_matched_run_final));

fprintf('\n=== For comparison, UNMATCHED (original, full-trial) counts ===\n');
fprintf('%-40s %10d %10d\n', 'SD-heuristic (original)', sum([allDotUnits.isResponsive_stat]), sum([allDotUnits.isResponsive_run]));
fprintf('%-40s %10d %10d\n', 'ANOVA-protected t-test (original)', sum([allDotUnits.isResponsive_ttest_stat]), sum([allDotUnits.isResponsive_ttest_run]));

%% ===================== plot: matched vs unmatched gap =====================
figure('Position', [100 100 700 400]);
categories = {'SD-heur (unmatched)', 'SD-heur (matched)', 't-test (unmatched)', 't-test (matched)'};
statCounts = [sum([allDotUnits.isResponsive_stat]), sum(isResp_matched_stat_final), ...
              sum([allDotUnits.isResponsive_ttest_stat]), sum(isRespT_matched_stat_final)];
runCounts  = [sum([allDotUnits.isResponsive_run]), sum(isResp_matched_run_final), ...
              sum([allDotUnits.isResponsive_ttest_run]), sum(isRespT_matched_run_final)];
bar([statCounts; runCounts]');
set(gca, 'XTickLabel', categories);
legend({'Stationary', 'Running'}, 'Location', 'best');
ylabel('Number of responsive boutons');
title('Does the running/stationary gap survive trial-count matching?');
