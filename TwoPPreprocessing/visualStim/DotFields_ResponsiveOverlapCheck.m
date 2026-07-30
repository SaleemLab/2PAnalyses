% DotFields_ResponsiveOverlapCheck.m
%
% Checks the overlap (by ROI/bouton identity) between boutons classified
% as visually-responsive during stationary trials vs. running trials.
% Answers: are the 261 stationary-responsive boutons mostly a SUBSET of
% the 962 running-responsive boutons ("running preserves + adds"), or are
% a substantial fraction of them NOT running-responsive ("different
% subpopulations visible in each state")?
%
% Assumes allDotUnits already has .isResponsive_stat and
% .isResponsive_run fields, as produced by
% DotFields_TuningCurveAnalysis_compareStatesV2_2PData.m
% (SD-heuristic: ANOVA[speeds+blank]+threshold, p<0.05, >1SD above blank)

%%
respStat = [allDotUnits.isResponsive_stat]';   % logical, n=nBoutonsTotal
respRun  = [allDotUnits.isResponsive_run]';    % logical, n=nBoutonsTotal

nBoth       = sum(respStat & respRun);
nStatOnly   = sum(respStat & ~respRun);
nRunOnly    = sum(~respStat & respRun);
nNeither    = sum(~respStat & ~respRun);
nStatTotal  = sum(respStat);
nRunTotal   = sum(respRun);

fprintf('\n=== Responsive-bouton overlap (SD-heuristic, ANOVA+threshold) ===\n');
fprintf('Responsive in BOTH states:        %4d\n', nBoth);
fprintf('Responsive in STATIONARY only:    %4d\n', nStatOnly);
fprintf('Responsive in RUNNING only:       %4d\n', nRunOnly);
fprintf('Responsive in NEITHER:            %4d\n', nNeither);
fprintf('Total (check):                    %4d (should equal nBoutonsTotal = %d)\n', ...
    nBoth+nStatOnly+nRunOnly+nNeither, numel(respStat));

fprintf('\nOf the %d stationary-responsive boutons, %d (%.1f%%) are ALSO running-responsive.\n', ...
    nStatTotal, nBoth, 100*nBoth/nStatTotal);
fprintf('Of the %d running-responsive boutons, %d (%.1f%%) were ALSO stationary-responsive.\n', ...
    nRunTotal, nBoth, 100*nBoth/nRunTotal);
fprintf('%d boutons (%.1f%% of all running-responsive) are running-ONLY (i.e. "recruited" by locomotion).\n', ...
    nRunOnly, 100*nRunOnly/nRunTotal);

%% simple 2x2 visualization
figure('Color','w','Position',[200 200 380 380]);
countsMatrix = [nBoth, nStatOnly; nRunOnly, nNeither];
imagesc(countsMatrix);
colormap(gca, 'parula'); colorbar;
set(gca, 'XTick', 1:2, 'XTickLabel', {'Running: resp','Running: not resp'}, ...
         'YTick', 1:2, 'YTickLabel', {'Stationary: resp','Stationary: not resp'});
for r = 1:2
    for c = 1:2
        text(c, r, num2str(countsMatrix(r,c)), 'HorizontalAlignment','center', ...
            'FontSize', 13, 'FontWeight','bold', 'Color', 'w');
    end
end
title(sprintf('Responsive overlap (n=%d boutons)', numel(respStat)));

%% (optional) same check restricted to the ANOVA-protected t-test / Wilcoxon definitions,
% since these are the metrics you're actually leading with (254/255 stat, 930/947 run)
if isfield(allDotUnits, 'isResponsive_ttest_stat')
    respStat_tt = [allDotUnits.isResponsive_ttest_stat]';
    respRun_tt  = [allDotUnits.isResponsive_ttest_run]';
    nBoth_tt     = sum(respStat_tt & respRun_tt);
    nStatOnly_tt = sum(respStat_tt & ~respRun_tt);
    nRunOnly_tt  = sum(~respStat_tt & respRun_tt);

    fprintf('\n=== Same check, ANOVA-protected Welch t-test definition ===\n');
    fprintf('Both: %d | Stationary-only: %d | Running-only: %d\n', nBoth_tt, nStatOnly_tt, nRunOnly_tt);
    fprintf('Of %d stationary-responsive (t-test), %d (%.1f%%) also running-responsive.\n', ...
        sum(respStat_tt), nBoth_tt, 100*nBoth_tt/sum(respStat_tt));
end
